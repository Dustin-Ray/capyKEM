use super::{field_element::FieldElement as F, ring_element::RingElement};
use crate::constants::{ml_kem_constants::Q, K_MOD_ROOTS, K_NTT_ROOTS};
use core::fmt;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use core::ops::{AddAssign, Mul};
#[derive(Clone, Copy)]
pub struct NttElement {
    ring: [F; 256],
}

impl NttElement {
    fn new(r: &mut RingElement) -> Self {
        let mut ntt_el = NttElement {
            ring: r.coefficients,
        };
        ntt_el.ntt();
        ntt_el
    }

    pub fn zero() -> Self {
        NttElement {
            ring: [F::default(); 256],
        }
    }

    pub fn byte_encode_12() {}

    /// byte stream b should be rho||i||j
    /// per algorithm 12
    pub fn sample(byte_stream: Vec<u8>) -> Self {
        // Get the XOF for the input
        let mut prf = Shake256::default();
        prf.update(&byte_stream);

        // why is this 486
        let mut b = [0_u8; 486];
        let mut reader = prf.finalize_xof();
        reader.read(&mut b);

        let mut i = 0;
        let mut j = 0;
        let mut a_hat: [F; 256] = [F::default(); 256];

        // Populate the matrix â
        while j < 256 {
            let d_1 = u16::from_le_bytes([b[i], b[i + 1]]) & 0b1111_1111_1111; // Masking to get 12 bits
            let d_2 = u16::from_le_bytes([b[i + 1], b[i + 2]]) >> 4; // Shift right to get the next 12 bits

            if d_1 < Q {
                a_hat[j] = F::new(d_1);
                j += 1;
            }
            if d_2 < Q && j < 256 {
                a_hat[j] = F::new(d_2);
                j += 1;
            }
            i += 3
        }

        // return t_hat
        NttElement::new(&mut RingElement::new(a_hat))
    }

    #[inline(always)]
    fn multiply_ntts(&self, other: Self) -> Self {
        let mut h_hat = NttElement::zero();

        // Iterate over `K_MOD_ROOTS` with their indices
        for (i, &k_mod_root) in K_MOD_ROOTS.iter().enumerate() {
            (h_hat.ring[2 * i], h_hat.ring[2 * i + 1]) = NttElement::base_case_multiply(
                self.ring[2 * i],
                self.ring[(2 * i) + 1],
                other.ring[2 * i],
                other.ring[(2 * i) + 1],
                k_mod_root,
            )
        }

        h_hat
    }

    fn base_case_multiply(a_0: F, a_1: F, b_0: F, b_1: F, gamma: u16) -> (F, F) {
        let c_0 = (a_0 * b_0) + (a_1 * b_1) * gamma;
        let c_1 = (a_0 * b_1) + (a_1 * b_0);
        (c_0, c_1)
    }

    #[inline(always)]
    fn ntt(&mut self) {
        let mut k = 1;
        let mut len = 128;
        while len >= 2 {
            for start in (0..256).step_by(2 * len) {
                let zeta = K_NTT_ROOTS[k];
                k += 1;

                for j in start..start + len {
                    let t = zeta * self.ring[j + len] % Q;
                    self.ring[j + len] = self.ring[j] - F::new(t);
                    self.ring[j] += F::new(t);
                }
            }
            len /= 2;
        }
    }

    #[inline(always)]
    fn ntt_inv(&mut self) -> RingElement {
        let mut k = 127;
        let mut len = 2;
        while len <= 128 {
            for start in (0..256).step_by(2 * len) {
                let zeta = K_NTT_ROOTS[k];
                k -= 1;

                for j in start..start + len {
                    let t = self.ring[j];
                    self.ring[j] = t + self.ring[j + len];
                    self.ring[j + len] = F::new(zeta * (self.ring[j + len] - t));
                }
            }
            len *= 2;
        }
        for item in self.ring.iter_mut() {
            *item = *item * 3303;
        }
        RingElement::new(self.ring)
    }
}

impl AddAssign for NttElement {
    fn add_assign(&mut self, other: Self) {
        for (lhs, rhs) in self.ring.iter_mut().zip(other.ring.iter()) {
            *lhs += *rhs;
        }
    }
}

// Implementing From<[FieldElement; 256]> for RingElement
impl From<RingElement> for NttElement {
    fn from(mut val: RingElement) -> Self {
        NttElement::new(&mut val)
    }
}

// Implementing From<[FieldElement; 256]> for RingElement
impl From<&mut NttElement> for RingElement {
    fn from(val: &mut NttElement) -> Self {
        val.ntt_inv()
    }
}

impl PartialEq for NttElement {
    fn eq(&self, other: &Self) -> bool {
        self.ring.iter().zip(other.ring.iter()).all(|(a, b)| a == b)
    }
}

impl Mul<NttElement> for NttElement {
    type Output = Self;

    fn mul(self, rhs: NttElement) -> Self::Output {
        self.multiply_ntts(rhs)
    }
}

impl fmt::Debug for NttElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (index, element) in self.ring.iter().enumerate() {
            write!(f, "{:<8}", element.val())?;
            // Adjust for row width
            if (index + 1) % 16 == 0 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    use super::*;

    // REMARKS:
    // axiom tests:
    // -[ ] commutative
    // -[x] zero identity
    // -[x] associativity
    // -[ ] distributivity
    // -[x] multiplication closure

    #[test]
    fn test_sample_ntt() {
        let byte_stream: Vec<u8> = vec![42_u8; 32];
        let ntt_element = NttElement::sample(byte_stream);
        println!(
            "NTT element sampled from byte stream seed: \n\n{:?}",
            ntt_element
        );
    }

    #[test]
    fn test_ntt() {
        // sample output is in NTT domain
        let mut byte_stream = NttElement::sample(vec![42_u8; 32]);
        let mut byte_stream_copy = byte_stream;
        byte_stream.ntt_inv();
        byte_stream_copy.ntt_inv();
        assert_eq!(byte_stream_copy.ring, byte_stream.ring)
    }

    #[test]
    fn test_ntt_from_poly_cbd_inverse_with_random_input() {
        // Generate a random byte stream using a seeded RNG for reproducibility
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();
        // Sample a ring element using the random byte stream
        let mut ring_element = RingElement::sample_poly_cbd(&bytes, 0xFF);
        let ring_element_copy = ring_element;

        // runs .ntt() on intstantiation
        let mut ntt_element = NttElement::new(&mut ring_element);
        ntt_element.ntt_inv();
        assert_eq!(ring_element_copy.coefficients, ntt_element.ring);
    }

    #[test]
    fn test_multiply_ntts_associative() {
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        let a = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
        let b = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xBB));
        let c = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xCC));

        // Test associativity (ab)c = a(bc)
        let ab_c = (a * b) * c;
        let a_bc = a * (b * c);
        assert_eq!(ab_c.ring, a_bc.ring);

        let a = NttElement::sample(bytes.clone());
        let b = NttElement::sample(bytes.clone());
        let c = NttElement::sample(bytes.clone());

        // Test associativity (ab)c = a(bc)
        let ab_c = (a * b) * c;
        let a_bc = a * (b * c);
        assert_eq!(ab_c.ring, a_bc.ring);
    }

    #[test]
    fn test_multiply_ntts_zero() {
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        let a = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
        let zero = NttElement::zero();

        // Test multiplicative identity
        let res = zero * a;
        assert_eq!(
            res,
            NttElement {
                ring: [F::zero(); 256]
            }
        );
        let res = a * zero;
        assert_eq!(
            res,
            NttElement {
                ring: [F::zero(); 256]
            }
        );
    }

    #[test]
    fn test_closure_under_multiplication() {
        for _ in 0..1000 {
            let bytes: Vec<u8> = (0..32)
                .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
                .collect();

            let a = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
            let b = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xBB));

            let result = a * b;
            assert!(
                result.ring.iter().all(|x| x.val() < 3329),
                "Result of multiplication must be valid NttElement"
            );
        }
    }
}
