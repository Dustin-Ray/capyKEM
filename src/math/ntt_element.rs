use super::{field_element::FieldElement as F, ring_element::RingElement};
use crate::constants::ml_kem_constants::{N, Q};
use crate::constants::parameter_sets::ParameterSet;
use crate::constants::{K_MOD_ROOTS, K_NTT_ROOTS};
use core::fmt;
use core::ops::{AddAssign, Mul};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};

#[derive(Clone, Copy)]
pub struct NttElement<P> {
    pub ring: [F<P>; 256],
}

impl<P: ParameterSet + Copy> NttElement<P> {
    pub fn new(r: &mut RingElement<P>) -> Self {
        let mut ntt_el = NttElement {
            ring: r.coefficients,
        };
        ntt_el.ntt();
        ntt_el
    }

    pub fn zero() -> Self {
        NttElement {
            ring: [F::zero(); 256],
        }
    }

    pub fn get_ring(&self) -> [F<P>; 256] {
        self.ring
    }

    pub fn sample_ntt(rho: &[u8], ii: usize, jj: usize) -> NttElement<P> {
        let mut hasher = Shake128::default();
        hasher.update(rho);
        hasher.update(&[ii.try_into().unwrap(), jj.try_into().unwrap()]);

        let mut reader = hasher.finalize_xof();

        let mut a = NttElement::zero();
        let mut j = 0usize;
        let mut buf = [0u8; 24];
        let mut off = 24usize;
        while j < N.into() {
            if off >= 24 {
                reader.read(&mut buf);
                off = 0;
            }

            let d1 = u16::from_le_bytes([buf[off], buf[off + 1]]) & 0x0FFF;
            let d2 = ((u16::from(buf[off + 1]) | (u16::from(buf[off + 2]) << 8)) >> 4) & 0x0FFF;

            off += 3;

            if d1 < Q {
                a.ring[j] = F::new(d1);
                j += 1;
            }
            if j >= N.into() {
                break;
            }

            if d2 < Q {
                a.ring[j] = F::new(d2);
                j += 1;
            }
        }
        a
    }

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
            );
        }

        h_hat
    }

    fn base_case_multiply(a_0: F<P>, a_1: F<P>, b_0: F<P>, b_1: F<P>, gamma: u16) -> (F<P>, F<P>) {
        let c_0 = (a_0 * b_0) + (a_1 * b_1) * gamma;
        let c_1 = (a_0 * b_1) + (a_1 * b_0);
        (c_0, c_1)
    }

    // This should only be used when converting to Tq
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

    // This should only be used when converting to Rq
    pub fn ntt_inv(&mut self) -> RingElement<P> {
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
        for item in &mut self.ring {
            *item = *item * 3303;
        }
        RingElement::new(self.ring)
    }

    // REMARKS:
    // make generic for NTT and Ring
    pub fn byte_encode(self) -> Vec<u8> {
        let mut out = Vec::with_capacity(256 * 12 / 8);

        for i in (0..self.ring.len()).step_by(2) {
            // Combine two 12-bit integers into a single 24-bit integer
            let x = u32::from(self.ring[i].val()) | (u32::from(self.ring[i + 1].val()) << 12);

            // Split the 24-bit integer into 3 bytes and append to the output vector
            out.push((x & 0xFF) as u8); // First 8 bits
            out.push(((x >> 8) & 0xFF) as u8); // Next 8 bits
            out.push(((x >> 16) & 0xFF) as u8); // Last 8 bits
        }
        out
    }
}

impl<P: ParameterSet + Copy> AddAssign for NttElement<P> {
    fn add_assign(&mut self, other: Self) {
        for (lhs, rhs) in self.ring.iter_mut().zip(other.ring.iter()) {
            *lhs += *rhs;
        }
    }
}

impl<P: ParameterSet + Copy> From<RingElement<P>> for NttElement<P> {
    fn from(mut val: RingElement<P>) -> Self {
        NttElement::new(&mut val)
    }
}

impl<P: ParameterSet + Copy> From<NttElement<P>> for RingElement<P> {
    fn from(mut val: NttElement<P>) -> Self {
        val.ntt_inv()
    }
}

impl<P: ParameterSet + Copy + core::cmp::PartialEq> PartialEq for NttElement<P> {
    fn eq(&self, other: &Self) -> bool {
        self.ring.iter().zip(other.ring.iter()).all(|(a, b)| a == b)
    }
}

impl<P: ParameterSet + Copy> Mul<NttElement<P>> for NttElement<P> {
    type Output = Self;

    fn mul(self, rhs: NttElement<P>) -> Self::Output {
        self.multiply_ntts(rhs)
    }
}

impl<P: ParameterSet + Copy> fmt::Debug for NttElement<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (index, element) in self.ring.iter().enumerate() {
            // adjust for spacing between rows
            write!(f, "{:<8}", element.val())?;
            // Adjust for modulus for row width
            if (index + 1) % 8 == 0 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}
