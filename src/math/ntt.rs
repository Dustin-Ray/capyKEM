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
    ring: [F<P>; 256],
}

impl<P: ParameterSet + Copy> NttElement<P> {
    fn new(r: &mut RingElement<P>) -> Self {
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
            )
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
    fn ntt_inv(&mut self) -> RingElement<P> {
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

    use crate::constants::parameter_sets::P768;

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
        let byte_stream = [0_u8; 0];

        let a = NttElement::sample_ntt(&byte_stream, 0, 1);

        // testing against the great Filippo Valsorda https://github.com/FiloSottile/mlkem768
        let result: [F<P768>; 256] = [
            2278, 18, 2449, 1376, 2453, 1346, 66, 738, 2100, 1008, 950, 2669, 2121, 3030, 880,
            2569, 3146, 1432, 1285, 2106, 1943, 895, 2326, 3255, 1301, 1752, 1281, 2500, 3149,
            1061, 959, 687, 199, 1817, 1651, 2069, 3091, 2864, 120, 2222, 3005, 1823, 2721, 3012,
            665, 1426, 386, 1639, 1632, 591, 1405, 756, 464, 1405, 2701, 3275, 76, 2137, 664, 2457,
            2216, 2352, 1994, 1521, 1944, 1753, 999, 2051, 3219, 2771, 1596, 2123, 527, 339, 2532,
            2079, 2994, 576, 1876, 2698, 1708, 119, 537, 2122, 3132, 285, 3198, 3131, 2761, 3187,
            1, 3082, 2809, 3140, 895, 356, 1653, 2663, 2856, 2290, 3166, 1245, 1876, 2355, 2746,
            3213, 619, 551, 3216, 2092, 966, 479, 3079, 2557, 2706, 380, 2388, 915, 4, 2336, 144,
            3220, 1807, 95, 1109, 2105, 1441, 2379, 2890, 2985, 2129, 1040, 1472, 1350, 1976, 927,
            862, 1556, 2188, 447, 856, 1458, 2372, 1254, 2132, 2618, 200, 2880, 2834, 1811, 505,
            124, 621, 2574, 2546, 2974, 1875, 1646, 618, 1867, 1394, 1059, 486, 1232, 2574, 563,
            2509, 2805, 2674, 1594, 782, 1147, 12, 1853, 459, 2718, 1861, 913, 2538, 1986, 346,
            2139, 1256, 3148, 830, 615, 676, 2220, 2638, 893, 977, 474, 1096, 1307, 3285, 462,
            3082, 2805, 1286, 2645, 2733, 2695, 2082, 3216, 414, 1376, 2636, 971, 2671, 1721, 746,
            516, 1620, 688, 1903, 2497, 2869, 1587, 819, 256, 2326, 943, 1733, 117, 2941, 2933,
            1852, 2753, 2057, 2585, 1042, 2572, 220, 3049, 558, 2617, 1975, 45, 2593, 757, 3202,
            1164, 1123, 1458, 1720, 2365, 148, 605, 2229, 760, 90, 3212, 3015, 1643, 1962, 2954,
        ]
        .map(|val| F::new(val));

        assert_eq!(a.ring, result);
    }

    #[test]
    fn test_ntt() {
        // sample output is in NTT domain
        let mut byte_stream: NttElement<P768> = NttElement::sample_ntt(&vec![42_u8; 32], 1, 1);
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
        let mut ring_element: RingElement<P768> = RingElement::sample_poly_cbd(&bytes, 0xFF);
        let ring_element_copy = ring_element;

        // runs .ntt() on intstantiation
        let ntt_element = NttElement::new(&mut ring_element);
        assert_eq!(
            ring_element_copy.coefficients,
            Into::<RingElement<P768>>::into(ntt_element).coefficients
        );
    }

    #[test]
    fn test_multiply_ntts_associative() {
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        let a: NttElement<P768> = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
        let b = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xBB));
        let c = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xCC));

        // Test associativity (ab)c = a(bc)
        let ab_c = (a * b) * c;
        let a_bc = a * (b * c);
        assert_eq!(ab_c.ring, a_bc.ring);

        let a: NttElement<P768> = NttElement::sample_ntt(&bytes.clone(), 0, 0);
        let b = NttElement::sample_ntt(&bytes.clone(), 0, 0);
        let c = NttElement::sample_ntt(&bytes.clone(), 0, 0);

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

        let a: NttElement<P768> = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
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

            let a: NttElement<P768> =
                NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
            let b = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xBB));

            let result = a * b;
            assert!(
                result.ring.iter().all(|x| x.val() < 3329),
                "Result of multiplication must be valid NttElement"
            );
        }
    }
}
