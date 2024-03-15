use super::{field_element::FieldElement as F, ring_element::RingElement};
use crate::constants::ml_kem_constants::{self, N, Q};
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
    pub coefficients: [F<P>; 256],
}

impl<P: ParameterSet + Copy> NttElement<P> {
    pub fn new(r: &mut RingElement<P>) -> Self {
        let mut ntt_el = NttElement {
            coefficients: r.coefficients,
        };
        ntt_el.ntt();
        ntt_el
    }

    pub fn zero() -> Self {
        NttElement {
            coefficients: [F::zero(); 256],
        }
    }

    pub fn get_ring(&self) -> [F<P>; 256] {
        self.coefficients
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
                a.coefficients[j] = F::new(d1);
                j += 1;
            }
            if j >= N.into() {
                break;
            }

            if d2 < Q {
                a.coefficients[j] = F::new(d2);
                j += 1;
            }
        }
        a
    }

    fn multiply_ntts(&self, other: Self) -> Self {
        let mut h_hat = NttElement::zero();

        // Iterate over `K_MOD_ROOTS` with their indices
        for (i, &k_mod_root) in K_MOD_ROOTS.iter().enumerate() {
            (h_hat.coefficients[2 * i], h_hat.coefficients[2 * i + 1]) =
                NttElement::base_case_multiply(
                    self.coefficients[2 * i],
                    self.coefficients[(2 * i) + 1],
                    other.coefficients[2 * i],
                    other.coefficients[(2 * i) + 1],
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
                    let t = zeta * self.coefficients[j + len] % Q;
                    self.coefficients[j + len] = self.coefficients[j] - F::new(t);
                    self.coefficients[j] += F::new(t);
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
                    let t = self.coefficients[j];
                    self.coefficients[j] = t + self.coefficients[j + len];
                    self.coefficients[j + len] = F::new(zeta * (self.coefficients[j + len] - t));
                }
            }
            len *= 2;
        }
        for item in &mut self.coefficients {
            *item = *item * 3303;
        }
        RingElement::new(self.coefficients)
    }

    // REMARKS:
    // make generic for NTT and Ring
    pub fn poly_byte_encode(self) -> Vec<u8> {
        let mut out = Vec::with_capacity(256 * 12 / 8);

        for i in (0..self.coefficients.len()).step_by(2) {
            // Combine two 12-bit integers into a single 24-bit integer
            let x = u32::from(self.coefficients[i].val())
                | (u32::from(self.coefficients[i + 1].val()) << 12);

            // Split the 24-bit integer into 3 bytes and append to the output vector
            out.push((x & 0xFF) as u8); // First 8 bits
            out.push(((x >> 8) & 0xFF) as u8); // Next 8 bits
            out.push(((x >> 16) & 0xFF) as u8); // Last 8 bits
        }
        out
    }

    pub fn poly_byte_decode(b: &[u8]) -> Result<Self, String> {
        const MASK_12: u32 = 0b1111_1111_1111;
        if b.len() != (ml_kem_constants::ENCODE_SIZE_12).into() {
            return Err("Invalid encoding length".to_owned());
        }

        let mut f = Vec::with_capacity(N.into());

        let mut i = 0;
        while i < b.len() {
            let d = u32::from(b[i]) | (u32::from(b[i + 1]) << 8) | (u32::from(b[i + 2]) << 16);
            let elem1 = F::new((d & MASK_12) as u16)
                .check_reduced()
                .map_err(|_| "Invalid polynomial encoding")?;

            let elem2 = F::new((d >> 12) as u16)
                .check_reduced()
                .map_err(|_| "Invalid polynomial encoding")?;

            f.push(elem1);
            f.push(elem2);

            i += 3;
        }
        let coefficients: [F<P>; 256] = f
            .try_into()
            .map_err(|_| "Conversion to fixed-size array failed")?;

        Ok(NttElement { coefficients })
    }
}

impl<P: ParameterSet + Copy> AddAssign for NttElement<P> {
    fn add_assign(&mut self, other: Self) {
        for (lhs, rhs) in self.coefficients.iter_mut().zip(other.coefficients.iter()) {
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
        self.coefficients
            .iter()
            .zip(other.coefficients.iter())
            .all(|(a, b)| a == b)
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
        for (index, element) in self.coefficients.iter().enumerate() {
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
