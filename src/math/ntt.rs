use super::field_element::FieldElement as F;
use crate::constants::{ml_kem_constants::Q, K_MOD_ROOTS, K_NTT_ROOTS};
use core::fmt;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
#[derive(Clone)]
pub struct NttElement([F; 256]);

impl fmt::Debug for NttElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (index, element) in self.0.iter().enumerate() {
            // adjust for space between columns
            write!(f, "{:<8}", element.val)?;
            // adust for row width
            if (index + 1) % 16 == 0 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

impl NttElement {
    fn default() -> Self {
        NttElement([F::default(); 256])
    }

    /// byte stream b should be rho||i||j
    fn sample(byte_stream: Vec<u8>) -> Self {
        // Get the XOF for the input
        let mut prf = Shake256::default();
        prf.update(&byte_stream);

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
                a_hat[j] = F::new(d_1).reduce_once();
                j += 1;
            }
            if d_2 < Q && j < 256 {
                a_hat[j] = F::new(d_2).reduce_once();
                j += 1;
            }
            i += 3
        }
        NttElement(a_hat)
    }

    fn multiply_ntts(&self, other: Self) -> Self {
        let mut h_hat = NttElement::default();
    
        // Iterate over `K_MOD_ROOTS` with their indices
        for (i, &k_mod_root) in K_MOD_ROOTS.iter().enumerate() {
            (h_hat.0[2 * i], h_hat.0[2 * i + 1]) = NttElement::base_case_multiply(
                self.0[2 * i],
                self.0[2 * i + 1],
                other.0[2 * i],
                other.0[2 * i + 1],
                k_mod_root,
            )
        }
    
        h_hat
    }

    fn base_case_multiply(a_0: F, a_1: F, b_0: F, b_1: F, gamma: u16) -> (F, F) {
        let c_0 = a_0 * b_0 + a_1 * b_1 * gamma;
        let c_1 = a_0 * b_1 + a_1 * b_0;
        (c_0, c_1)
    }

    fn ntt(&mut self) {
        let mut k = 1;
        let mut len = 128;
        while len >= 2 {
            for start in (0..256).step_by(2 * len) {
                let zeta = K_NTT_ROOTS[k];
                k += 1;

                for j in start..start + len {
                    let t = zeta * self.0[j + len] % Q;
                    self.0[j + len] = self.0[j] - F::new(t).reduce_once();
                    self.0[j] = self.0[j] + F::new(t).reduce_once();
                }
            }
            len /= 2;
        }
    }

    fn ntt_inv(&mut self) {
        let mut k = 127;
        let mut len = 2;
        while len <= 128 {
            for start in (0..256).step_by(2 * len) {
                let zeta = K_NTT_ROOTS[k];
                k -= 1;

                for j in start..start + len {
                    let t = self.0[j];
                    self.0[j] = t + self.0[j + len].reduce_once();
                    self.0[j + len] = F::new(zeta * (self.0[j + len] - t)).reduce_once();
                }
            }
            len *= 2;
        }
        for item in self.0.iter_mut() {
            *item = (*item * 3303).reduce_once();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let mut byte_stream = NttElement::sample(vec![42_u8; 32]);
        let byte_stream_copy = byte_stream.clone();
        byte_stream.ntt();
        byte_stream.ntt_inv();
        assert_eq!(byte_stream_copy.0, byte_stream.0)
    }
}
