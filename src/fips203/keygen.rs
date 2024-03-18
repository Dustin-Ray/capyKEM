use crate::{
    constants::{
        ml_kem_constants::{self, D_PKE_KEYSIZE, E_PKE_KEYSIZE},
        parameter_sets::ParameterSet,
    },
    math::{ntt_element::NttElement, ring_element::RingElement},
    Secret,
};
use alloc::vec::Vec;
use sha3::{Digest, Sha3_512};

impl<P: ParameterSet + Copy> Secret<P> {
    pub fn k_pke_keygen(&self, d: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
        let mut hasher = Sha3_512::default();
        hasher.update(d);
        let binding = hasher.finalize();
        let b = binding.as_slice();

        // (ρ, σ ) <- G(d)
        let rho: &[u8] = &b[0..32];
        let sigma = &b[32..64];

        const k: usize = ml_kem_constants::k;
        let mut n = 0;

        // Generate the matrix a_hat
        let mut a_hat: [NttElement<P>; 9] = [NttElement::zero(); k * k];
        for i in 0..k {
            for j in 0..k {
                // see: https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/s-C-zIAeKfE/m/eZJmXYsSAQAJ?
                a_hat[i * k + j] = NttElement::sample_ntt(rho, j, i);
            }
        }

        // generate s
        let mut s_hat = [NttElement::zero(); k];
        for s_elem in s_hat.iter_mut().take(k) {
            *s_elem = RingElement::sample_poly_cbd(sigma, n).into();
            n += 1;
        }

        // generate e
        let mut e_hat = [NttElement::zero(); k];
        for e_elem in e_hat.iter_mut().take(k) {
            *e_elem = RingElement::sample_poly_cbd(sigma, n).into();
            n += 1;
        }

        // t_hat = A o s_hat + e_hat
        let mut t = [NttElement::zero(); k];
        for i in 0..t.len() {
            t[i] = e_hat[i];
            for j in 0..s_hat.len() {
                t[i] += a_hat[i * k + j] * s_hat[j];
            }
        }

        // ByteEncode12(t_hat||rho)
        let mut ek_pke: Vec<u8> = Vec::with_capacity(E_PKE_KEYSIZE);
        for &item in t.iter() {
            ek_pke = item.byte_encode_12(ek_pke);
        }
        ek_pke.append(&mut rho.into());

        let mut dk_pke: Vec<u8> = Vec::with_capacity(D_PKE_KEYSIZE);
        for &item in s_hat.iter() {
            dk_pke = item.byte_encode_12(dk_pke);
        }
        // pretty_print_vec_u8(&dk_pke);
        (ek_pke, dk_pke)
    }
}
