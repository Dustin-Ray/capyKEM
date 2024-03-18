use alloc::vec::Vec;

use crate::{
    constants::{
        ml_kem_constants::{self, CIPHERTEXT_SIZE, ENCODE_12},
        parameter_sets::ParameterSet,
    },
    math::{ntt_element::NttElement, ring_element::RingElement},
    Secret,
};

impl<P: ParameterSet + Copy> Secret<P> {
    pub fn k_pke_encrypt(&self, ek_pke: &[u8], rand: &[u8; 32]) -> Vec<u8> {
        // TODO: parameterize this
        const k: usize = ml_kem_constants::k;
        let mut n = 0;
        // handle this error
        let mut t_hat = [NttElement::<P>::zero(); k];

        // TODO: parameterize 384 as encodesize12
        for i in 0..t_hat.len() {
            t_hat[i] = NttElement::<P>::byte_decode_12(&ek_pke[i * ENCODE_12..(i + 1) * ENCODE_12])
                .unwrap();
        }

        let rho: &[u8] = &ek_pke[ENCODE_12 * k..(ENCODE_12 * k) + 32];

        // Generate the matrix a_hat^T
        let mut a_hat_transpose = [NttElement::<P>::zero(); k * k];
        for i in 0..k {
            for j in 0..k {
                a_hat_transpose[i * k + j] = NttElement::sample_ntt(rho, i, j);
            }
        }

        // generate r, run ntt k times
        let mut r_hat = [NttElement::<P>::zero(); k];
        for r_elem in r_hat.iter_mut().take(k) {
            *r_elem = RingElement::sample_poly_cbd(rand, n).into();
            n += 1;
        }

        // generate e1
        let mut e_1 = [RingElement::<P>::zero(); k];
        for e_elem in e_1.iter_mut().take(k) {
            *e_elem = RingElement::sample_poly_cbd(rand, n);
            n += 1;
        }

        // sample e2
        let e2: RingElement<P> = RingElement::sample_poly_cbd(rand, n);

        let u: Vec<RingElement<P>> = e_1
            .iter()
            .enumerate()
            .map(|(i, e1_elem)| {
                let sum: RingElement<P> = (0..k)
                    .map(|j| (a_hat_transpose[i * k + j] * r_hat[j]).into())
                    .sum();
                *e1_elem + sum
            })
            .collect();

        // TODO: handle this error
        let mu = RingElement::<P>::decode_decompress_1(&self.m).unwrap();

        let mut v = NttElement::<P>::zero();
        for i in 0..t_hat.len() {
            v += t_hat[i] * r_hat[i];
        }
        let mut v = v.ntt_inv();
        v += e2;
        v += mu;

        let mut c: Vec<u8> = Vec::with_capacity(CIPHERTEXT_SIZE);
        for f in u.iter() {
            c = RingElement::compress_and_encode_10(c, *f);
        }
        c = RingElement::compress_and_encode_4(c, v);

        c
    }
}
