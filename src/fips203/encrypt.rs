use crate::{
    constants::{
        ml_kem_constants::{self, ENCODE_12},
        parameter_sets::ParameterSet,
    },
    math::{
        encoding::{Compress, Encode},
        ntt_element::NttElement,
        ring_element::RingElement,
    },
    Secret,
};
use alloc::vec::Vec;
use typenum::U1;

impl<P: ParameterSet + Copy> Secret<P> {
    pub fn k_pke_encrypt(&self, ek_pke: &[u8], rand: &[u8; 32]) -> Vec<u8> {
        // TODO: parameterize this
        const k: usize = ml_kem_constants::k;
        let mut n = 0;
        // handle this error
        let mut t_hat = [NttElement::zero(); k];

        // TODO: parameterize 384 as encodesize12
        for i in 0..t_hat.len() {
            t_hat[i] =
                NttElement::byte_decode_12(&ek_pke[i * ENCODE_12..(i + 1) * ENCODE_12]).unwrap();
            // Encode::<U12>::decode(&ek_pke[i * ENCODE_12..(i + 1) * ENCODE_12]);
            // t_hat[i].decompress::<U12>();
        }

        let rho: &[u8] = &ek_pke[ENCODE_12 * k..(ENCODE_12 * k) + 32];

        // Generate the matrix a_hat^T
        let mut a_hat_transpose = [NttElement::zero(); k * k];
        for i in 0..k {
            for j in 0..k {
                a_hat_transpose[i * k + j] = NttElement::sample_ntt(rho, i, j);
            }
        }

        // generate r, run ntt k times
        let mut r_hat = [NttElement::zero(); k];
        for r_elem in r_hat.iter_mut().take(k) {
            *r_elem = RingElement::sample_poly_cbd(rand, n).into();
            n += 1;
        }

        // generate e1
        let mut e_1 = [RingElement::zero(); k];
        for e_elem in e_1.iter_mut().take(k) {
            *e_elem = RingElement::sample_poly_cbd(rand, n);
            n += 1;
        }

        // sample e2
        let e2: RingElement = RingElement::sample_poly_cbd(rand, n);

        let mut u: Vec<RingElement> = e_1
            .iter()
            .enumerate()
            .map(|(i, e1_elem)| {
                let sum: RingElement = (0..k)
                    .map(|j| (a_hat_transpose[i * k + j] * r_hat[j]).into())
                    .sum();
                *e1_elem + sum
            })
            .collect();

        let mut mu: RingElement = Encode::<U1>::decode(&self.m);
        mu.decompress::<U1>();

        let mut v = NttElement::zero();
        for i in 0..t_hat.len() {
            v += t_hat[i] * r_hat[i];
        }
        let mut v = v.ntt_inv();
        v += e2;
        v += mu;

        let mut c: Vec<u8> = Vec::new();

        for ring in u.iter_mut() {
            let bytes = &mut Encode::<P::Du>::encode(ring.compress::<P::Du>());
            c.append(bytes);
        }

        c.append(&mut Encode::<P::Dv>::encode(v.compress::<P::Dv>()));

        c
    }
}
