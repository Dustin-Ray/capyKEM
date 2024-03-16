use alloc::vec::Vec;
use alloc::{string::ToString, vec};

use crate::{
    constants::parameter_sets::ParameterSet,
    math::{ntt_element::NttElement, ring_element::RingElement},
    Secret,
};

impl<P: ParameterSet + Copy> Secret<P> {
    fn k_pke_encrypt(&self, ek_pke: &[u8], rand: &[u8; 32]) -> Vec<u8> {
        let k = P::k as usize;
        let mut n = 0;
        // handle this error
        let mut t_hat: Vec<NttElement<P>> = Vec::with_capacity(k);

        // Process each chunk of 384 bytes to decode it into an NttElement
        // TODO: parameterize 384
        let decoded_ek_pke: Result<Vec<_>, _> = ek_pke
            .chunks_exact(384) // Create an iterator that yields slices of ek_pke in chunks of 384 bytes.
            .take(k) // Ensure we only take k elements
            .map(|chunk| NttElement::<P>::byte_decode_12(chunk).map_err(|e| e.to_string()))
            .collect(); // Collect results into a vector, handling errors automatically.

        t_hat.extend(decoded_ek_pke.unwrap());

        let rho: &[u8] = &ek_pke[384 * k..(384 * k) + 32];

        // Generate the matrix a_hat^T
        let mut a_hat_transpose = vec![NttElement::<P>::zero(); k * k];
        for i in 0..k {
            for j in 0..k {
                a_hat_transpose[i * k + j] = NttElement::sample_ntt(rho, i, j);
            }
        }

        // generate r, run ntt k times
        let mut r_hat = vec![NttElement::<P>::zero(); k];
        for r_elem in r_hat.iter_mut().take(k) {
            *r_elem = RingElement::sample_poly_cbd(rand, n).into();
            n += 1;
        }

        // generate e1
        let mut e_1 = vec![RingElement::<P>::zero(); k];
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

        let mut c: Vec<u8> = Vec::with_capacity(1088);
        for f in u.iter() {
            c = RingElement::compress_and_encode_10(c, *f);
        }
        c = RingElement::compress_and_encode_4(c, v);

        c
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        constants::{kpke_encrypt_result, parameter_sets::P768},
        Secret,
    };

    // fn pretty_print_vec_u8(vec: &Vec<u8>) {
    //     for (index, &element) in vec.iter().enumerate() {
    //         print!("{:<8}", element);
    //         if (index + 1) % 8 == 0 {
    //             println!();
    //         }
    //     }
    //     if !vec.is_empty() && vec.len() % 16 != 0 {
    //         println!();
    //     }
    // }

    #[test]
    fn smoke_test_instantiate_over_parameter_sets() {
        let s: Secret<P768> = Secret::new([4_u8; 32]);
        let ek: &[u8] = &[0_u8; 1184];
        let r: &[u8; 32] = &[0_u8; 32];

        let res = s.k_pke_encrypt(ek, r);
        assert_eq!(res, kpke_encrypt_result)
    }
}
