use alloc::vec;
use alloc::vec::Vec;

use crate::{
    constants::parameter_sets::ParameterSet,
    math::{ntt_element::NttElement, ring_element::RingElement},
    Secret,
};

impl<P: ParameterSet + Copy> Secret<P> {
    fn k_pke_encrypt(&self, ek_pke: &[u8], rand: &[u8; 32]) -> Vec<u8> {
        let k = P::K as usize;
        let mut n = 0;
        // handle this error
        let mut t_hat = Vec::with_capacity(k);
        let mut ek_slice = &ek_pke[..];

        for _ in 0..k {
            let decoded_element =
                NttElement::<P>::poly_byte_decode(&ek_slice[0..384]).map_err(|err| err.to_string());

            // Append the decoded element to t
            t_hat.push(decoded_element);
            ek_slice = &ek_slice[384..];
        }

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

        let mut u = vec![RingElement::zero(); k];
        for i in 0..k {
            u[i] = e_1[i];
            for j in 0..k {
                u[i] += (a_hat_transpose[i * k + j] * r_hat[j]).into(); // into ring = NTT^(-1)
            }
        }
        println!("{:?}", u);
        vec![]
    }
}

#[cfg(test)]
mod tests {
    use crate::{constants::parameter_sets::P768, Secret};

    fn pretty_print_vec_u8(vec: &Vec<u8>) {
        for (index, &element) in vec.iter().enumerate() {
            print!("{:<8}", element);
            if (index + 1) % 8 == 0 {
                println!();
            }
        }
        if !vec.is_empty() && vec.len() % 16 != 0 {
            println!();
        }
    }

    #[test]
    fn smoke_test_instantiate_over_parameter_sets() {
        // let m_512: Secret<P512> = Secret::new([0_u8; 32]);
        // let (ek, dk) = m_512.k_pke_keygen(&[0_u8; 32]);

        let s: Secret<P768> = Secret::new([0_u8; 32]);
        let ek: &[u8] = &[0_u8; 1184];
        let r: &[u8; 32] = &[0_u8; 32];

        s.k_pke_encrypt(ek, r);

        // let m_1024: Secret<P1024> = Secret::new([0_u8; 32]);
        // let (ek, dk) = m_1024.k_pke_keygen(&[0_u8; 32]);
    }
}
