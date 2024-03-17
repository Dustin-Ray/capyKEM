use crate::{
    constants::parameter_sets::ParameterSet,
    math::{ntt_element::NttElement, ring_element::RingElement},
    Secret,
};
use alloc::vec;
use alloc::vec::Vec;
use sha3::{Digest, Sha3_512};

// NOTES:
// TODO: determine difference between PRF and XOF here
impl<P: ParameterSet + Copy> Secret<P> {
    // should be private
    pub fn k_pke_keygen(&self, d: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
        let mut hasher = Sha3_512::default();
        hasher.update(d);
        let binding = hasher.finalize();
        let b = binding.as_slice();

        // (ρ, σ ) <- G(d)
        let rho: &[u8] = &b[0..32];
        let sigma = &b[32..64];

        let k = P::k as usize;
        let mut n = 0;

        // Generate the matrix a_hat
        let mut a_hat = vec![NttElement::zero(); k * k];
        for i in 0..k {
            for j in 0..k {
                // see: https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/s-C-zIAeKfE/m/eZJmXYsSAQAJ?
                a_hat[i * k + j] = NttElement::sample_ntt(rho, j, i);
            }
        }

        // generate s
        let mut s_hat: Vec<NttElement<P>> = vec![NttElement::zero(); k];
        for s_elem in s_hat.iter_mut().take(k) {
            *s_elem = RingElement::sample_poly_cbd(sigma, n).into();
            n += 1;
        }

        // dbg!(s_hat.clone());

        // generate e
        let mut e_hat = vec![NttElement::zero(); k];
        for e_elem in e_hat.iter_mut().take(k) {
            *e_elem = RingElement::sample_poly_cbd(sigma, n).into();
            n += 1;
        }

        // t_hat = A o s_hat + e_hat
        let mut t = vec![NttElement::zero(); k];
        for i in 0..t.len() {
            t[i] = e_hat[i];
            for j in 0..s_hat.len() {
                t[i] += a_hat[i * k + j] * s_hat[j];
            }
        }

        // ByteEncode12(t_hat||rho)
        let mut ek_pke: Vec<u8> = Vec::with_capacity(1184);
        for &item in t.iter() {
            ek_pke = item.byte_encode_12(ek_pke);
        }
        ek_pke.append(&mut rho.into());

        let mut dk_pke: Vec<u8> = Vec::with_capacity(1152);

        for &item in s_hat.iter() {
            dk_pke = item.byte_encode_12(dk_pke);
        }
        // pretty_print_vec_u8(&dk_pke);
        (ek_pke, dk_pke)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        constants::{
            parameter_sets::P768,
            pke_keygen_dk_pke, pke_keygen_ek_pke,
        },
        Secret,
    };

    #[test]
    fn smoke_test_instantiate_over_parameter_set() {
        let m_768: Secret<P768> = Secret::new([4_u8; 32]);
        let (ek, dk) = m_768.k_pke_keygen(&[4_u8; 32]);

        assert_eq!(ek, pke_keygen_ek_pke);
        assert_eq!(dk, pke_keygen_dk_pke);
    }
}

fn pretty_print_vec_u8(vec: &[u8]) {
    for (index, &element) in vec.iter().enumerate() {
        print!("{:<8}", element);
        if (index + 1) % 8 == 0 {
            println!();
        }
    }
    // Handle the case where the Vec doesn't end exactly at a row boundary
    if !vec.is_empty() && vec.len() % 16 != 0 {
        println!(); // Ensure there's a newline at the end if needed
    }
}

fn compare_arrays(a: &[u8], b: &[u8]) {
    if a.len() != b.len() {
        println!("Arrays have different lengths");
    }

    for (i, (&val1, &val2)) in a.iter().zip(b.iter()).enumerate() {
        if val1 != val2 {
            println!("Difference found at index {}: {} vs {}", i, val1, val2);
        }
    }
}
