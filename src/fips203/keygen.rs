use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::{
    constants::parameter_sets::{ParameterSet, P768},
    math::{ntt::NttElement, ring_element::RingElement},
    Message,
};

// NOTES:
// TODO: determine difference between PRF and XOF here
impl Message {
    fn k_pke_keygen<P: ParameterSet>(&self) -> (Vec<u8>, Vec<u8>) {
        let mut xof = Shake256::default(); // I think this should be 512 -> 256 bits of security
        let bytes: Vec<u8> = (0..32).map(|_| ChaCha20Rng::from_entropy().gen()).collect();

        xof.update(&bytes);
        let mut b = [0_u8; 64];
        let mut reader = xof.finalize_xof();
        reader.read(&mut b);

        // (ρ, σ ) <- G(d)
        let rho: &[u8] = &b[32..64];
        let sigma = &b[32..64];
        let k = self.k as usize;
        let mut n = 0;

        // Generate the matrix a_hat
        let mut a_hat = vec![NttElement::zero(); k * k];
        for i in 0..k {
            for j in 0..k {
                a_hat[i * k + j] = NttElement::sample_ntt(rho, i, j);
            }
        }

        // generate s
        let mut s = vec![NttElement::zero(); k];
        for s_elem in s.iter_mut().take(k) {
            *s_elem = RingElement::sample_poly_cbd::<P>(sigma, n).into();
            n += 1;
        }

        // generate e
        let mut e = vec![NttElement::zero(); k];
        for e_elem in e.iter_mut().take(k) {
            *e_elem = RingElement::sample_poly_cbd::<P>(sigma, n).into();
            n += 1;
        }

        // t_hat = A o s_hat + e_hat
        let mut t = vec![NttElement::zero(); k];
        for i in 0..t.len() {
            t[i] = e[i];
            for j in 0..s.len() {
                // Is addition happening in NTT domain as well?
                t[i] += a_hat[i * k + j] * s[j] // This accounts for FIPS 203 revision and will likely change
            }
        }

        // ByteEncode12(t_hat||rho)
        let ek_pke: Vec<u8> = t
            .iter_mut()
            .take(k)
            .flat_map(|t_elem| {
                let mut bytes = Into::<RingElement>::into(*t_elem).byte_encode();
                bytes.extend_from_slice(rho);
                bytes
            })
            .collect();

        let dk_pke: Vec<u8> = s
            .iter_mut()
            .take(k)
            .flat_map(|s_elem| Into::<RingElement>::into(*s_elem).byte_encode())
            .collect();

        (ek_pke, dk_pke)
    }
}

#[cfg(test)]
mod tests {
    use crate::{constants::parameter_sets::P768, Message};

    #[test]
    fn test_keygen() {
        let m = Message {
            m: vec![],
            k: 3,
            du: 10,
            dv: 4,
            eta_1: 2,
            eta_2: 2,
        };
        let (ek, dk) = m.k_pke_keygen::<P768>();

        // Convert to hexadecimal string and print
        println!("Encryption Key: {}", hex::encode(&ek));
        println!("Decryption Key: {}", hex::encode(&dk));
    }
}