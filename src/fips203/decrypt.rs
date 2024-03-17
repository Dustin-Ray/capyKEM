use alloc::vec::Vec;

use crate::{constants::parameter_sets::ParameterSet, math::ring_element::RingElement, Secret};

impl<P: ParameterSet + Copy> Secret<P> {
    // TODO: make this return an error
    pub fn k_pke_decrypt(&self, dk_pke: &[u8], mut c: &[u8]) -> Vec<u8> {
        let mut u: Vec<RingElement<P>> = Vec::new();
        let slice = c;

        let encoding_size_10: usize = 320;

        for _ in 0..P::k {
            // Attempt to decode and decompress the first portion of `c`.
            // `decode_and_decompress_10` is assumed to be a function that returns a `RingElement`
            // and consumes exactly `ENCODING_SIZE_10` bytes from the beginning of `c`.
            let (current, next) = c.split_at(320);
            let f = RingElement::decode_and_decompress_10(current).unwrap();

            u.push(f);
            c = next; // Move to the next segment of `c`.
        }
        // println!("{:?}", u);
        vec![]
    }
}

#[cfg(test)]
mod tests {
    use crate::{constants::parameter_sets::P768, Secret};

    #[test]
    fn test_decrypt() {
        let m_768: Secret<P768> = Secret::new([4_u8; 32]);
        let rand: &[u8; 32] = &[4_u8; 32];

        let (ek, dk) = m_768.k_pke_keygen(rand);
        // pretty_print_vec_u8(&ek);

        let c = m_768.k_pke_encrypt(&ek, rand);
        // pretty_print_vec_u8(&c);
        let m = m_768.k_pke_decrypt(&dk, &c);
    }

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
}
