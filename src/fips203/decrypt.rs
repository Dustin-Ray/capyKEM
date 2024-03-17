use alloc::vec::Vec;

use crate::{
    constants::parameter_sets::ParameterSet,
    math::{ntt_element::NttElement, ring_element::RingElement},
    Secret,
};

impl<P: ParameterSet + Copy> Secret<P> {
    // TODO: make this return an error and handle them internally
    pub fn k_pke_decrypt(&self, dk_pke: &[u8], c: &[u8]) -> Vec<u8> {
        // TODO: make this parameter
        let encoding_size_10: usize = 320;
        let mut slice = c;
        let mut u: Vec<RingElement<P>> = Vec::new();
        for _ in 0..P::k {
            let (current, next) = slice.split_at(encoding_size_10);
            let f = RingElement::decode_and_decompress_10(current).unwrap();
            u.push(f);
            slice = next;
        }
        
        // TODO: make this parameter
        let encoding_size_12 = 384;
        let mut slice = dk_pke;
        let mut s_hat: Vec<NttElement<P>> = Vec::new();
        for _ in 0..P::k {
            let (current, next) = slice.split_at(encoding_size_12);
            let f = NttElement::byte_decode_12(current).unwrap();
            s_hat.push(f);
            slice = next;
        }

        //TODO: handle this error
        let v = RingElement::<P>::decode_and_decompress_4(&c[c.len() - 128..c.len()]).unwrap();

        let mut y = RingElement::<P>::zero();
        for i in 0..s_hat.len() {
            y += (s_hat[i] * u[i].into()).into();
        }

        let w = v - y;
        w.compress_and_encode_1()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        constants::
            parameter_sets::P768
        ,
        Secret,
    };

    #[test]
    fn roundtrip() {
        let s: Secret<P768> = Secret::new([4_u8; 32]);
        let r: &[u8; 32] = &[4_u8; 32];
        let (ek, dk_pke) = s.k_pke_keygen(&[4_u8; 32]);
        let c = s.k_pke_encrypt(&ek, r);
        let dec = s.k_pke_decrypt(&dk_pke, &c);
        assert_eq!(dec, s.m);
    }
}
