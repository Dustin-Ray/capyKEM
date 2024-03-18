use crate::{
    constants::{
        ml_kem_constants::{C2_SIZE, ENCODE_10, ENCODE_12},
        parameter_sets::ParameterSet,
    },
    math::{ntt_element::NttElement, ring_element::RingElement},
    Secret,
};
use alloc::vec::Vec;
use typenum::Unsigned;

impl<P: ParameterSet + Copy> Secret<P> {
    // TODO: make this return an error and handle them internally
    pub fn k_pke_decrypt(&self, dk_pke: &[u8], c: &[u8]) -> Vec<u8> {
        let mut slice = c;
        let mut u: Vec<RingElement<P>> = Vec::new();
        for _ in 0..P::K::to_usize() {
            let (current, next) = slice.split_at(ENCODE_10);
            // TODO: this should be generic for du
            let f = RingElement::decode_and_decompress_10(current).unwrap();
            u.push(f);
            slice = next;
        }

        let mut slice = dk_pke;
        let mut s_hat: Vec<NttElement<P>> = Vec::new();
        for _ in 0..P::K::to_usize() {
            let (current, next) = slice.split_at(ENCODE_12);
            let f = NttElement::byte_decode_12(current).unwrap();
            s_hat.push(f);
            slice = next;
        }

        // TODO: handle this error
        // TODO: this should be generic for dv but it might be already
        let v = RingElement::<P>::decode_and_decompress_4(&c[c.len() - C2_SIZE..c.len()]).unwrap();

        let mut y = RingElement::<P>::zero();
        for i in 0..s_hat.len() {
            y += (s_hat[i] * u[i].into()).into();
        }

        let w = v - y;
        let mut s = [0_u8; 32];
        w.compress_and_encode_1(&mut s);
        s.to_vec()
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        constants::parameter_sets::{KEM_1024, KEM_512, KEM_768},
        Secret,
    };

    #[test]
    fn roundtrip() {
        let s: Secret<KEM_512> = Secret::new([4_u8; 32]);
        let r = &[42_u8; 32];
        let (ek, dk_pke) = s.k_pke_keygen(&r);
        let c = s.k_pke_encrypt(&ek, r);
        let dec = s.k_pke_decrypt(&dk_pke, &c);
        assert_eq!(dec, s.m);
    }
}
