use crate::{
    constants::{
        ml_kem_constants::{C2_SIZE, ENCODE_10, ENCODE_12},
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
use typenum::Unsigned;
use typenum::U1;

impl<P: ParameterSet + Copy> Secret<P> {
    // TODO: make this return an error and handle them internally
    pub fn k_pke_decrypt(&self, dk_pke: &[u8], c: &[u8]) -> Vec<u8> {
        let mut slice = c;
        let mut u: Vec<RingElement> = Vec::new();
        for _ in 0..P::K::to_usize() {
            let (current, next) = slice.split_at(ENCODE_10);
            let mut f: RingElement = Encode::<P::Du>::decode(current);
            f.decompress::<P::Du>();

            u.push(f);
            slice = next;
        }

        let mut slice = dk_pke;
        let mut s_hat: Vec<NttElement> = Vec::new();
        for _ in 0..P::K::to_usize() {
            let (current, next) = slice.split_at(ENCODE_12);
            // let mut f: NttElement = Encode::<U12>::decode(current);
            // f.decompress::<U12>();
            let f = NttElement::byte_decode_12(current).unwrap();
            s_hat.push(f);
            slice = next;
        }

        let mut v: RingElement = Encode::<P::Dv>::decode(&c[c.len() - C2_SIZE..c.len()]);
        v.decompress::<P::Dv>();

        let mut y = RingElement::zero();
        for i in 0..s_hat.len() {
            y += (s_hat[i] * u[i].into()).into();
        }

        let mut w = v - y;
        let s = Encode::<U1>::encode(w.compress::<U1>());
        s.to_vec()
    }
}
