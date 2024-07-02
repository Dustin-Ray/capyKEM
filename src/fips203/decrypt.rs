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
};
use alloc::vec::Vec;
use sha3::{Digest, Sha3_512};
use typenum::Unsigned;
use typenum::U1;

use super::encrypt::k_pke_encrypt;

pub fn mlkem_decaps<P: ParameterSet>(c: &[u8], dk: &[u8]) -> Vec<u8> {
    let k = P::K::to_usize();

    // extract (from KEM decaps key) the PKE decryption key
    let dk_pke = &dk[0..ENCODE_12 * k];

    // extract PKE encryption key
    let ek_pke = &dk[ENCODE_12 * k..768 * k + 32];

    // extract hash of PKE encryption key
    let h = &dk[768 * k + 32..768 * k + 64];

    // extract implicit rejection value
    let z = &dk[768 * k + 64..768 * k + 96];

    // decrypt ciphertex
    let m_prime = k_pke_decrypt::<P>(dk_pke, c);

    // (K', r') ← G(m ∥ h)
    let mut hasher = Sha3_512::default();
    let mut m_prime_h = m_prime.to_vec();

    m_prime_h.extend_from_slice(h);
    hasher.update(&m_prime_h);
    let binding = hasher.finalize();
    let binding = binding.as_slice().to_vec();
    let (mut k_prime, r_prime) = binding.split_at(32);

    // K̄ ← J(z∥c, 32)
    let mut hasher = Sha3_512::default();
    z.to_vec().extend_from_slice(c);
    hasher.update(z);
    let binding = hasher.finalize();
    let k_bar = &binding.as_slice()[0..32];

    // re-encrypt using the derived randomness r′
    let c_prime = k_pke_encrypt::<P>(ek_pke, &m_prime, r_prime);
    if c != c_prime {
        k_prime = k_bar; // if ciphertexts do not match, “implicitly reject”
    }
    k_prime.to_vec()
}

// TODO: make this return an error and handle them internally
fn k_pke_decrypt<P: ParameterSet>(dk_pke: &[u8], c: &[u8]) -> Vec<u8> {
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
