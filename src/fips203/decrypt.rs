use super::encrypt::k_pke_encrypt;
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
use alloc::{string::String, vec::Vec};
use sha3::{Digest, Sha3_512};
use typenum::{Unsigned, U1};

/// FIPS 203 Section 6.3, Algorithm 17
/// Uses the decapsulation key to produce a shared key from a ciphertext.
pub fn mlkem_decaps<P: ParameterSet>(c: &[u8], dk: &[u8]) -> Result<Vec<u8>, String> {
    // Unpack the key based on parameter k
    let (dk_pke, ek_pke, h, z) = unpack_dk::<P>(dk);

    // Decrypt ciphertext
    let m_prime = k_pke_decrypt::<P>(dk_pke, c)?;

    // Derive K' and r' from m_prime and h
    let (mut k_prime, r_prime) = derive_keys(&m_prime, h);

    // Compute K̄ from z and c
    let k_bar = compute_k_bar(z, c);

    // Re-encrypt using derived randomness r' and check ciphertext match
    let c_prime = k_pke_encrypt::<P>(ek_pke, &m_prime, &r_prime)?;
    if c != c_prime {
        k_prime = k_bar; // If ciphertexts do not match, "implicitly reject"
    }

    Ok(k_prime)
}

// Extracts keys from dk based on the size multiplier k
fn unpack_dk<P: ParameterSet>(dk: &[u8]) -> (&[u8], &[u8], &[u8], &[u8]) {
    let k = P::K::to_usize();
    let dk_pke = &dk[0..ENCODE_12 * k];
    let ek_pke = &dk[ENCODE_12 * k..768 * k + 32];
    let h = &dk[768 * k + 32..768 * k + 64];
    let z = &dk[768 * k + 64..768 * k + 96];
    (dk_pke, ek_pke, h, z)
}

// Derive K' and r' using SHA3-512 hasher
fn derive_keys(m_prime: &[u8], h: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut hasher = Sha3_512::default();
    hasher.update(m_prime);
    hasher.update(h);
    let binding = hasher.finalize().as_slice().to_vec();
    let (k_prime, r_prime) = binding.split_at(32);
    (k_prime.to_vec(), r_prime.to_vec())
}

// Compute K̄ using SHA3-512 hasher
fn compute_k_bar(z: &[u8], c: &[u8]) -> Vec<u8> {
    let mut hasher = Sha3_512::default();
    hasher.update(z);
    hasher.update(c);
    hasher.finalize().as_slice()[0..32].to_vec()
}

// FIPS 203 Section 5.3 Algorithm 14
// Uses the decryption key to decrypt a ciphertext.
fn k_pke_decrypt<P: ParameterSet>(dk_pke: &[u8], c: &[u8]) -> Result<Vec<u8>, String> {
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
        let f = NttElement::byte_decode_12(current)?;
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
    Ok(s)
}
