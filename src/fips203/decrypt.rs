use super::encrypt::k_pke_encrypt;
use crate::{
    constants::{ml_kem_constants::ENCODE_12, parameter_sets::ParameterSet},
    error::Result,
    math::{
        encoding::{Compress, Encode, EncodingSize},
        ntt_element::NttElement,
        ring_element::RingElement,
    },
};
use alloc::vec::Vec;
use hybrid_array::{typenum::Unsigned, Array};
use sha3::{Digest, Sha3_512};
use subtle::ConstantTimeEq;
use typenum::U1;
use zeroize::Zeroize;

/// FIPS 203 Section 6.3, Algorithm 17
/// Uses the decapsulation key to produce a shared key from a ciphertext.
///
/// # Security
///
/// This function uses constant-time comparison to prevent timing attacks
/// and zeroizes sensitive intermediate values.
pub fn mlkem_decaps<P: ParameterSet>(c: &[u8], dk: &[u8]) -> Result<Vec<u8>> {
    // Unpack the key based on parameter k
    let (dk_pke, ek_pke, h, z) = unpack_dk::<P>(dk);

    // Decrypt ciphertext
    let mut m_prime = k_pke_decrypt::<P>(dk_pke, c)?;

    // Derive K' and r' from m_prime and h
    let (mut k_prime, mut r_prime) = derive_keys(&m_prime, h);

    // Compute K̄ from z and c
    let k_bar = compute_k_bar(z, c);

    // Re-encrypt using derived randomness r' and check ciphertext match
    let c_prime = k_pke_encrypt::<P>(ek_pke, &m_prime, &r_prime)?;
    
    // Constant-time comparison to prevent timing attacks
    let comparison = c.ct_eq(&c_prime);
    if comparison.unwrap_u8() != 1 {
        k_prime = k_bar; // If ciphertexts do not match, "implicitly reject"
    }

    // Zeroize sensitive intermediate values
    m_prime.zeroize();
    r_prime.zeroize();

    Ok(k_prime)
}

// Extracts keys from dk based on the size multiplier k
fn unpack_dk<P: ParameterSet>(dk: &[u8]) -> (&[u8], &[u8], &[u8], &[u8]) {
    let k = P::K::to_usize();
    let dk_pke_size = ENCODE_12 * k;
    let ek_pke_size = ENCODE_12 * k + 32;
    let dk_pke = &dk[0..dk_pke_size];
    let ek_pke = &dk[dk_pke_size..dk_pke_size + ek_pke_size];
    let h = &dk[dk_pke_size + ek_pke_size..dk_pke_size + ek_pke_size + 32];
    let z = &dk[dk_pke_size + ek_pke_size + 32..dk_pke_size + ek_pke_size + 64];
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
fn k_pke_decrypt<P: ParameterSet>(dk_pke: &[u8], c: &[u8]) -> Result<Vec<u8>> {
    let encode_du_size = <<P as ParameterSet>::Du as EncodingSize>::EncodedPolynomialSize::USIZE;
    let mut slice = c;
    let mut u = Array::<RingElement, P::K>::default();
    for i in 0..P::K::to_usize() {
        let (current, next) = slice.split_at(encode_du_size);
        let mut f: RingElement = Encode::<P::Du>::decode(current);
        f.decompress::<P::Du>();

        u[i] = f;
        slice = next;
    }

    let mut slice = dk_pke;
    let mut s_hat = Array::<NttElement, P::K>::default();
    for i in 0..P::K::to_usize() {
        let (current, next) = slice.split_at(ENCODE_12);
        let f = NttElement::byte_decode_12(current)?;
        s_hat[i] = f;
        slice = next;
    }

    let c2_size = <<P as ParameterSet>::Dv as EncodingSize>::EncodedPolynomialSize::USIZE;
    let mut v: RingElement = Encode::<P::Dv>::decode(&c[c.len() - c2_size..c.len()]);
    v.decompress::<P::Dv>();

    let mut y = RingElement::zero();
    for i in 0..s_hat.len() {
        y += (s_hat[i] * u[i].into()).into();
    }

    let mut w = v - y;
    let s = Encode::<U1>::encode(w.compress::<U1>());
    Ok(s)
}
