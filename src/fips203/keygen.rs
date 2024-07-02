use crate::{
    constants::{
        ml_kem_constants::{k, D_PKE_KEYSIZE, E_PKE_KEYSIZE},
        parameter_sets::ParameterSet,
    },
    math::{ntt_element::NttElement, ring_element::RingElement},
};
use alloc::vec::Vec;
use serde::{Deserialize, Serialize};
use sha3::{Digest, Sha3_512};
use rand::{thread_rng, RngCore};

/// Represents a private key for Key Encapsulation Mechanism (KEM).
///
/// This structure holds the private decryption key (`dk`)
/// necessary for the decryption process in KEM.
///
/// ## Fields
/// * `dk: Vec<u8>` - The private decryption key data,
/// essential for decrypting the KEM ciphertext.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct KEMPrivateKey {
    pub dk: Vec<u8>,
}

/// Represents a public key for Key Encapsulation Mechanism (KEM).
///
/// This structure holds the public encryption key (`ek`) and a set of random bytes
/// (`rand_bytes`) used to initialize or seed certain operations within the KEM.
/// The public key is used in the encryption process,
/// encapsulating data in such a way that only someone with the
/// corresponding private key can decrypt it.
///
/// ## Fields
/// * `rand_bytes: [u8; 32]` - Random bytes used to seed KEM
/// operations, ensuring the uniqueness and security of the public key.
/// * `ek: Vec<u8>` - The public encryption key data,
/// used to encrypt data in the KEM scheme.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct KEMPublicKey {
    pub ek: Vec<u8>,
}

/// Generates a public-private key pair for use with the Key Encapsulation Mechanism (KEM).
///
/// This function interfaces with [`capy_kem`] for partial ML-KEM-768 support.(Partial because
/// the other parameter sets are a work in progress) to generate a compatible key pair.
///
/// It initializes the necessary randomness and calls the library-specific key generation
/// function to produce both encryption and decryption keys.
pub fn ml_kem_keygen<P: ParameterSet>() -> (KEMPublicKey, KEMPrivateKey) {
    let mut rng = thread_rng();
    let mut z = [0u8; 32];

    // Generate randomness for the KEM
    rng.fill_bytes(&mut z);
    let (ek, mut dk) = k_pke_keygen::<P>(&z);

    let h_ek = hash_ek(&ek);

    // Concatenate dk, ek, h_ek, and z into a single Vec<u8>
    pack_dk(&mut dk, &ek, &h_ek, &z);

    (KEMPublicKey { ek }, KEMPrivateKey { dk })
}

/// Hashes the encryption key and returns the first 32 bytes of the hash
fn hash_ek(ek: &[u8]) -> Vec<u8> {
    let mut hasher = Sha3_512::default();
    hasher.update(ek);
    hasher.finalize().as_slice()[0..32].to_vec()
}

/// Concatenates dk, ek, h_ek, and z into dk
fn pack_dk(dk: &mut Vec<u8>, ek: &[u8], h_ek: &[u8], z: &[u8]) {
    dk.extend_from_slice(ek);
    dk.extend_from_slice(h_ek);
    dk.extend_from_slice(z);
}

// TODO: parameterize this
#[allow(clippy::extra_unused_type_parameters)]
fn k_pke_keygen<P: ParameterSet>(d: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
    let mut hasher = Sha3_512::default();
    hasher.update(d);
    let binding = hasher.finalize();
    let b = binding.as_slice();

    // (ρ, σ ) <- G(d)
    let rho: &[u8] = &b[0..32];
    let sigma = &b[32..64];

    let mut n = 0;

    // Generate the matrix a_hat
    let mut a_hat: [NttElement; 9] = [NttElement::zero(); k * k];
    for i in 0..k {
        for j in 0..k {
            // see: https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/s-C-zIAeKfE/m/eZJmXYsSAQAJ?
            a_hat[i * k + j] = NttElement::sample_ntt(rho, j, i);
        }
    }

    // generate s
    let mut s_hat = [NttElement::zero(); k];
    for s_elem in s_hat.iter_mut().take(k) {
        *s_elem = RingElement::sample_poly_cbd(sigma, n).into();
        n += 1;
    }

    // generate e
    let mut e_hat = [NttElement::zero(); k];
    for e_elem in e_hat.iter_mut().take(k) {
        *e_elem = RingElement::sample_poly_cbd(sigma, n).into();
        n += 1;
    }

    // t_hat = A o s_hat + e_hat
    let mut t = [NttElement::zero(); k];
    for i in 0..t.len() {
        t[i] = e_hat[i];
        for j in 0..s_hat.len() {
            t[i] += a_hat[i * k + j] * s_hat[j];
        }
    }

    // ByteEncode12(t_hat||rho)
    let mut ek_pke: Vec<u8> = Vec::with_capacity(E_PKE_KEYSIZE);
    for &item in t.iter() {
        ek_pke = item.byte_encode_12(ek_pke);
    }
    ek_pke.append(&mut rho.into());

    let mut dk_pke: Vec<u8> = Vec::with_capacity(D_PKE_KEYSIZE);
    for &item in s_hat.iter() {
        dk_pke = item.byte_encode_12(dk_pke);
    }
    // pretty_print_vec_u8(&dk_pke);
    (ek_pke, dk_pke)
}
