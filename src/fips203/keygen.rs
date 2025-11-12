use crate::{
    constants::{ml_kem_constants::ENCODE_12, parameter_sets::ParameterSet},
    math::{ntt_element::NttElement, ring_element::RingElement},
};
use alloc::vec::Vec;
use core::fmt;
use hybrid_array::{typenum::Unsigned, Array};
use rand_core::{CryptoRng, RngCore};
use serde::{Deserialize, Serialize};
use sha3::{Digest, Sha3_512};
use zeroize::{Zeroize, ZeroizeOnDrop};

/// Represents a private key for Key Encapsulation Mechanism (KEM).
///
/// This structure holds the private decryption key (`dk`)
/// necessary for the decryption process in KEM.
///
/// ## Fields
///
/// * `dk: Vec<u8>` - The private decryption key data,
///   essential for decrypting the KEM ciphertext.
///
/// ## Security
///
/// This type implements `ZeroizeOnDrop` to ensure that private key material
/// is securely erased from memory when dropped. The Debug and Display
/// implementations are redacted to prevent accidental leakage of secret material.
#[derive(Serialize, Deserialize, Clone, Zeroize, ZeroizeOnDrop)]
pub struct KEMPrivateKey {
    pub dk: Vec<u8>,
}

impl fmt::Debug for KEMPrivateKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("KEMPrivateKey")
            .field("dk", &"<redacted>")
            .finish()
    }
}

impl fmt::Display for KEMPrivateKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "KEMPrivateKey {{ dk: <redacted> }}")
    }
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
///
/// * `rand_bytes: [u8; 32]` - Random bytes used to seed KEM
///   operations, ensuring the uniqueness and security of the public key.
/// * `ek: Vec<u8>` - The public encryption key data,
///   used to encrypt data in the KEM scheme.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct KEMPublicKey {
    pub ek: Vec<u8>,
}

/// Generates a public-private key pair for use with the Key Encapsulation Mechanism (KEM).
///
/// This function generates a ML-KEM key pair for the specified parameter set.
/// All three parameter sets (ML-KEM-512, ML-KEM-768, ML-KEM-1024) are supported.
///
/// # Arguments
///
/// * `rng` - A cryptographically secure random number generator
///
/// # Security
///
/// This function zeroizes sensitive intermediate values. The RNG must implement
/// `CryptoRng` for security.
///
/// # Examples
///
/// ```ignore
/// use rand::thread_rng;
/// use capy_kem::fips203::keygen::ml_kem_keygen;
/// use capy_kem::constants::parameter_sets::KEM_768;
///
/// let mut rng = thread_rng();
/// let (pk, sk) = ml_kem_keygen::<KEM_768>(&mut rng);
/// ```
pub fn ml_kem_keygen<P: ParameterSet, R: RngCore + CryptoRng>(
    rng: &mut R,
) -> (KEMPublicKey, KEMPrivateKey) {
    let mut z = [0u8; 32];

    // Generate randomness for the KEM
    rng.fill_bytes(&mut z);
    let (ek, mut dk) = k_pke_keygen::<P>(&z);

    let h_ek = hash_ek(&ek);

    // Concatenate dk, ek, h_ek, and z into a single Vec<u8>
    pack_dk(&mut dk, &ek, &h_ek, &z);

    // Zeroize sensitive intermediate value
    z.zeroize();

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

fn k_pke_keygen<P: ParameterSet>(d: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
    let k = P::K::to_usize();
    let mut hasher = Sha3_512::default();
    hasher.update(d);
    let binding = hasher.finalize();
    let b = binding.as_slice();

    // (ρ, σ ) <- G(d)
    let rho: &[u8] = &b[0..32];
    let sigma = &b[32..64];

    let mut n = 0;

    // Generate the matrix a_hat (k * k elements)
    let mut a_hat = Array::<NttElement, P::KSquared>::default();
    for i in 0..k {
        for j in 0..k {
            // see: https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/s-C-zIAeKfE/m/eZJmXYsSAQAJ?
            a_hat[i * k + j] = NttElement::sample_ntt(rho, j, i);
        }
    }

    // generate s (uses EtaOne)
    let mut s_hat = Array::<NttElement, P::K>::default();
    for s_elem in s_hat.iter_mut() {
        *s_elem = RingElement::sample_poly_cbd::<P::EtaOne>(sigma, n).into();
        n += 1;
    }

    // generate e (uses EtaOne)
    let mut e_hat = Array::<NttElement, P::K>::default();
    for e_elem in e_hat.iter_mut() {
        *e_elem = RingElement::sample_poly_cbd::<P::EtaOne>(sigma, n).into();
        n += 1;
    }

    // t_hat = A o s_hat + e_hat
    let mut t = Array::<NttElement, P::K>::default();
    for i in 0..t.len() {
        t[i] = e_hat[i];
        for j in 0..s_hat.len() {
            t[i] += a_hat[i * k + j] * s_hat[j];
        }
    }

    // ByteEncode12(t_hat||rho)
    let ek_pke_size = ENCODE_12 * k + 32;
    let mut ek_pke: Vec<u8> = Vec::with_capacity(ek_pke_size);
    for item in t.iter() {
        ek_pke = item.byte_encode_12(ek_pke);
    }
    ek_pke.extend_from_slice(rho);

    let dk_pke_size = ENCODE_12 * k;
    let mut dk_pke: Vec<u8> = Vec::with_capacity(dk_pke_size);
    for item in s_hat.iter() {
        dk_pke = item.byte_encode_12(dk_pke);
    }
    (ek_pke, dk_pke)
}
