use crate::{
    constants::{ml_kem_constants::ENCODE_12, parameter_sets::ParameterSet},
    error::{KemError, Result},
    math::{
        encoding::{Compress, Encode},
        ntt_element::NttElement,
        ring_element::RingElement,
    },
};
use alloc::vec::Vec;
use hybrid_array::{typenum::Unsigned, Array};
use rand_core::{CryptoRng, RngCore};
use sha3::{Digest, Sha3_512};
use subtle::ConstantTimeEq;
use typenum::U1;
use zeroize::Zeroize;

/// Encapsulation with provided RNG
///
/// # Security
/// 
/// This function performs constant-time comparisons and zeroizes sensitive
/// intermediate values. The RNG must implement `CryptoRng` for security.
#[allow(non_snake_case)]
pub fn mlkem_encaps<P: ParameterSet, R: RngCore + CryptoRng>(
    ek: &[u8],
    rng: &mut R,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let k = P::K::to_usize();
    let ek_pke_size = ENCODE_12 * k;

    // Step 1. (Type check) Validate the key length
    if ek.len() != ek_pke_size + 32 {
        return Err(KemError::InvalidInput);
    }

    // Step 2. modulus check ek~ <- ByteEncode12(ByteDecode12(ek))
    // Using constant-time comparison to prevent timing attacks
    let mut ek_reencoded = Vec::with_capacity(ek_pke_size);
    for i in 0..k {
        let poly_slice = &ek[i * ENCODE_12..(i + 1) * ENCODE_12];
        let decoded = NttElement::byte_decode_12(poly_slice)?;
        ek_reencoded = decoded.byte_encode_12(ek_reencoded);
    }

    // Constant-time comparison
    let comparison = ek_reencoded.ct_eq(&ek[0..ek_pke_size]);
    if comparison.unwrap_u8() != 1 {
        // Zeroize before returning error
        ek_reencoded.zeroize();
        return Err(KemError::InvalidInput);
    }
    ek_reencoded.zeroize();

    // Step 3. Generate 32 random bytes (see Section 3.3)
    let mut m = [0_u8; 32];
    rng.fill_bytes(&mut m);

    // Step 4. Compute hash of encryption key
    let h_ek = hash_to_slice(ek, 32);

    // Step 5. Concatenate m and h_ek, and hash to derive K and r
    let (K, mut r) = derive_keys(&m, &h_ek);

    // Step 6. Encrypt the message
    let c = k_pke_encrypt::<P>(ek, &m, &r)?;

    // Zeroize sensitive intermediate values
    m.zeroize();
    r.zeroize();

    Ok((K, c))
}

fn hash_to_slice(data: &[u8], slice_size: usize) -> Vec<u8> {
    let mut hasher = Sha3_512::default();
    hasher.update(data);
    hasher.finalize().as_slice()[0..slice_size].to_vec()
}

#[allow(non_snake_case)]
fn derive_keys(m: &[u8; 32], h_ek: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut hasher = Sha3_512::default();
    hasher.update(m);
    hasher.update(h_ek);
    let binding = hasher.finalize();
    let (K, r) = binding.as_slice().split_at(32);
    (K.to_vec(), r.to_vec())
}

pub(crate) fn k_pke_encrypt<P: ParameterSet>(
    ek_pke: &[u8],
    m: &[u8],
    rand: &[u8],
) -> Result<Vec<u8>> {
    let k = P::K::to_usize();
    let mut n = 0;
    let mut t_hat = Array::<NttElement, P::K>::default();

    for i in 0..t_hat.len() {
        t_hat[i] = NttElement::byte_decode_12(&ek_pke[i * ENCODE_12..(i + 1) * ENCODE_12])?;
    }

    let rho: &[u8] = &ek_pke[ENCODE_12 * k..(ENCODE_12 * k) + 32];

    // Generate the matrix a_hat^T
    let mut a_hat_transpose = Array::<NttElement, P::KSquared>::default();
    for i in 0..k {
        for j in 0..k {
            a_hat_transpose[i * k + j] = NttElement::sample_ntt(rho, i, j);
        }
    }

    // generate r, run ntt k times (uses EtaTwo)
    let mut r_hat = Array::<NttElement, P::K>::default();
    for r_elem in r_hat.iter_mut() {
        *r_elem = RingElement::sample_poly_cbd::<P::EtaTwo>(rand, n).into();
        n += 1;
    }

    // generate e1 (uses EtaTwo)
    let mut e_1 = Array::<RingElement, P::K>::default();
    for e_elem in e_1.iter_mut() {
        *e_elem = RingElement::sample_poly_cbd::<P::EtaTwo>(rand, n);
        n += 1;
    }

    // sample e2 (uses EtaTwo)
    let e2: RingElement = RingElement::sample_poly_cbd::<P::EtaTwo>(rand, n);

    let mut u: Vec<RingElement> = e_1
        .iter()
        .enumerate()
        .map(|(i, e1_elem)| {
            let sum: RingElement = (0..k)
                .map(|j| (a_hat_transpose[i * k + j] * r_hat[j]).into())
                .sum();
            *e1_elem + sum
        })
        .collect();

    let mut mu: RingElement = Encode::<U1>::decode(m);
    mu.decompress::<U1>();

    let mut v = NttElement::zero();
    for i in 0..t_hat.len() {
        v += t_hat[i] * r_hat[i];
    }
    let mut v = v.ntt_inv();
    v += e2;
    v += mu;

    let mut c: Vec<u8> = Vec::new();

    for ring in u.iter_mut() {
        let bytes = &mut Encode::<P::Du>::encode(ring.compress::<P::Du>());
        c.append(bytes);
    }

    c.append(&mut Encode::<P::Dv>::encode(v.compress::<P::Dv>()));

    Ok(c)
}
