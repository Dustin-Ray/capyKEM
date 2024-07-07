use crate::{
    constants::{
        ml_kem_constants::{k, ENCODE_12},
        parameter_sets::ParameterSet,
    },
    math::{
        encoding::{byte_decode, byte_encode, Compress, Encode},
        ntt_element::NttElement,
        ring_element::RingElement,
    },
};
use alloc::{string::String, vec::Vec};
use rand::{thread_rng, RngCore};
use sha3::{Digest, Sha3_512};
use typenum::U1;

#[allow(non_snake_case)]
pub fn mlkem_encaps<P: ParameterSet>(ek: &[u8]) -> Result<(Vec<u8>, Vec<u8>), String> {
    // Step 1. (Type check) Validate the key length
    if ek.len() != ENCODE_12 * k + 32 {
        return Err("Key length validation failed".into());
    }

    // Step 2. modulus check ek~ <- ByteEncode12(ByteDecode12(ek))
    // TODO: this is wrong, the entire things should be equal, not portions
    // TODO: make this fixed-time
    let ek_decoded = byte_decode::<P::Encode12>(ek);
    let ek_encoded = byte_encode::<P::Encode12>(&ek_decoded.coefs);
    assert_eq!(ek[0..384], ek_encoded);

    if ek_encoded != ek[0..384] {
        // Compare byte-wise
        return Err("Modulus check failed: Key is not consistent after encode-decode cycle".into());
    }

    // Step 3. Generate 32 random bytes (see Section 3.3)
    let mut rng = thread_rng();
    let mut m = [0_u8; 32];
    rng.fill_bytes(&mut m);

    // Step 4. Compute hash of encryption key
    let h_ek = hash_to_slice(ek, 32);

    // Step 5. Concatenate m and h_ek, and hash to derive K and r
    let (K, r) = derive_keys(&m, &h_ek);

    // Step 6. Encrypt the message
    let c = k_pke_encrypt::<P>(ek, &m, &r)?;

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
) -> Result<Vec<u8>, String> {
    let mut n = 0;
    let mut t_hat = [NttElement::zero(); k];

    for i in 0..t_hat.len() {
        t_hat[i] = NttElement::byte_decode_12(&ek_pke[i * ENCODE_12..(i + 1) * ENCODE_12])?;
    }

    let rho: &[u8] = &ek_pke[ENCODE_12 * k..(ENCODE_12 * k) + 32];

    // Generate the matrix a_hat^T
    let mut a_hat_transpose = [NttElement::zero(); k * k];
    for i in 0..k {
        for j in 0..k {
            a_hat_transpose[i * k + j] = NttElement::sample_ntt(rho, i, j);
        }
    }

    // generate r, run ntt k times
    let mut r_hat = [NttElement::zero(); k];
    for r_elem in r_hat.iter_mut().take(k) {
        *r_elem = RingElement::sample_poly_cbd(rand, n).into();
        n += 1;
    }

    // generate e1
    let mut e_1 = [RingElement::zero(); k];
    for e_elem in e_1.iter_mut().take(k) {
        *e_elem = RingElement::sample_poly_cbd(rand, n);
        n += 1;
    }

    // sample e2
    let e2: RingElement = RingElement::sample_poly_cbd(rand, n);

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
