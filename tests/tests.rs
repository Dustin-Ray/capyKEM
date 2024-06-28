#[cfg(test)]
mod tests {

    use capy_kem::{
        constants::parameter_sets::KEM_768,
        fips203::{decrypt::k_pke_decrypt, encrypt::k_pke_encrypt, keygen::k_pke_keygen},
    };
    use rand::{thread_rng, RngCore};

    #[test]
    fn roundtrip() {
        let mut secret = [0u8; 32]; // Generate a random secret
        let mut rand_bytes = [0u8; 32]; // Generate randomness for the KEM

        let mut rng = thread_rng(); // Get a thread-local RNG
        rng.fill_bytes(&mut rand_bytes); // Generate secure random bytes
        rng.fill_bytes(&mut secret); // Generate secure random bytes

        let (ek, dk_pke) = k_pke_keygen::<KEM_768>(&rand_bytes); // Generate key pair

        // Encrypt using the encrypt method from the Encrypt trait
        let c = k_pke_encrypt::<KEM_768>(&secret, &ek, &rand_bytes);

        // Decrypt using the decrypt method from the Decrypt trait
        let dec = k_pke_decrypt::<KEM_768>(&dk_pke, &c);
        // Assert that the decrypted message matches the original message
        assert_eq!(dec, secret);
    }
}
