#[cfg(test)]
mod tests {
    use rand::{thread_rng, RngCore};
    use capykem::{constants::parameter_sets::KEM_768, Decrypt, Encrypt, Secret};

    #[test]
    fn roundtrip() {
        let s = Secret::<KEM_768>::new([4_u8; 32]); // Create a new Secret with the KEM_768 parameter

        let mut rng = thread_rng(); // Get a thread-local RNG
        let mut rand_bytes = [0u8; 32]; // Create a buffer for the randomness
        rng.fill_bytes(&mut rand_bytes); // Generate secure random bytes

        let (ek, dk_pke) = s.k_pke_keygen(&rand_bytes); // Generate key pair

        // Encrypt using the encrypt method from the Encrypt trait
        let c = s.encrypt(&ek, &rand_bytes);

        // Decrypt using the decrypt method from the Decrypt trait
        let dec = s.decrypt(&dk_pke, &c);

        // Assert that the decrypted message matches the original message
        assert_eq!(dec, s.s);
    }
}