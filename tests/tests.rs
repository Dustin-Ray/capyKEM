#[cfg(test)]
mod tests {

    use capy_kem::{
        constants::parameter_sets::{KEM_1024, KEM_512, KEM_768},
        fips203::{decrypt::mlkem_decaps, encrypt::mlkem_encaps, keygen::ml_kem_keygen},
    };
    use rand::thread_rng;

    #[test]
    #[allow(non_snake_case)]
    fn roundtrip_768() {
        let mut rng = thread_rng();
        let (ek_pke, dk_pke) = ml_kem_keygen::<KEM_768, _>(&mut rng); // Generate key pair

        let (K, c) = mlkem_encaps::<KEM_768, _>(&ek_pke.ek, &mut rng).unwrap();

        let dec = mlkem_decaps::<KEM_768>(&c, &dk_pke.dk).unwrap();
        // Assert that the decrypted message matches the original message
        assert_eq!(dec, K);
    }

    #[test]
    #[allow(non_snake_case)]
    fn roundtrip_512() {
        let mut rng = thread_rng();
        let (ek_pke, dk_pke) = ml_kem_keygen::<KEM_512, _>(&mut rng); // Generate key pair

        let (K, c) = mlkem_encaps::<KEM_512, _>(&ek_pke.ek, &mut rng).unwrap();

        let dec = mlkem_decaps::<KEM_512>(&c, &dk_pke.dk).unwrap();
        // Assert that the decrypted message matches the original message
        assert_eq!(dec, K);
    }

    #[test]
    #[allow(non_snake_case)]
    fn roundtrip_1024() {
        let mut rng = thread_rng();
        let (ek_pke, dk_pke) = ml_kem_keygen::<KEM_1024, _>(&mut rng); // Generate key pair

        let (K, c) = mlkem_encaps::<KEM_1024, _>(&ek_pke.ek, &mut rng).unwrap();

        let dec = mlkem_decaps::<KEM_1024>(&c, &dk_pke.dk).unwrap();
        // Assert that the decrypted message matches the original message
        assert_eq!(dec, K);
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_private_key_redaction() {
        let mut rng = thread_rng();
        let (_, dk) = ml_kem_keygen::<KEM_768, _>(&mut rng);

        // Verify that Debug output doesn't leak secret material
        let debug_output = format!("{:?}", dk);
        assert!(debug_output.contains("<redacted>"));
        assert!(!debug_output.contains(&format!("{:?}", dk.dk)));

        // Verify that Display output doesn't leak secret material
        let display_output = format!("{}", dk);
        assert!(display_output.contains("<redacted>"));
    }
}
