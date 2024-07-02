#[cfg(test)]
mod tests {

    use capy_kem::{
        constants::parameter_sets::KEM_768,
        fips203::{decrypt::mlkem_decaps, encrypt::mlkem_encaps, keygen::ml_kem_keygen},
    };

    #[test]
    #[allow(non_snake_case)]
    fn roundtrip() {
        let (ek_pke, dk_pke) = ml_kem_keygen::<KEM_768>(); // Generate key pair

        let (K, c) = mlkem_encaps::<KEM_768>(&ek_pke.ek).unwrap();

        let dec = mlkem_decaps::<KEM_768>(&c, &dk_pke.dk);
        // Assert that the decrypted message matches the original message
        assert_eq!(dec, K);
    }
}
