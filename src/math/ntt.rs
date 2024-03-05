use hex_literal::hex;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::constants::ml_kem_constants::N;

use super::{field_element::FieldElement, ring_element::RingElement};

fn sample_poly_cbd(s: &[u8], b: u8) -> RingElement {
    let mut prf = Shake256::default();
    prf.update(s);
    prf.update(&[b]);

    let mut b = [0u8; (N / 2) as usize];
    let mut reader = prf.finalize_xof();
    reader.read(&mut b);

    // dbg!(b.clone());

    let mut f = [FieldElement::new(0); N as usize];
    for i in 0..N {
        let b = b[(i / 2) as usize];
        let bits = [
            (b >> 7) & 1,
            (b >> 6) & 1,
            (b >> 5) & 1,
            (b >> 4) & 1,
            (b >> 3) & 1,
            (b >> 2) & 1,
            (b >> 1) & 1,
            b & 1,
        ];
        
        // The i-th coefficient is based on the first four bits
        // The (i+1)-th coefficient is based on the second four bits
        if i % 2 == 0 {
            f[i as usize] = FieldElement::new((bits[0] + bits[1]).into())
                - FieldElement::new((bits[2] + bits[3]).into());
        } else {
            f[i as usize] = FieldElement::new((bits[4] + bits[5]).into())
                - FieldElement::new((bits[2] + bits[3]).into());
        }
    }
    // dbg!(f.clone());
    RingElement::new(f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_poly_cbd() {
        // Define a known seed and byte to initialize the PRF
        let seed: &[u8] = &[0x12, 0x34, 0x56, 0x78]; // Example seed
        let byte: u8 = 0xAB; // Example additional byte

        // Call the function to test
        let result = sample_poly_cbd(seed, byte);

        // Check the size of the resulting RingElement
        assert_eq!(result.val.len(), N.into());

        // Perform additional checks
        // Since the function is random, we can't check for specific values without a fixed PRF output.
        // However, we can check properties that should always hold, such as the range of values.
        for coef in result.val.iter() {
            assert!(coef.val >= 0 && coef.val < 3329); // Assuming the coefficients are mod 3329 and FieldElement has a tuple struct with the value as the first element.
        }

        // To make a more precise test, you might want to mock the PRF to produce a predictable output.
        // This would involve creating a version of the function where you can inject the PRF or its output.
    }
}