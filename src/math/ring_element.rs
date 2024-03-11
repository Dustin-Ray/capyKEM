use core::ops::{Add, Sub};

use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::{
    constants::ml_kem_constants::{self, N},
    math::field_element::FieldElement as F,
};

/// A polynomial is an element of the ring R. It is an array of 256 coefficients
/// which themselves are [F].
#[derive(Clone, Copy, Debug)]
pub struct RingElement {
    pub coefficients: [F; 256],
}

impl RingElement {
    // Create a new RingElement from a vector of FieldElements
    pub fn new(val: [F; 256]) -> Self {
        RingElement { coefficients: val }
    }

    pub fn zero() -> Self {
        [F::new(0); 256].into()
    }

    pub fn byte_encode(self) -> Vec<u8> {
        let mut out = Vec::with_capacity(256 * 12 / 8); // Preallocate the output vector

        for i in (0..self.coefficients.len()).step_by(2) {
            // Combine two 12-bit integers into a single 24-bit integer
            let x = u32::from(self.coefficients[i].val())
                | (u32::from(self.coefficients[i + 1].val()) << 12);

            // Split the 24-bit integer into 3 bytes and append to the output vector
            out.push((x & 0xFF) as u8); // First 8 bits
            out.push(((x >> 8) & 0xFF) as u8); // Next 8 bits
            out.push(((x >> 16) & 0xFF) as u8); // Last 8 bits
        }
        out
    }

    fn byte_decode(b: &[u8]) -> Result<Vec<F>, &'static str> {
        if b.len() != ml_kem_constants::ENCODE_SIZE_12.into() {
            return Err("Invalid encoding length");
        }

        let mut f = Vec::with_capacity(ml_kem_constants::N.into());
        let mut i = 0;
        while i < b.len() {
            let d = u32::from(b[i]) | (u32::from(b[i + 1]) << 8) | (u32::from(b[i + 2]) << 16);
            const MASK_12: u32 = 0b1111_1111_1111;

            let elem1 = F::new((d & MASK_12) as u16)
                .check_reduced()
                .map_err(|_| "Invalid polynomial encoding")?;

            let elem2 = F::new((d >> 12) as u16)
                .check_reduced()
                .map_err(|_| "Invalid polynomial encoding")?;

            f.push(elem1);
            f.push(elem2);

            i += 3;
        }

        Ok(f)
    }

    pub fn sample_poly_cbd(s: &[u8], b: u8) -> RingElement {
        let mut prf = Shake256::default();
        prf.update(s);
        prf.update(&[b]);

        let mut b = [0u8; (N / 2) as usize];
        let mut reader = prf.finalize_xof();
        reader.read(&mut b);

        let mut f = [F::new(0); N as usize];

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
                f[i as usize] =
                    F::new((bits[0] + bits[1]).into()) - F::new((bits[2] + bits[3]).into());
            } else {
                f[i as usize] =
                    F::new((bits[4] + bits[5]).into()) - F::new((bits[2] + bits[3]).into());
            }
        }
        RingElement::new(f)
    }
}

// Implementing From<[F; 256]> for RingElement
impl From<[F; 256]> for RingElement {
    fn from(val: [F; 256]) -> Self {
        RingElement::new(val)
    }
}

impl Add for RingElement {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        assert_eq!(
            self.coefficients.len(),
            other.coefficients.len(),
            "RingElements must be of the same length"
        );

        let mut result = [F::default(); 256]; // Assuming F has a default value
        for (i, item) in self.coefficients.iter().enumerate().take(256) {
            result[i] = *item + other.coefficients[i];
        }
        RingElement::new(result)
    }
}

impl Sub for RingElement {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        assert_eq!(
            self.coefficients.len(),
            other.coefficients.len(),
            "RingElements must be of the same length"
        );

        let mut result = [F::default(); 256]; // Assuming F has a default value
        for (i, item) in self.coefficients.iter().enumerate().take(256) {
            result[i] = *item - other.coefficients[i];
        }
        RingElement::new(result)
    }
}

impl PartialEq for RingElement {
    fn eq(&self, other: &Self) -> bool {
        if self.coefficients.len() != other.coefficients.len() {
            return false;
        }
        self.coefficients
            .iter()
            .zip(other.coefficients.iter())
            .all(|(a, b)| a == b)
    }
}

#[cfg(test)]
mod tests {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    use super::*;

    // REMARKS:
    // axiom tests:
    // -[x] closure under addition
    // -[x] commutative
    // -[x] identity
    // -[ ] associativity
    // -[ ] distributivity

    #[test]
    fn test_additive_commutativity() {
        // Initialize a CRNG
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        // Step 1: Initialize a RingElement with random values
        let a = RingElement::sample_poly_cbd(&bytes, 0xAA);
        let b = RingElement::sample_poly_cbd(&bytes, 0xBB);

        // Perform the addition in both orders
        let a_plus_b = a + b;
        let b_plus_a = b + a;

        // Assert that a + b equals b + a
        assert_eq!(a_plus_b, b_plus_a, "Addition should be commutative.");
    }

    #[test]
    fn test_additive_closure() {
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        let a = RingElement::sample_poly_cbd(&bytes, 0xAA);
        let b = RingElement::sample_poly_cbd(&bytes, 0xBB);

        let result = a + b;

        assert!(
            result.coefficients.iter().all(|x| x.val() < 3329),
            "Summed elements are not reduced!"
        );
    }

    #[test]
    fn test_additive_identity() {
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        // Step 1: Initialize a RingElement with random values
        let a = RingElement::sample_poly_cbd(&bytes, 0xAA);
        let zero = RingElement::zero();

        // Adding the additive identity to a should result in a
        let a_plus_zero = a.clone() + zero;
        assert_eq!(
            a, a_plus_zero,
            "Adding the identity element should not change the element."
        );

        // Additionally, adding a to the additive identity should also result in a
        let zero_plus_a = zero + a.clone();
        assert_eq!(
            a, zero_plus_a,
            "Adding an element to the identity element should not change the element."
        );
    }

    #[test]
    fn test_additive_inverse() {
        // Create an example RingElement a
        let a_vals = [F::new(123); 256];
        let a = RingElement::new(a_vals);
        let mut inverse_vals = [F::default(); 256];
        for (i, val) in a.coefficients.iter().enumerate() {
            inverse_vals[i] = -*val; // Negate each coefficient
        }
        let b = RingElement::new(inverse_vals);

        // Perform the addition of a and its inverse
        let result = a + b;

        // The result should be the additive identity, i.e., all coefficients are 0
        // Verify that result is the additive identity of the ring
        let zero = RingElement::zero(); // The additive identity
        assert_eq!(result, zero, "a + additive inverse of a should be zero.");
    }
    #[test]
    fn test_encode_decode() {
        // Initialize a CRNG
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        // Initialize a RingElement with random values
        let a = RingElement::sample_poly_cbd(&bytes, 0xAA);

        // Encode the RingElement to bytes
        let encoded_bytes = a.byte_encode();

        // Decode the bytes back into a vector of F
        let decoded_elements = RingElement::byte_decode(&encoded_bytes).expect("Decoding failed");

        // Verify the decoded vector matches the original RingElement's array
        assert_eq!(
            decoded_elements.len(),
            a.coefficients.len(),
            "Length mismatch"
        );

        // Check each value for equality after encoding and decoding
        for (decoded_elem, original_elem) in decoded_elements.iter().zip(a.coefficients.iter()) {
            assert_eq!(decoded_elem.val(), original_elem.val(), "Value mismatch");
        }
    }
}
