use core::fmt;
use core::ops::AddAssign;
use core::ops::{Add, Sub};

use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use alloc::vec::Vec;

use crate::constants::ml_kem_constants::{self, N};
use crate::constants::parameter_sets::ParameterSet;
use crate::math::field_element::FieldElement as F;

/// A polynomial is an element of the ring R. It is an array of 256 coefficients
/// which themselves are [F].
#[derive(Clone, Copy)]
pub struct RingElement<P> {
    pub coefficients: [F<P>; 256],
}

impl<P: ParameterSet + Copy> RingElement<P> {
    // Create a new RingElement from a vector of FieldElements
    pub fn new(val: [F<P>; 256]) -> Self {
        RingElement { coefficients: val }
    }

    pub fn zero() -> Self {
        [F::new(0); 256].into()
    }

    // REMARKS:
    // TODO: parameterize du and dv
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

    pub fn byte_decode(b: &[u8]) -> Result<Self, &'static str> {
        if b.len() != ml_kem_constants::ENCODE_SIZE_12.into() {
            return Err("Invalid encoding length");
        }

        let mut f = Vec::with_capacity(N.into());
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
        let array: [F<P>; 256] = f
            .try_into()
            .map_err(|_| "Conversion to fixed-size array failed")?;

        Ok(RingElement::new(array))
    }

    // REMARKS:
    // TODO: parameterize eta1 and eta 2
    pub fn sample_poly_cbd(s: &[u8], b: u8) -> RingElement<P> {
        let mut prf = Shake256::default();
        prf.update(s);
        prf.update(&[b]);

        let mut b = [0u8; (N / 2) as usize];
        let mut reader = prf.finalize_xof();
        reader.read(&mut b);

        let mut f = [F::new(0); N as usize];

        for i in (0..N).step_by(2) {
            // Iterate through indices, stepping by 2.
            let b = b[(i / 2) as usize];
            let b_7 = (b >> 7) & 1;
            let b_6 = (b >> 6) & 1;
            let b_5 = (b >> 5) & 1;
            let b_4 = (b >> 4) & 1;
            let b_3 = (b >> 3) & 1;
            let b_2 = (b >> 2) & 1;
            let b_1 = (b >> 1) & 1;
            let b_0 = b & 1;

            f[i as usize] = F::new((b_0 + b_1).into()) - F::new((b_2 + b_3).into());
            // Ensure i+1 doesn't go out of bounds, relevant if N is odd.
            if i + 1 < N {
                f[(i + 1) as usize] = F::new((b_4 + b_5).into()) - F::new((b_6 + b_7).into());
            }
        }
        RingElement::new(f)
    }
}

impl<P: ParameterSet + Copy> fmt::Debug for RingElement<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (index, element) in self.coefficients.iter().enumerate() {
            write!(f, "{:<8}", element.val())?;
            // Adjust for row width
            if (index + 1) % 16 == 0 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

impl<P: ParameterSet + Copy> From<[F<P>; 256]> for RingElement<P> {
    fn from(val: [F<P>; 256]) -> Self {
        RingElement::new(val)
    }
}

impl<P: ParameterSet + Copy> From<&[u8]> for RingElement<P> {
    fn from(slice: &[u8]) -> Self {
        RingElement::byte_decode(slice).expect("failed to decode bytes, check reduction?")
    }
}

impl<P: ParameterSet + Copy> AddAssign for RingElement<P> {
    fn add_assign(&mut self, other: Self) {
        for (lhs, rhs) in self.coefficients.iter_mut().zip(other.coefficients.iter()) {
            *lhs += *rhs;
        }
    }
}

impl<P: ParameterSet + Copy> Add for RingElement<P> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        assert_eq!(
            self.coefficients.len(),
            other.coefficients.len(),
            "RingElements must be of the same length"
        );

        let mut result = [F::zero(); 256];
        for (i, item) in self.coefficients.iter().enumerate().take(256) {
            result[i] = *item + other.coefficients[i];
        }
        RingElement::new(result)
    }
}

impl<P: ParameterSet + Copy> Sub for RingElement<P> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        assert_eq!(
            self.coefficients.len(),
            other.coefficients.len(),
            "RingElements must be of the same length"
        );

        let mut result = [F::zero(); 256];
        for (i, item) in self.coefficients.iter().enumerate().take(256) {
            result[i] = *item - other.coefficients[i];
        }
        RingElement::new(result)
    }
}

impl<P: ParameterSet + Copy + PartialEq> PartialEq for RingElement<P> {
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

    use crate::constants::parameter_sets::P768;

    use super::*;

    // REMARKS:
    // axiom tests:
    // -[x] closure under addition
    // -[x] commutative
    // -[x] identity
    // -[ ] associativity
    // -[ ] distributivity

    #[test]
    fn test_sample_poly_cbd() {
        let bytes = "example input bytes".as_bytes();

        let a: RingElement<P768> = RingElement::sample_poly_cbd(&bytes, 0x01);

        // testing against the great Filippo Valsorda https://github.com/FiloSottile/mlkem768
        let result = [
            3328, 3328, 3327, 3328, 3328, 2, 0, 3328, 3328, 0, 0, 0, 0, 3328, 2, 1, 0, 0, 3328, 1,
            0, 3328, 1, 3328, 3328, 0, 0, 0, 1, 3328, 1, 3328, 1, 0, 2, 0, 1, 0, 3328, 0, 2, 3328,
            3328, 2, 1, 3328, 1, 0, 0, 0, 3327, 0, 1, 2, 3328, 1, 0, 0, 3328, 0, 0, 0, 0, 1, 1, 1,
            1, 3328, 2, 3328, 3327, 3328, 0, 1, 0, 0, 0, 3328, 1, 3328, 0, 0, 0, 1, 1, 0, 3328, 0,
            1, 1, 0, 3328, 3328, 3328, 1, 1, 0, 0, 3328, 0, 3327, 0, 2, 1, 0, 3328, 1, 0, 3328,
            3328, 3328, 3328, 0, 1, 0, 0, 1, 2, 1, 0, 1, 0, 3328, 3328, 0, 2, 3328, 3328, 3328, 0,
            3328, 3328, 3328, 1, 3328, 1, 3328, 1, 0, 3328, 1, 3328, 2, 3327, 3328, 0, 3328, 2, 1,
            1, 0, 3328, 0, 3328, 0, 1, 3327, 0, 0, 1, 3328, 1, 3328, 1, 1, 3328, 2, 3328, 0, 0,
            3328, 1, 3327, 0, 3327, 1, 3328, 3327, 0, 0, 3328, 1, 0, 1, 1, 3328, 1, 3328, 0, 0, 0,
            2, 3328, 2, 3328, 0, 3328, 2, 3328, 0, 0, 1, 0, 3328, 2, 0, 0, 2, 3328, 3328, 0, 3327,
            3328, 1, 0, 0, 1, 0, 3328, 3328, 1, 3328, 3327, 0, 3327, 3328, 0, 0, 0, 0, 3328, 0, 1,
            0, 0, 0, 2, 3327, 1, 1, 0, 0, 1, 0, 0, 3327, 1, 3327, 0, 0, 0, 0, 1, 3327, 1, 1,
        ]
        .map(|val| F::new(val));

        let b = RingElement::new(result);
        assert_eq!(a.coefficients, b.coefficients);
    }

    #[test]
    fn test_additive_commutativity() {
        // Initialize a CRNG
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        // Step 1: Initialize a RingElement with random values
        let a: RingElement<P768> = RingElement::sample_poly_cbd(&bytes, 0xAA);
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

        let a: RingElement<P768> = RingElement::sample_poly_cbd(&bytes, 0xAA);
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
        let a: RingElement<P768> = RingElement::sample_poly_cbd(&bytes, 0xAA);
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
        let a_vals: [F<P768>; 256] = [F::new(123); 256];
        let a = RingElement::new(a_vals);
        let mut inverse_vals = [F::zero(); 256];
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
        let a: RingElement<P768> = RingElement::sample_poly_cbd(&bytes, 0xAA);

        // Encode the RingElement to bytes
        let encoded_bytes = a.byte_encode();

        // Decode the bytes back into a vector of F
        let decoded_elements: RingElement<P768> =
            RingElement::byte_decode(&encoded_bytes).expect("Decoding failed");

        // Verify the decoded vector matches the original RingElement's array
        assert_eq!(
            decoded_elements.coefficients.len(),
            a.coefficients.len(),
            "Length mismatch"
        );

        // Check each value for equality after encoding and decoding
        for (decoded_elem, original_elem) in decoded_elements
            .coefficients
            .iter()
            .zip(a.coefficients.iter())
        {
            assert_eq!(decoded_elem.val(), original_elem.val(), "Value mismatch");
        }
    }
}
