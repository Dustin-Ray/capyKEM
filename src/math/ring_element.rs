use std::ops::{Add, Sub};

use crate::{constants::ml_kem_constants, math::field_element::FieldElement};

/// A polynomial is an element of the ring R. It is an array of 256 coefficients
/// which themselves are [FieldElement].
#[derive(Clone, Copy, Debug)]
pub struct RingElement {
    pub val: [FieldElement; 256],
}

impl RingElement {
    // Create a new RingElement from a vector of FieldElements
    pub fn new(val: [FieldElement; 256]) -> Self {
        RingElement { val }
    }

    fn byte_encode(self) -> Vec<u8> {
        let mut out = Vec::with_capacity(256 * 12 / 8); // Preallocate the output vector

        for i in (0..self.val.len()).step_by(2) {
            // Combine two 12-bit integers into a single 24-bit integer
            let x = u32::from(self.val[i].val) | (u32::from(self.val[i + 1].val) << 12);

            // Split the 24-bit integer into 3 bytes and append to the output vector
            out.push((x & 0xFF) as u8); // First 8 bits
            out.push(((x >> 8) & 0xFF) as u8); // Next 8 bits
            out.push(((x >> 16) & 0xFF) as u8); // Last 8 bits
        }
        out
    }

    fn byte_decode(b: &[u8]) -> Result<Vec<FieldElement>, &'static str> {
        if b.len() != ml_kem_constants::ENCODE_SIZE_12.into() {
            return Err("Invalid encoding length");
        }

        let mut f = Vec::with_capacity(ml_kem_constants::N.into());
        let mut i = 0;
        while i < b.len() {
            let d = u32::from(b[i]) | (u32::from(b[i + 1]) << 8) | (u32::from(b[i + 2]) << 16);
            const MASK_12: u32 = 0b1111_1111_1111;

            let elem1 = FieldElement::new((d & MASK_12) as u16)
                .check_reduced()
                .map_err(|_| "Invalid polynomial encoding")?;

            let elem2 = FieldElement::new((d >> 12) as u16)
                .check_reduced()
                .map_err(|_| "Invalid polynomial encoding")?;

            f.push(elem1);
            f.push(elem2);

            i += 3;
        }

        Ok(f)
    }
}

impl Add for RingElement {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        assert_eq!(
            self.val.len(),
            other.val.len(),
            "RingElements must be of the same length"
        );

        let mut result = [FieldElement::default(); 256]; // Assuming FieldElement has a default value
        for i in 0..256 {
            result[i] = self.val[i] + other.val[i];
        }
        RingElement::new(result)
    }
}

impl Sub for RingElement {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        assert_eq!(
            self.val.len(),
            other.val.len(),
            "RingElements must be of the same length"
        );

        let mut result = [FieldElement::default(); 256]; // Assuming FieldElement has a default value
        for i in 0..256 {
            result[i] = self.val[i] - other.val[i];
        }
        RingElement::new(result)
    }
}

impl PartialEq for RingElement {
    fn eq(&self, other: &Self) -> bool {
        if self.val.len() != other.val.len() {
            return false;
        }
        self.val.iter().zip(other.val.iter()).all(|(a, b)| a == b)
    }
}

#[cfg(test)]
mod tests {
    use rand::{thread_rng, Rng};

    use super::*;

    #[test]
    fn test_encode_decode() {
        // Initialize a RNG
        let mut rng = thread_rng();

        // Step 1: Initialize a RingElement with random values
        let mut original_values = [FieldElement { val: 0 }; 256];
        for elem in original_values.iter_mut() {
            // Generate a random value within the 12-bit range and assign it
            elem.val = rng.gen_range(0..4096) % 3029;
        }
        let original_ring_element = RingElement {
            val: original_values,
        };

        // Step 2: Encode the RingElement to bytes
        let encoded_bytes = original_ring_element.byte_encode();

        // Step 3: Decode the bytes back into a vector of FieldElement
        let decoded_elements = RingElement::byte_decode(&encoded_bytes).expect("Decoding failed");

        // Step 4: Verify the decoded vector matches the original RingElement's array
        assert_eq!(
            decoded_elements.len(),
            original_ring_element.val.len(),
            "Length mismatch"
        );

        for (decoded_elem, original_elem) in decoded_elements
            .iter()
            .zip(original_ring_element.val.iter())
        {
            assert_eq!(decoded_elem.val, original_elem.val, "Value mismatch");
        }
    }
}
