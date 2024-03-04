use std::ops::{Add, Sub};

use crate::field_element::FieldElement;

/// A polynomial is an element of the ring R. It is an array of 256 coefficients
/// which themselves are [FieldElement].
#[derive(Clone, Debug)]
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
