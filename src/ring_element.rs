use std::ops::{Add, Sub};

use crate::field_element::FieldElement;

/// A polynomial is an element of the ring R. It is an array of 256 coefficients
/// which themselves are [FieldElement].
#[derive(Clone, Debug)]
pub struct RingElement {
    pub val: Vec<FieldElement>,
}

impl RingElement {
    // Create a new RingElement from a vector of FieldElements
    pub fn new(val: Vec<FieldElement>) -> Self {
        RingElement { val }
    }

    fn byte_encode(&mut self, d: u8) -> Vec<u8> {
        assert!(d <= 12, "d must be less than or equal to 12");
        assert_eq!(self.val.len(), 256, "Input array must have 256 elements");

        let mut b = vec![0u8; 0];

        todo!()
    }
}

impl Add for RingElement {
    type Output = Self;

    // Simple polynomial addition for two ring elements of the same length
    fn add(self, other: Self) -> Self {
        assert_eq!(
            self.val.len(),
            other.val.len(),
            "RingElements must be of the same length"
        );
        RingElement::new(
            self.val
                .iter()
                .zip(other.val.iter())
                .map(|(x, y)| *x + *y)
                .collect(),
        )
    }
}

impl Sub for RingElement {
    type Output = Self;

    // Simple polynomial subtraction for two ring elements of the same length
    fn sub(self, other: Self) -> Self {
        assert_eq!(
            self.val.len(),
            other.val.len(),
            "RingElements must be of the same length"
        );

        RingElement::new(
            self.val
                .iter()
                .zip(other.val.iter())
                .map(|(x, y)| *x - *y)
                .collect(),
        )
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
    use super::*;

    #[test]
    fn test_addition() {
        let a = RingElement::new(vec![FieldElement::new(1), FieldElement::new(2)]);
        let b = RingElement::new(vec![FieldElement::new(3), FieldElement::new(4)]);
        let expected = RingElement::new(vec![FieldElement::new(4), FieldElement::new(6)]);
        assert_eq!(a + b, expected);
    }

    #[test]
    fn test_subtraction() {
        let a = RingElement::new(vec![FieldElement::new(5), FieldElement::new(7)]);
        let b = RingElement::new(vec![FieldElement::new(2), FieldElement::new(3)]);
        let expected = RingElement::new(vec![FieldElement::new(3), FieldElement::new(4)]);
        assert_eq!(a - b, expected);
    }

    #[test]
    #[should_panic(expected = "RingElements must be of the same length")]
    fn test_different_length() {
        let a = RingElement::new(vec![FieldElement::new(1)]);
        let b = RingElement::new(vec![FieldElement::new(2), FieldElement::new(3)]);
        let _ = a + b; // This should panic
    }
}
