// use crypto_bigint::modular::ConstMontyParams;
// use crypto_bigint::{const_monty_form, impl_modulus, Uint, U64};
// use crypto_bigint::{rand_core::OsRng, NonZero, RandomMod};
// use rand::Rng;

use core::ops::{Add, AddAssign, Mul, Neg, Sub};

use crate::constants::{
    barrett_constants::{MULTIPLIER as bar_mul, SHIFT as bar_shift},
    ml_kem_constants::Q,
};

pub enum OperationError {
    UnreducedFieldElementError,
}

/// An integer modulo q
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FieldElement {
    val: u16,
}

impl FieldElement {
    pub fn new(val: u16) -> Self {
        let mut f = FieldElement { val };
        f.reduce_once();
        f
    }

    pub fn zero() -> Self {
        FieldElement::new(0)
    }

    // Change to modify `self` in place instead of consuming it
    fn reduce_once(&mut self) {
        let mut x = self.val.wrapping_sub(Q);
        x = x.wrapping_add((x >> 15).wrapping_mul(Q));
        self.val = x; // Apply the modification directly to `self`
    }

    pub fn default() -> Self {
        FieldElement { val: 0 }
    }

    pub fn check_reduced(self) -> Result<Self, OperationError> {
        if self.val > Q {
            Err(OperationError::UnreducedFieldElementError)
        } else {
            Ok(self)
        }
    }

    pub fn val(self) -> u16 {
        self.val
    }

    // FIPS 203 (DRAFT), Definition 4.5.
    // REMARKS:
    // parameterize du and dv
    pub fn compress(&self, d: u8) -> u16 {
        let dividend = (self.val as u64).wrapping_shl(d.into());
        let quotient = dividend
            .wrapping_mul(bar_mul.into())
            .wrapping_shr(bar_shift as u32);
        let remainder = dividend.wrapping_sub(quotient.wrapping_mul(Q.into()));
        let mut adjusted_quotient = quotient;
        if remainder > (Q / 2).into() {
            adjusted_quotient = adjusted_quotient.wrapping_add(1);
        }
        if remainder > (Q + Q / 2).into() {
            adjusted_quotient = adjusted_quotient.wrapping_add(1);
        }
        let mask = (1u64 << d) - 1;
        (adjusted_quotient & mask) as u16
    }

    /// FIPS 203 (DRAFT), Definition 4.6
    pub fn decompress(&self, d: u8) -> u16 {
        let dividend = self.val as u32;
        let dividend = dividend.wrapping_mul(Q.into());
        let mut quotient = dividend.wrapping_shr(d.into());
        quotient = quotient.wrapping_add((dividend.wrapping_shr((d - 1).into())) & 1);
        quotient as u16
    }

    fn barrett_reduce(product: u32) -> Self {
        let quotient: u32 = ((product as u64 * bar_mul as u64) >> bar_shift) as u32;
        Self::new((product - quotient * Q as u32) as u16)
    }
}

impl AddAssign for FieldElement {
    fn add_assign(&mut self, other: Self) {
        self.val = (self.val + other.val) % Q;
    }
}

impl Neg for FieldElement {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let a = self.val;
        FieldElement::new(3329 - a)
    }
}

impl Mul<u16> for FieldElement {
    type Output = Self;

    fn mul(self, other: u16) -> Self {
        // Perform multiplication in a larger integer type to handle overflow
        let product = (self.val as u32) * (other as u32);
        Self::barrett_reduce(product)
    }
}

impl Mul<FieldElement> for u16 {
    type Output = Self;

    fn mul(self, other: FieldElement) -> Self {
        // Perform multiplication in a larger integer type to handle overflow
        let product = (other.val as u32) * (self as u32);
        FieldElement::barrett_reduce(product).val
    }
}

impl Mul<FieldElement> for FieldElement {
    type Output = Self;

    fn mul(self, other: FieldElement) -> Self {
        // Perform multiplication in a larger integer type to handle overflow
        let product = (other.val as u32) * (self.val as u32);
        FieldElement::barrett_reduce(product)
    }
}

impl Add for FieldElement {
    type Output = Self;

    // a + b % q
    fn add(self, other: Self) -> Self {
        Self::new(self.val + other.val)
    }
}

impl Sub for FieldElement {
    type Output = Self;

    // a - b % q
    fn sub(self, other: Self) -> Self {
        // If `self.val` is less than `other.val`, adding `q`
        // ensures the result stays positive and wraps around correctly.
        let result = if self.val < other.val {
            self.val + Q - other.val
        } else {
            self.val - other.val
        };
        Self::new(result)
    }
}

#[cfg(test)]
mod tests {
    use super::Q;
    use super::*;

    #[test]
    fn exhaustive_test_reduce_once() {
        for i in Q + 1..=2 * Q {
            let element = FieldElement::new(i);
            assert!(
                element.val <= Q,
                "Value should be reduced within [0, q] range for input: {i}, but got: {}",
                element.val
            );
        }
    }

    #[test]
    fn exhaustive_test_addition() {
        for i in 0..3329 {
            for j in 0..3329 {
                let a = FieldElement::new(i);
                let b = FieldElement::new(j);
                let result = a + b;
                assert_eq!(result.val, (i + j) % 3329);
            }
        }
    }

    #[test]
    fn exhaustive_test_subtraction() {
        for i in 0..3329 {
            for j in 0..3329 {
                let a = FieldElement::new(i);
                let b = FieldElement::new(j);
                let result = a - b;
                let expected = (i as isize - j as isize + 3329 as isize) % 3329;

                assert_eq!(result.val, expected as u16, "Failed at i = {i} and j = {j}");
            }
        }
    }

    #[test]
    fn exhaustive_test_multiplication() {
        for i in 0..3329 {
            for j in 0..3329 {
                let a = FieldElement::new(i);
                let b = j;
                let result = a * b;
                let expected = (i as u32 * j as u32) % 3329;
                assert_eq!(result.val, expected as u16, "Failed at i = {i} and j = {j}");
            }
        }
    }

    #[test]
    fn test_multiplication_with_potential_overflow() {
        let a = FieldElement::new(3000);
        let b = 3000;
        let result = a * b;
        let expected = (3000u32 * 3000u32) % 3329u32;
        assert_eq!(
            result.val, expected as u16,
            "Multiplication resulted in overflow or incorrect handling"
        );
    }

    #[test]
    fn test_check_reduced_ok() {
        assert!(FieldElement::new(Q - 1).check_reduced().is_ok());
    }

    #[test]
    fn test_check_reduced_err() {
        assert!(FieldElement { val: Q + 1 }.check_reduced().is_err());
    }

    // Test that verifies compression into a range with d = 10, where q is assumed to be 3329.
    #[test]
    fn test_compress() {
        let d = 10; // Example bit width
        let test_cases = vec![
            // (input value, expected compressed output)
            (0, 0),      // A value that maps directly to 0
            (1664, 512), // A value around q/2 should map near the middle of the range
            (3328, 0),   // A value around q should map to the upper end of the range
        ];

        for (val, expected) in test_cases {
            let fe = FieldElement { val };
            let compressed = fe.compress(d);
            assert_eq!(compressed, expected, "Compression of {} failed", val);
        }
    }

    // Test that verifies compression into a range with d = 10, where q is assumed to be 3329.
    #[test]
    fn test_decompress() {
        let d = 10; // Example bit width
        let test_cases = vec![
            // (input value, expected compressed output)
            (0, 0),        // A value that maps directly to 0
            (1664, 5410),  // A value around q/2 should map near the middle of the range
            (3328, 10819), // A value around q should map to the upper end of the range
        ];

        for (val, expected) in test_cases {
            let fe = FieldElement { val };
            let compressed = fe.decompress(d);
            assert_eq!(compressed, expected, "Compression of {} failed", val);
        }
    }

    #[test]
    fn test_compress_with_mask() {
        let fe = FieldElement { val: 12345 };
        let compressed = fe.compress(10);
        assert!(
            compressed < 1024,
            "Compressed value should be within the mask limit"
        );
    }

    #[test]
    fn test_decompress_then_compress() {
        let d = 10;
        // Example original values within the range of 0 to 2^d - 1
        let original_values = vec![0, 1, 1023];

        for &original in &original_values {
            let fe = FieldElement::new(original);
            let decompressed = fe.decompress(d);
            let decompressed_fe = FieldElement::new(decompressed);
            let compressed = decompressed_fe.compress(d);
            assert_eq!(
                compressed, original,
                "Original value: {}, Decompressed then Compressed value: {}",
                original, compressed
            );
        }
    }
}
