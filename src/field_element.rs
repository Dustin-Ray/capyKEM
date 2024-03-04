// use crypto_bigint::modular::ConstMontyParams;
// use crypto_bigint::{const_monty_form, impl_modulus, Uint, U64};
// use crypto_bigint::{rand_core::OsRng, NonZero, RandomMod};
// use rand::Rng;

use std::ops::{Add, Mul, Sub};

use crate::constants::{
    barrett_constants::MULTIPLIER as bar_mul, barrett_constants::SHIFT as bar_shift,
    ml_kem_constants::q,
};

pub enum OperationError {
    UnreducedFieldElementError,
}

/// An integer modulo q
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FieldElement {
    pub val: u16,
}

impl FieldElement {
    pub fn new(val: u16) -> Self {
        Self { val }
    }

    // Checks whether a given field element is less than q and thus reduced
    fn check_reduced(self) -> Result<Self, OperationError> {
        if self.val > q {
            Err(OperationError::UnreducedFieldElementError)
        } else {
            Ok(self)
        }
    }

    /// FIPS 203 (DRAFT), Definition 4.5.
    pub fn compress(&self, d: u8) -> u16 {
        // Use wrapping shift left for x * 2^d, ensuring no overflow occurs
        let dividend = (self.val as u64).wrapping_shl(d.into()); // Ensure `d` is converted to u32 for the shift

        // Use wrapping arithmetic for multiplication and shifting
        let quotient = dividend
            .wrapping_mul(bar_mul.into())
            .wrapping_shr(bar_shift as u32);

        // Calculate remainder using wrapping subtraction to ensure no underflow
        let remainder = dividend.wrapping_sub(quotient.wrapping_mul(q.into()));

        // Adjusting quotient based on the remainder's span, using wrapping arithmetic for addition
        let mut adjusted_quotient = quotient;
        if remainder > (q / 2).into() {
            adjusted_quotient = adjusted_quotient.wrapping_add(1);
        }
        if remainder > (q + q / 2).into() {
            adjusted_quotient = adjusted_quotient.wrapping_add(1);
        }

        // Apply mask to handle potential overflow using bitwise AND
        let mask = (1u64 << d) - 1;
        (adjusted_quotient & mask) as u16
    }

    /// FIPS 203 (DRAFT), Definition 4.6.
    pub fn decompress(&self, d: u8) -> u16 {
        // Use wrapping multiplication to avoid overflow panic
        let dividend = self.val as u32; // Cast to u32 for larger range
        let dividend = dividend.wrapping_mul(q.into());

        // Use wrapping shift right for the division part
        let mut quotient = dividend.wrapping_shr(d.into()); // Ensure `d` is cast to u32

        // Adjust quotient based on the rounding rule, using wrapping arithmetic
        quotient = quotient.wrapping_add((dividend.wrapping_shr((d - 1).into())) & 1);

        // The final result is cast back to u16
        quotient as u16
    }

    fn reduce_once(self) -> Self {
        let mut x = self.val.wrapping_sub(q);
        // If x underflowed, then x >= 2^16 - q > 2^15, so the top bit is set.
        // The wrapping_sub method ensures that underflow will wrap around in
        // a manner similar to unsigned arithmetic in C, without panicking.
        x = x.wrapping_add((x >> 15).wrapping_mul(q));
        Self::new(x)
    }

    // reduces self.val modulo q^2 avoiding potentially variable-time division
    fn barrett_reduce(self) -> Self {
        let quotient: u32 = ((self.val as u64 * bar_mul as u64) >> bar_shift) as u32;
        Self::new((self.val as u32 - quotient * q as u32) as u16).reduce_once()
    }
}

impl Add for FieldElement {
    type Output = Self;

    // a + b % q
    fn add(self, other: Self) -> Self {
        Self::new(self.val + other.val).reduce_once()
    }
}

impl Sub for FieldElement {
    type Output = Self;

    // a - b % q
    fn sub(self, other: Self) -> Self {
        // If `self.val` is less than `other.val`, adding `q`
        // ensures the result stays positive and wraps around correctly.
        let result = if self.val < other.val {
            self.val + q - other.val
        } else {
            self.val - other.val
        };
        Self::new(result).reduce_once()
    }
}

impl Mul<u16> for FieldElement {
    type Output = Self;
    fn mul(self, other: u16) -> Self {
        Self::new(self.val * other).barrett_reduce()
    }
}

#[cfg(test)]
mod tests {
    use super::q;
    use super::*;

    #[test]
    fn test_addition_within_bounds() {
        let a = FieldElement::new(1000);
        let b = FieldElement::new(1000);
        let result = a + b;
        assert_eq!(result.val, 2000 % q);
    }

    #[test]
    fn test_addition_overflow() {
        let a = FieldElement::new(q); // 3329
        let b = FieldElement::new(1);
        let result = a + b;
        // Since the addition will overflow beyond q, the reduction should bring it back to 0
        assert_eq!(result.val, (q + 1) % q); // effectively, q % q, which is 0
    }

    #[test]
    fn test_subtraction_within_bounds() {
        let a = FieldElement::new(2000);
        let b = FieldElement::new(1000);
        let result = a - b;
        // Expecting the direct subtraction result since 2000 - 1000 does not underflow
        assert_eq!(result.val, ((2000 + q) - 1000) % q);
    }

    #[test]
    fn test_subtraction_underflow() {
        let a = FieldElement::new(0);
        let b = FieldElement::new(1);
        let result = a - b;
        // Expecting q - 1 since 0 - 1 underflows and should wrap around modulo q
        assert_eq!(result.val, (q - 1) % q);
    }

    #[test]
    fn test_check_reduced_ok() {
        assert!(FieldElement::new(q - 1).check_reduced().is_ok());
    }

    #[test]
    fn test_check_reduced_err() {
        assert!(FieldElement::new(q + 1).check_reduced().is_err());
    }

    #[test]
    fn reduce_large_value() {
        let large_val = q as u32 * 2; // A value clearly larger than q to test reduction.
        let element = FieldElement::new(large_val as u16);
        let reduced = element.barrett_reduce();
        assert!(reduced.val < q); // Ensure the reduced value is less than q.
    }

    #[test]
    fn reduce_value_just_over_q() {
        let just_over_q = q as u32 + 1; // Just one more than q to see if it reduces to 1.
        let element = FieldElement::new(just_over_q as u16);
        let reduced = element.barrett_reduce();
        assert_eq!(reduced.val, just_over_q as u16 % q); // Should be equivalent to 1.
    }

    #[test]
    fn reduce_value_under_q() {
        // A value less than q should remain unchanged after reduction.
        let under_q = q as u32 - 1;
        let element = FieldElement::new(under_q as u16);
        let reduced = element.barrett_reduce();
        assert_eq!(reduced.val, under_q as u16); // Expect the original value, as it's already less than q.
    }

    #[test]
    fn reduce_exact_multiple_of_q() {
        // An exact multiple of q should reduce to 0.
        let multiple_q = q as u32 * 3; // Three times q for testing.
        let element = FieldElement::new(multiple_q as u16);
        let reduced = element.barrett_reduce();
        assert_eq!(reduced.val, (multiple_q % q as u32) as u16); // Should be 0 if perfectly reduced.
    }

    // Test that verifies compression into a range with `d=10`, where q is assumed to be 3329.
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

    // Test that verifies compression into a range with `d=10`, where q is assumed to be 3329.
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
        // Test to ensure masking works correctly
        let fe = FieldElement { val: 12345 }; // A larger value to test masking
        let compressed = fe.compress(10); // `d` value that results in a specific mask
                                          // The expected value depends on the mask applied, adjust this based on your logic
        assert!(
            compressed < 1024,
            "Compressed value should be within the mask limit"
        );
    }

    #[test]
    fn test_decompress_then_compress() {
        let d = 10; // Bit width for the operations
        let original_values = vec![0, 1, 1023]; // Example original values within the range of 0 to 2^d - 1

        for &original in &original_values {
            // Simulate a FieldElement from the original value
            let fe = FieldElement::new(original);

            // Decompress the value
            let decompressed = fe.decompress(d);

            // Create a FieldElement from the decompressed value for compression
            let decompressed_fe = FieldElement::new(decompressed);

            // Compress the decompressed value
            let compressed = decompressed_fe.compress(d);

            // Check if the compressed value matches the original value
            assert_eq!(
                compressed, original,
                "Original value: {}, Decompressed then Compressed value: {}",
                original, compressed
            );
        }
    }
}
