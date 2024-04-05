use core::ops::{Add, AddAssign, Mul, Neg, Sub};

use crate::constants::barrett_constants::{MULTIPLIER as bar_mul, SHIFT as bar_shift};
use crate::constants::ml_kem_constants::q;

use super::encoding::{Compress, CompressionFactor};

pub enum OperationError {
    UnreducedFieldElementError,
}

#[derive(Copy, Clone, Debug, PartialEq)]
/// An integer modulo Q
pub struct FieldElement(pub u16);

impl FieldElement {
    pub fn new(val: u16) -> Self {
        let mut f = FieldElement(val);
        f.reduce_once();
        f
    }

    /// This should ONLY be used when certain that
    /// val is < q
    pub fn from(val: u16) -> Self {
        FieldElement(val)
    }

    pub fn zero() -> Self {
        FieldElement(0)
    }

    pub fn set(&mut self, val: u16) {
        self.0 = val;
    }

    pub fn reduce_once(&mut self) {
        let mut x = self.val().wrapping_sub(q);
        x = x.wrapping_add((x >> 15).wrapping_mul(q));
        self.0 = x;
    }

    pub fn check_reduced(self) -> Result<Self, OperationError> {
        if self.val() > q {
            Err(OperationError::UnreducedFieldElementError)
        } else {
            Ok(self)
        }
    }

    pub fn val(self) -> u16 {
        self.0
    }

    // FIPS 203 (DRAFT), Definition 4.5.
    // TODO: sometimes these might need to be called with
    // values of du/dv that are different from param defs
    pub fn compress<const d: u16>(&self) -> u16 {
        let dividend = u64::from(self.val()).wrapping_shl(d.into());
        let quotient = dividend
            .wrapping_mul(bar_mul.into())
            .wrapping_shr(u32::from(bar_shift));
        let remainder = dividend.wrapping_sub(quotient.wrapping_mul(q.into()));
        let mut adjusted_quotient = quotient;
        if remainder > (q / 2).into() {
            adjusted_quotient = adjusted_quotient.wrapping_add(1);
        }
        if remainder > (q + q / 2).into() {
            adjusted_quotient = adjusted_quotient.wrapping_add(1);
        }
        let mask = (1u64 << d) - 1;
        (adjusted_quotient & mask) as u16
    }

    /// FIPS 203 (DRAFT), Definition 4.6
    pub fn decompress<const d: u16>(y: u16) -> Self {
        let dividend = u32::from(y);
        let dividend = dividend.wrapping_mul(q.into());
        let mut quotient = dividend.wrapping_shr(d.into());
        quotient = quotient.wrapping_add((dividend.wrapping_shr(u32::from(d) - 1)) & 1);
        Self::from(quotient as u16)
    }

    fn barrett_reduce(product: u32) -> Self {
        let quotient: u32 = ((u64::from(product) * u64::from(bar_mul)) >> bar_shift) as u32;
        Self::new((product - quotient * u32::from(q)) as u16)
    }
}

impl Compress for FieldElement {
    // FIPS 203 (DRAFT), Definition 4.5.
    // TODO: sometimes these might need to be called with
    // values of du/dv that are different from param defs
    fn compress<D: CompressionFactor>(&mut self) -> &FieldElement {
        let dividend = u64::from(self.val()).wrapping_shl(D::USIZE as u32);
        let quotient = dividend
            .wrapping_mul(bar_mul.into())
            .wrapping_shr(u32::from(bar_shift));
        let remainder = dividend.wrapping_sub(quotient.wrapping_mul(q.into()));
        let mut adjusted_quotient = quotient;
        if remainder > (q / 2).into() {
            adjusted_quotient = adjusted_quotient.wrapping_add(1);
        }
        if remainder > (q + q / 2).into() {
            adjusted_quotient = adjusted_quotient.wrapping_add(1);
        }
        let mask = (1u64 << D::USIZE) - 1;
        self.0 = (adjusted_quotient & mask) as u16;

        self
    }

    /// FIPS 203 (DRAFT), Definition 4.6
    fn decompress<D: CompressionFactor>(&mut self) -> &FieldElement {
        let dividend = u32::from(self.val());
        let dividend = dividend.wrapping_mul(q.into());
        let mut quotient = dividend.wrapping_shr(D::USIZE as u32);
        quotient = quotient.wrapping_add((dividend.wrapping_shr(D::USIZE as u32 - 1)) & 1);
        self.0 = quotient as u16;
        self
    }
}

impl AddAssign for FieldElement {
    fn add_assign(&mut self, other: Self) {
        self.0 = (self.val() + other.val()) % q;
    }
}

impl Neg for FieldElement {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let a = self.val();
        FieldElement::new(3329 - a)
    }
}

impl Mul<u16> for FieldElement {
    type Output = Self;

    fn mul(self, other: u16) -> Self {
        let product = u32::from(self.val()) * u32::from(other);
        Self::barrett_reduce(product)
    }
}

impl Mul<FieldElement> for u16 {
    type Output = Self;

    fn mul(self, other: FieldElement) -> Self {
        let product = u32::from(other.val()) * u32::from(self);
        FieldElement::barrett_reduce(product).val()
    }
}

impl Mul<FieldElement> for FieldElement {
    type Output = Self;

    fn mul(self, other: FieldElement) -> Self {
        let product = u32::from(other.val()) * u32::from(self.val());
        FieldElement::barrett_reduce(product)
    }
}

impl Add for FieldElement {
    type Output = Self;

    // a + b % Q
    fn add(self, other: Self) -> Self {
        Self::new(self.val() + other.val())
    }
}

impl Sub for FieldElement {
    type Output = Self;

    // a - b % Q
    fn sub(self, other: Self) -> Self {
        // If `self.val()` is less than `other.val()`, adding `Q`
        // ensures the result stays positive and wraps around correctly.
        let result = if self.val() < other.val() {
            self.val() + q - other.val()
        } else {
            self.val() - other.val()
        };
        Self::new(result)
    }
}

#[cfg(test)]
mod tests {

    use crate::{constants::ml_kem_constants::q, math::field_element::FieldElement as F};

    // REMARK:
    // because the field Q is so small, it is actually possible
    // to test arithmetic operations on every single element
    // of the field

    #[test]
    fn exhaustive_test_reduce_once() {
        for i in q + 1..=2 * q {
            let element = F::new(i);
            assert!(
                element.val() <= q,
                "Value should be reduced within [0, Q] range for input: {i}, but got: {}",
                element.val()
            );
        }
    }

    #[test]
    fn exhaustive_test_addition() {
        for i in 0..3329 {
            for j in 0..3329 {
                let a = F::new(i);
                let b = F::new(j);
                let result = a + b;
                assert_eq!(result.val(), (i + j) % 3329);
            }
        }
    }

    #[test]
    fn exhaustive_test_subtraction() {
        for i in 0..3329 {
            for j in 0..3329 {
                let a = F::new(i);
                let b = F::new(j);
                let result = a - b;
                let expected = (i as isize - j as isize + 3329_isize) % 3329;

                assert_eq!(
                    result.val(),
                    expected as u16,
                    "Failed at i = {i} and j = {j}"
                );
            }
        }
    }

    #[test]
    fn exhaustive_test_multiplication() {
        for i in 0..3329 {
            for j in 0..3329 {
                let a = F::new(i);
                let b = j;
                let result = a * b;
                let expected = (i as u32 * j as u32) % 3329;
                assert_eq!(
                    result.val(),
                    expected as u16,
                    "Failed at i = {i} and j = {j}"
                );
            }
        }
    }

    #[test]
    fn test_multiplication_with_potential_overflow() {
        let a = F::new(3000);
        let b = 3000;
        let result = a * b;
        let expected = (3000u32 * 3000u32) % 3329u32;
        assert_eq!(
            result.val(),
            expected as u16,
            "Multiplication resulted in overflow or incorrect handling"
        );
    }

    #[test]
    fn test_check_reduced_ok() {
        assert!(F::new(q - 1).check_reduced().is_ok());
    }

    #[test]
    fn test_check_reduced_err() {
        assert!(F(q + 1).check_reduced().is_err());
    }

    // Test that verifies compression into a range with d = 10, where Q is assumed to be 3329.
    #[test]
    fn test_compress() {
        use alloc::vec;
        let test_cases = vec![
            // (input value, expected compressed output)
            (0, 0),      // A value that maps directly to 0
            (1664, 512), // A value around Q/2 should map near the middle of the range
            (3328, 0),   // A value around Q should map to the upper end of the range
        ];

        for (val, expected) in test_cases {
            let fe = F(val);
            let compressed = fe.compress::<10>();
            assert_eq!(compressed, expected, "Compression of {} failed", val);
        }
    }

    // Test that verifies compression into a range with d = 10, where Q is assumed to be 3329.
    #[test]
    fn test_decompress() {
        use alloc::vec;
        let test_cases = vec![
            // (input value, expected compressed output)
            (0, 0),        // A value that maps directly to 0
            (1664, 18536), // A value around Q/2 should map near the middle of the range
            (3328, 37072), // A value around Q should map to the upper end of the range
        ];

        for (val, expected) in test_cases {
            let compressed = F::decompress::<4>(val);
            assert_eq!(compressed.val(), expected, "Compression of {} failed", val);
        }
    }

    #[test]
    fn test_compress_with_mask() {
        let fe = F(12345);
        let compressed = fe.compress::<10>();
        assert!(
            compressed < 1024,
            "Compressed value should be within the mask limit"
        );
    }
}
