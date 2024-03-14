use core::marker::PhantomData;
use core::ops::{Add, AddAssign, Mul, Neg, Sub};

use crate::constants::ml_kem_constants::Q;
use crate::constants::{
    barrett_constants::{MULTIPLIER as bar_mul, SHIFT as bar_shift},
    parameter_sets::ParameterSet,
};

pub enum OperationError {
    UnreducedFieldElementError,
}

#[derive(Copy, Clone, Debug, PartialEq)]
/// An integer modulo Q
pub struct FieldElement<P> {
    val: u16,
    _marker: PhantomData<P>,
}

impl<P: ParameterSet> FieldElement<P> {
    pub fn new(val: u16) -> Self {
        let mut f = FieldElement {
            val,
            _marker: PhantomData,
        };
        f.reduce_once();
        f
    }

    pub fn zero() -> Self {
        FieldElement {
            val: 0,
            _marker: PhantomData,
        }
    }

    fn reduce_once(&mut self) {
        let mut x = self.val.wrapping_sub(Q);
        x = x.wrapping_add((x >> 15).wrapping_mul(Q));
        self.val = x;
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
    // TODO: sometimes these might need to be called with
    // values of du/dv that are different from param defs
    pub fn compress(&self) -> u16 {
        let dividend = (self.val as u64).wrapping_shl(P::Du.into());
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
        let mask = (1u64 << P::Du) - 1;
        (adjusted_quotient & mask) as u16
    }

    /// FIPS 203 (DRAFT), Definition 4.6
    pub fn decompress(&self) -> u16 {
        let dividend = self.val as u32;
        let dividend = dividend.wrapping_mul(Q.into());
        let mut quotient = dividend.wrapping_shr(P::Dv.into());
        quotient = quotient.wrapping_add((dividend.wrapping_shr(P::Dv as u32 - 1)) & 1);
        quotient as u16
    }

    fn barrett_reduce(product: u32) -> Self {
        let quotient: u32 = ((product as u64 * bar_mul as u64) >> bar_shift) as u32;
        Self::new((product - quotient * Q as u32) as u16)
    }
}

impl<P: ParameterSet + Copy> AddAssign for FieldElement<P> {
    fn add_assign(&mut self, other: Self) {
        self.val = (self.val + other.val) % Q;
    }
}

impl<P: ParameterSet> Neg for FieldElement<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let a = self.val;
        FieldElement::new(3329 - a)
    }
}

impl<P: ParameterSet> Mul<u16> for FieldElement<P> {
    type Output = Self;

    fn mul(self, other: u16) -> Self {
        let product = (self.val as u32) * (other as u32);
        Self::barrett_reduce(product)
    }
}

impl<P: ParameterSet> Mul<FieldElement<P>> for u16 {
    type Output = Self;

    fn mul(self, other: FieldElement<P>) -> Self {
        let product = (other.val as u32) * (self as u32);
        FieldElement::<P>::barrett_reduce(product).val
    }
}

impl<P: ParameterSet> Mul<FieldElement<P>> for FieldElement<P> {
    type Output = Self;

    fn mul(self, other: FieldElement<P>) -> Self {
        let product = (other.val as u32) * (self.val as u32);
        FieldElement::barrett_reduce(product)
    }
}

impl<P: ParameterSet> Add for FieldElement<P> {
    type Output = Self;

    // a + b % Q
    fn add(self, other: Self) -> Self {
        Self::new(self.val + other.val)
    }
}

impl<P: ParameterSet> Sub for FieldElement<P> {
    type Output = Self;

    // a - b % Q
    fn sub(self, other: Self) -> Self {
        // If `self.val` is less than `other.val`, adding `Q`
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
    use core::marker::PhantomData;

    use crate::{
        constants::{ml_kem_constants::Q, parameter_sets::P768},
        math::field_element::FieldElement as F,
    };

    // REMARK:
    // because the field Q is so small, it is actually possible
    // to test arithmetic operations on every single element
    // of the field

    #[test]
    fn exhaustive_test_reduce_once() {
        for i in Q + 1..=2 * Q {
            let element: F<P768> = F::new(i);
            assert!(
                element.val <= Q,
                "Value should be reduced within [0, Q] range for input: {i}, but got: {}",
                element.val
            );
        }
    }

    #[test]
    fn exhaustive_test_addition() {
        for i in 0..3329 {
            for j in 0..3329 {
                let a: F<P768> = F::new(i);
                let b: F<P768> = F::new(j);
                let result = a + b;
                assert_eq!(result.val, (i + j) % 3329);
            }
        }
    }

    #[test]
    fn exhaustive_test_subtraction() {
        for i in 0..3329 {
            for j in 0..3329 {
                let a: F<P768> = F::new(i);
                let b: F<P768> = F::new(j);
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
                let a: F<P768> = F::new(i);
                let b = j;
                let result = a * b;
                let expected = (i as u32 * j as u32) % 3329;
                assert_eq!(result.val, expected as u16, "Failed at i = {i} and j = {j}");
            }
        }
    }

    #[test]
    fn test_multiplication_with_potential_overflow() {
        let a: F<P768> = F::new(3000);
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
        assert!(F::<P768>::new(Q - 1).check_reduced().is_ok());
    }

    #[test]
    fn test_check_reduced_err() {
        assert!(F {
            val: Q + 1,
            _marker: PhantomData::<P768>
        }
        .check_reduced()
        .is_err());
    }

    // Test that verifies compression into a range with d = 10, where Q is assumed to be 3329.
    #[test]
    fn test_compress() {
        let test_cases = vec![
            // (input value, expected compressed output)
            (0, 0),      // A value that maps directly to 0
            (1664, 512), // A value around Q/2 should map near the middle of the range
            (3328, 0),   // A value around Q should map to the upper end of the range
        ];

        for (val, expected) in test_cases {
            let fe: F<P768> = F {
                val,
                _marker: PhantomData,
            };
            let compressed = fe.compress();
            assert_eq!(compressed, expected, "Compression of {} failed", val);
        }
    }

    // Test that verifies compression into a range with d = 10, where Q is assumed to be 3329.
    #[test]
    fn test_decompress() {
        let test_cases = vec![
            // (input value, expected compressed output)
            (0, 0),        // A value that maps directly to 0
            (1664, 18536), // A value around Q/2 should map near the middle of the range
            (3328, 37072), // A value around Q should map to the upper end of the range
        ];

        for (val, expected) in test_cases {
            let fe: F<P768> = F {
                val,
                _marker: PhantomData,
            };
            let compressed = fe.decompress();
            assert_eq!(compressed, expected, "Compression of {} failed", val);
        }
    }

    #[test]
    fn test_compress_with_mask() {
        let fe: F<P768> = F {
            val: 12345,
            _marker: PhantomData,
        };
        let compressed = fe.compress();
        assert!(
            compressed < 1024,
            "Compressed value should be within the mask limit"
        );
    }
}
