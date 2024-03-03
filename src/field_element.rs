// use crypto_bigint::modular::ConstMontyParams;
// use crypto_bigint::{const_monty_form, impl_modulus, Uint, U64};
// use crypto_bigint::{rand_core::OsRng, NonZero, RandomMod};
// use rand::Rng;

use std::ops::{Add, Mul, Sub};

use crate::constants::{barrett_constants, ml_kem_constants};

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
        if self.val > ml_kem_constants::q {
            Err(OperationError::UnreducedFieldElementError)
        } else {
            Ok(self)
        }
    }

    // See FIPS 203 (DRAFT), Definition 4.5.
    // "Compresses" a field element into a byte representation defined as the following:

    // fn compress(self, d: u16) -> Self {

    // }

    fn reduce_once(self) -> Self {
        let q = ml_kem_constants::q;
        let mut x = self.val.wrapping_sub(q);
        // If x underflowed, then x >= 2^16 - q > 2^15, so the top bit is set.
        // The wrapping_sub method ensures that underflow will wrap around in
        // a manner similar to unsigned arithmetic in C, without panicking.
        x = x.wrapping_add((x >> 15).wrapping_mul(q));
        Self::new(x)
    }

    // reduces self.val modulo q^2 avoiding potentially variable-time division
    fn barrett_reduce(self) -> Self {
        let quotient: u32 = ((self.val as u64 * barrett_constants::MULTIPLIER as u64)
            >> barrett_constants::SHIFT) as u32;
        Self::new((self.val as u32 - quotient * ml_kem_constants::q as u32) as u16).reduce_once()
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
            self.val + ml_kem_constants::q - other.val
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
    use super::ml_kem_constants::q;
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
        let large_val = ml_kem_constants::q as u32 * 2; // A value clearly larger than q to test reduction.
        let element = FieldElement::new(large_val as u16);
        let reduced = element.barrett_reduce();
        assert!(reduced.val < ml_kem_constants::q); // Ensure the reduced value is less than q.
    }

    #[test]
    fn reduce_value_just_over_q() {
        let just_over_q = ml_kem_constants::q as u32 + 1; // Just one more than q to see if it reduces to 1.
        let element = FieldElement::new(just_over_q as u16);
        let reduced = element.barrett_reduce();
        assert_eq!(reduced.val, just_over_q as u16 % ml_kem_constants::q); // Should be equivalent to 1.
    }

    #[test]
    fn reduce_value_under_q() {
        // A value less than q should remain unchanged after reduction.
        let under_q = ml_kem_constants::q as u32 - 1;
        let element = FieldElement::new(under_q as u16);
        let reduced = element.barrett_reduce();
        assert_eq!(reduced.val, under_q as u16); // Expect the original value, as it's already less than q.
    }

    #[test]
    fn reduce_exact_multiple_of_q() {
        // An exact multiple of q should reduce to 0.
        let multiple_q = ml_kem_constants::q as u32 * 3; // Three times q for testing.
        let element = FieldElement::new(multiple_q as u16);
        let reduced = element.barrett_reduce();
        assert_eq!(
            reduced.val,
            (multiple_q % ml_kem_constants::q as u32) as u16
        ); // Should be 0 if perfectly reduced.
    }
}
