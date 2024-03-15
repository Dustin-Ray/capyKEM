#[cfg(test)]
mod tests {

    use alloc::vec::Vec;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    use crate::{
        constants::{parameter_sets::P768, sample_poly_cbd_result},
        math::{field_element::FieldElement as F, ring_element::RingElement},
    };

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
        let b = RingElement::new(sample_poly_cbd_result.map(|val| F::new(val)));
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
}
