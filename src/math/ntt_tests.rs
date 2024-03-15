#[cfg(test)]
mod tests {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha20Rng;

    extern crate alloc;
    use crate::{
        constants::{parameter_sets::P768, sample_ntt_result},
        math::{field_element::FieldElement as F, ntt_element::NttElement, ring_element::RingElement},
    };
    use alloc::vec::Vec;

    // REMARKS:
    // axiom tests:
    // -[ ] commutative
    // -[x] zero identity
    // -[x] associativity
    // -[ ] distributivity
    // -[x] multiplication closure

    #[test]
    fn test_sample_ntt() {
        let byte_stream = [0_u8; 0];
        let a: NttElement<P768> = NttElement::sample_ntt(&byte_stream, 0, 1);
        assert_eq!(a.ring, sample_ntt_result.map(|val| F::new(val)));
    }

    #[test]
    fn test_ntt() {
        use alloc::vec;
        // sample output is in NTT domain
        let mut byte_stream: NttElement<P768> = NttElement::sample_ntt(&vec![42_u8; 32], 1, 1);
        let mut byte_stream_copy = byte_stream;
        byte_stream.ntt_inv();
        byte_stream_copy.ntt_inv();
        assert_eq!(byte_stream_copy.ring, byte_stream.ring)
    }

    #[test]
    fn test_ntt_from_poly_cbd_inverse_with_random_input() {
        // Generate a random byte stream using a seeded RNG for reproducibility
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();
        // Sample a ring element using the random byte stream
        let mut ring_element: RingElement<P768> = RingElement::sample_poly_cbd(&bytes, 0xFF);
        let ring_element_copy = ring_element;

        // runs .ntt() on intstantiation
        let ntt_element = NttElement::new(&mut ring_element);
        assert_eq!(
            ring_element_copy.coefficients,
            Into::<RingElement<P768>>::into(ntt_element).coefficients
        );
    }

    #[test]
    fn test_multiply_ntts_associative() {
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        let a: NttElement<P768> = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
        let b = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xBB));
        let c = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xCC));

        // Test associativity (ab)c = a(bc)
        let ab_c = (a * b) * c;
        let a_bc = a * (b * c);
        assert_eq!(ab_c.ring, a_bc.ring);

        let a: NttElement<P768> = NttElement::sample_ntt(&bytes.clone(), 0, 0);
        let b = NttElement::sample_ntt(&bytes.clone(), 0, 0);
        let c = NttElement::sample_ntt(&bytes.clone(), 0, 0);

        // Test associativity (ab)c = a(bc)
        let ab_c = (a * b) * c;
        let a_bc = a * (b * c);
        assert_eq!(ab_c.ring, a_bc.ring);
    }

    #[test]
    fn test_multiply_ntts_zero() {
        let bytes: Vec<u8> = (0..32)
            .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
            .collect();

        let a: NttElement<P768> = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
        let zero = NttElement::zero();

        // Test multiplicative identity
        let res = zero * a;
        assert_eq!(
            res,
            NttElement {
                ring: [F::zero(); 256]
            }
        );
        let res = a * zero;
        assert_eq!(
            res,
            NttElement {
                ring: [F::zero(); 256]
            }
        );
    }

    #[test]
    fn test_closure_under_multiplication() {
        for _ in 0..1000 {
            let bytes: Vec<u8> = (0..32)
                .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
                .collect();

            let a: NttElement<P768> =
                NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
            let b = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xBB));

            let result = a * b;
            assert!(
                result.ring.iter().all(|x| x.val() < 3329),
                "Result of multiplication must be valid NttElement"
            );
        }
    }
}
