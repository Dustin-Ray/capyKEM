// #[cfg(test)]
// mod tests {
//     use rand::{Rng, SeedableRng};
//     use rand_chacha::ChaCha20Rng;

//     extern crate alloc;
//     use crate::{
//         constants::{parameter_sets::P768, sample_ntt_result},
//         math::{
//             field_element::FieldElement as F, ntt_element::NttElement, ring_element::RingElement,
//         },
//     };
//     use alloc::vec::Vec;

//     // REMARKS:
//     // axiom tests:
//     // -[ ] commutative
//     // -[x] zero identity
//     // -[x] associativity
//     // -[ ] distributivity
//     // -[x] multiplication closure

//     #[test]
//     fn test_sample_ntt() {
//         let byte_stream = [0_u8; 0];
//         let a: NttElement<P768> = NttElement::sample_ntt(&byte_stream, 0, 1);
//         assert_eq!(a.coefficients, sample_ntt_result.map(|val| F::new(val)));
//     }

//     #[test]
//     fn test_ntt() {
//         use alloc::vec;
//         // sample output is in NTT domain
//         let mut byte_stream: NttElement<P768> = NttElement::sample_ntt(&vec![42_u8; 32], 1, 1);
//         let mut byte_stream_copy = byte_stream;
//         byte_stream.ntt_inv();
//         byte_stream_copy.ntt_inv();
//         assert_eq!(byte_stream_copy.coefficients, byte_stream.coefficients)
//     }

//     #[test]
//     fn test_ntt_from_poly_cbd_inverse_with_random_input() {
//         // Generate a random byte stream using a seeded RNG for reproducibility
//         let bytes: Vec<u8> = (0..32)
//             .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
//             .collect();
//         // Sample a ring element using the random byte stream
//         let mut ring_element: RingElement<P768> = RingElement::sample_poly_cbd(&bytes, 0xFF);
//         let ring_element_copy = ring_element;

//         // runs .ntt() on intstantiation
//         let ntt_element = NttElement::new(&mut ring_element);
//         assert_eq!(
//             ring_element_copy.coefs,
//             Into::<RingElement<P768>>::into(ntt_element).coefs
//         );
//     }

//     #[test]
//     fn test_multiply_ntts_associative() {
//         let bytes: Vec<u8> = (0..32)
//             .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
//             .collect();

//         let a: NttElement<P768> = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
//         let b = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xBB));
//         let c = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xCC));

//         // Test associativity (ab)c = a(bc)
//         let ab_c = (a * b) * c;
//         let a_bc = a * (b * c);
//         assert_eq!(ab_c.coefficients, a_bc.coefficients);

//         let a: NttElement<P768> = NttElement::sample_ntt(&bytes.clone(), 0, 0);
//         let b = NttElement::sample_ntt(&bytes.clone(), 0, 0);
//         let c = NttElement::sample_ntt(&bytes.clone(), 0, 0);

//         // Test associativity (ab)c = a(bc)
//         let ab_c = (a * b) * c;
//         let a_bc = a * (b * c);
//         assert_eq!(ab_c.coefficients, a_bc.coefficients);
//     }

//     #[test]
//     fn test_multiply_ntts_zero() {
//         let bytes: Vec<u8> = (0..32)
//             .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
//             .collect();

//         let a: NttElement<P768> = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
//         let zero = NttElement::zero();

//         // Test multiplicative identity
//         let res = zero * a;
//         assert_eq!(
//             res,
//             NttElement {
//                 coefficients: [F::zero(); 256]
//             }
//         );
//         let res = a * zero;
//         assert_eq!(
//             res,
//             NttElement {
//                 coefficients: [F::zero(); 256]
//             }
//         );
//     }

//     #[test]
//     fn test_closure_under_multiplication() {
//         for _ in 0..1000 {
//             let bytes: Vec<u8> = (0..32)
//                 .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
//                 .collect();

//             let a: NttElement<P768> =
//                 NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xAA));
//             let b = NttElement::new(&mut RingElement::sample_poly_cbd(&bytes, 0xBB));

//             let result = a * b;
//             assert!(
//                 result.coefficients.iter().all(|x| x.val() < 3329),
//                 "Result of multiplication must be valid NttElement"
//             );
//         }
//     }

//     #[test]
//     fn test_encode_decode() {
//         // Initialize a CRNG
//         let bytes: Vec<u8> = (0..384)
//             .map(|_| ChaCha20Rng::seed_from_u64(0x7FFFFFFFFFFFFFFF).gen())
//             .collect();

//         // Initialize a RingElement with random values
//         let a: NttElement<P768> = NttElement::sample_ntt(&bytes, 1, 1);

//         // Encode the RingElement to bytes
//         let encoded_bytes = {
//             let mut out = Vec::with_capacity(256 * 12 / 8);

//             for i in (0..a.coefficients.len()).step_by(2) {
//                 // Combine two 12-bit integers into a single 24-bit integer
//                 let x = u32::from(a.coefficients[i].val())
//                     | (u32::from(a.coefficients[i + 1].val()) << 12);

//                 // Split the 24-bit integer into 3 bytes and append to the output vector
//                 out.push((x & 0xFF) as u8); // First 8 bits
//                 out.push(((x >> 8) & 0xFF) as u8); // Next 8 bits
//                 out.push(((x >> 16) & 0xFF) as u8); // Last 8 bits
//             }
//             out
//         };

//         // Decode the bytes back into a vector of F
//         let decoded_elements: NttElement<P768> =
//             NttElement::byte_decode_12(&encoded_bytes).expect("Decoding failed");

//         // Verify the decoded vector matches the original RingElement's array
//         assert_eq!(
//             decoded_elements.coefficients.len(),
//             a.coefficients.len(),
//             "Length mismatch"
//         );

//         // Check each value for equality after encoding and decoding
//         for (decoded_elem, original_elem) in decoded_elements
//             .coefficients
//             .iter()
//             .zip(a.coefficients.iter())
//         {
//             assert_eq!(decoded_elem.val(), original_elem.val(), "Value mismatch");
//         }
//     }
// }
