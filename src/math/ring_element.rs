use crate::constants::ml_kem_constants::{n, q, ENCODE_10};
use crate::math::field_element::FieldElement as F;
use crate::ParameterSet;
use alloc::vec::Vec;
use core::fmt;
use core::iter::Sum;
use core::ops::AddAssign;
use core::ops::{Add, Sub};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use hybrid_array::Array;
use hybrid_array::ArrayOps;
use typenum::{Unsigned, U0, U6};

/// A polynomial is an element of the ring R. It is an array of 256 coefficients
/// which themselves are [F].
#[derive(Clone, Copy)]
pub struct RingElement<P> {
    pub coefs: [F<P>; n],
}

impl<P: ParameterSet + Copy> RingElement<P> {
    // Create a new RingElement from a vector of FieldElements
    pub fn new(val: [F<P>; n]) -> Self {
        RingElement { coefs: val }
    }

    pub fn zero() -> Self {
        [F::new(0); n].into()
    }

    pub fn decode_decompress_1(b: &[u8]) -> Result<Self, &'static str> {
        // TODO: check fips203 for length check requirement
        let coefs: [F<P>; n] = (0..n)
            .map(|i| {
                let bit = (b[i / 8] >> (i % 8)) & 1; // Extract the i-th bit
                F::new(bit as u16 * q / 2)
            })
            .collect::<Vec<F<P>>>()
            .try_into()
            .unwrap_or_else(|_| panic!("Incorrect vector size, expected 256 elements"));

        Ok(RingElement { coefs })
    }

    pub fn decode_and_decompress_10(b: &[u8]) -> Result<Self, &'static str> {
        let mut f = RingElement {
            coefs: [F::zero(); n],
        };
        let mut x: u64;
        let mut slice = b;

        for i in (0..n).step_by(4) {
            x = u64::from(slice[0])
                | u64::from(slice[1]) << 8
                | u64::from(slice[2]) << 16
                | u64::from(slice[3]) << 24
                | u64::from(slice[4]) << 32;
            slice = &slice[5..]; // Move the slice window

            f.coefs[i] = F::decompress::<10>((x & 0x3FF) as u16);
            f.coefs[i + 1] = F::decompress::<10>((x >> 10 & 0x3FF) as u16);
            f.coefs[i + 2] = F::decompress::<10>((x >> 20 & 0x3FF) as u16);
            f.coefs[i + 3] = F::decompress::<10>((x >> 30 & 0x3FF) as u16);
        }

        Ok(f)
    }

    pub fn decode_and_decompress_4(b: &[u8]) -> Result<Self, &'static str> {
        if b.len() != 128 {
            return Err("invalid encoding length");
        }

        let mut f = RingElement {
            coefs: [F::zero(); n],
        }; // Assuming n and F::zero() are defined
        let mut index = 0;

        for chunk in b.iter() {
            if index >= f.coefs.len() {
                break; // Prevents out-of-bounds access if n is not exactly double b's length
            }

            // Decompress the lower 4 bits
            f.coefs[index] = F::decompress::<4>(((*chunk) & 0x0F) as u16);
            index += 1;

            // Check if we should also process the upper 4 bits
            if index < f.coefs.len() {
                // Decompress the upper 4 bits
                f.coefs[index] = F::decompress::<4>(((*chunk) >> 4) as u16);
                index += 1;
            }
        }

        Ok(f)
    }

    pub fn compress_and_encode_10(mut s: Vec<u8>, f: Self) -> Vec<u8> {
        s.reserve(ENCODE_10);

        for i in (0..f.coefs.len()).step_by(4) {
            let mut x: u64 = 0;
            for j in 0..4 {
                if i + j < f.coefs.len() {
                    let shift = j * 10;
                    x |= (F::<P>::compress::<10>(&f.coefs[i + j]) as u64) << shift;
                }
            }
            for shift in (0..40).step_by(8) {
                s.push(((x >> shift) & 0xFF) as u8);
            }
        }

        s
    }

    pub fn compress_and_encode_4(mut s: Vec<u8>, f: Self) -> Vec<u8> {
        s.reserve(128);

        for i in (0..n).step_by(2) {
            let compressed_pair = F::<P>::compress::<4>(&f.coefs[i]) as u8
                | (F::<P>::compress::<4>(&f.coefs[i + 1]) as u8) << 4;
            s.push(compressed_pair);
        }
        s
    }

    pub fn compress_and_encode_1(&self, s: &mut [u8]) {
        for (i, coef) in self.coefs.iter().enumerate() {
            let compressed = F::<P>::compress::<1>(coef) as u8;
            let byte_index = i / 8;
            let bit_index = i % 8;

            s[byte_index] |= compressed << bit_index;
        }
    }
    // TODO: make eta const/generic
    pub fn sample_poly_cbd(s: &[u8], b: u8) -> RingElement<P> {
        let mut prf = Shake256::default();
        prf.update(s);
        prf.update(&[b]);

        // this should be 64 * eta
        let mut b = [0u8; 128];
        let mut reader = prf.finalize_xof();
        reader.read(&mut b);

        let mut f = [F::new(0); n];

        for i in (0..n).step_by(2) {
            // Iterate through indices, stepping by 2.
            let b = b[i / 2];
            let b_7 = (b >> 7) & 1;
            let b_6 = (b >> 6) & 1;
            let b_5 = (b >> 5) & 1;
            let b_4 = (b >> 4) & 1;
            let b_3 = (b >> 3) & 1;
            let b_2 = (b >> 2) & 1;
            let b_1 = (b >> 1) & 1;
            let b_0 = b & 1;

            f[i] = F::new((b_0 + b_1).into()) - F::new((b_2 + b_3).into());
            // Ensure i+1 doesn't go out of bounds, relevant if N is odd.
            if i + 1 < n {
                f[i + 1] = F::new((b_4 + b_5).into()) - F::new((b_6 + b_7).into());
            }
        }
        RingElement::new(f)
    }
}

impl<P: ParameterSet + Copy> fmt::Debug for RingElement<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (index, element) in self.coefs.iter().enumerate() {
            write!(f, "{:<8}", element.val())?;
            // Adjust for row width
            if (index + 1) % 8 == 0 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

impl<P: ParameterSet + Copy> From<[F<P>; n]> for RingElement<P> {
    fn from(val: [F<P>; n]) -> Self {
        RingElement::new(val)
    }
}

impl<P: ParameterSet + Copy> AddAssign for RingElement<P> {
    fn add_assign(&mut self, other: Self) {
        for (lhs, rhs) in self.coefs.iter_mut().zip(other.coefs.iter()) {
            *lhs += *rhs;
        }
    }
}

impl<P: ParameterSet + Copy> Add for RingElement<P> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        assert_eq!(
            self.coefs.len(),
            other.coefs.len(),
            "RingElements must be of the same length"
        );

        let mut result = [F::zero(); n];
        for (i, item) in self.coefs.iter().enumerate().take(n) {
            result[i] = *item + other.coefs[i];
        }
        RingElement::new(result)
    }
}

impl<P: ParameterSet + Copy> Sum for RingElement<P> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(RingElement::<P>::zero(), |acc, elem| acc + elem)
    }
}

impl<P: ParameterSet + Copy> Sub for RingElement<P> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        assert_eq!(
            self.coefs.len(),
            other.coefs.len(),
            "RingElements must be of the same length"
        );

        let mut result = [F::zero(); n];
        for (i, item) in self.coefs.iter().enumerate().take(n) {
            result[i] = *item - other.coefs[i];
        }
        RingElement::new(result)
    }
}

impl<P: ParameterSet + Copy + PartialEq> PartialEq for RingElement<P> {
    fn eq(&self, other: &Self) -> bool {
        if self.coefs.len() != other.coefs.len() {
            return false;
        }
        self.coefs
            .iter()
            .zip(other.coefs.iter())
            .all(|(a, b)| a == b)
    }
}
