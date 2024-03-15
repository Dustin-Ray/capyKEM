use core::fmt;
use core::ops::AddAssign;
use core::ops::{Add, Sub};

use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use alloc::vec::Vec;

use crate::constants::ml_kem_constants::N;
use crate::constants::parameter_sets::ParameterSet;
use crate::math::field_element::FieldElement as F;

/// A polynomial is an element of the ring R. It is an array of 256 coefficients
/// which themselves are [F].
#[derive(Clone, Copy)]
pub struct RingElement<P> {
    pub coefficients: [F<P>; 256],
}

impl<P: ParameterSet + Copy> RingElement<P> {
    // Create a new RingElement from a vector of FieldElements
    pub fn new(val: [F<P>; 256]) -> Self {
        RingElement { coefficients: val }
    }

    pub fn zero() -> Self {
        [F::new(0); 256].into()
    }

    // REMARKS:
    // make generic for NTT and Ring
    pub fn byte_encode(self) -> Vec<u8> {
        let mut out = Vec::with_capacity(256 * 12 / 8); // Preallocate the output vector

        for i in (0..self.coefficients.len()).step_by(2) {
            // Combine two 12-bit integers into a single 24-bit integer
            let x = u32::from(self.coefficients[i].val())
                | (u32::from(self.coefficients[i + 1].val()) << 12);

            // Split the 24-bit integer into 3 bytes and append to the output vector
            out.push((x & 0xFF) as u8); // First 8 bits
            out.push(((x >> 8) & 0xFF) as u8); // Next 8 bits
            out.push(((x >> 16) & 0xFF) as u8); // Last 8 bits
        }
        out
    }

    // REMARKS:
    // TODO: parameterize eta1 and eta 2
    pub fn sample_poly_cbd(s: &[u8], b: u8) -> RingElement<P> {
        let mut prf = Shake256::default();
        prf.update(s);
        prf.update(&[b]);

        let mut b = [0u8; (N / 2) as usize];
        let mut reader = prf.finalize_xof();
        reader.read(&mut b);

        let mut f = [F::new(0); N as usize];

        for i in (0..N).step_by(2) {
            // Iterate through indices, stepping by 2.
            let b = b[(i / 2) as usize];
            let b_7 = (b >> 7) & 1;
            let b_6 = (b >> 6) & 1;
            let b_5 = (b >> 5) & 1;
            let b_4 = (b >> 4) & 1;
            let b_3 = (b >> 3) & 1;
            let b_2 = (b >> 2) & 1;
            let b_1 = (b >> 1) & 1;
            let b_0 = b & 1;

            f[i as usize] = F::new((b_0 + b_1).into()) - F::new((b_2 + b_3).into());
            // Ensure i+1 doesn't go out of bounds, relevant if N is odd.
            if i + 1 < N {
                f[(i + 1) as usize] = F::new((b_4 + b_5).into()) - F::new((b_6 + b_7).into());
            }
        }
        RingElement::new(f)
    }
}

impl<P: ParameterSet + Copy> fmt::Debug for RingElement<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (index, element) in self.coefficients.iter().enumerate() {
            write!(f, "{:<8}", element.val())?;
            // Adjust for row width
            if (index + 1) % 8 == 0 {
                writeln!(f)?;
            }
        }
        Ok(())
    }
}

impl<P: ParameterSet + Copy> From<[F<P>; 256]> for RingElement<P> {
    fn from(val: [F<P>; 256]) -> Self {
        RingElement::new(val)
    }
}

impl<P: ParameterSet + Copy> AddAssign for RingElement<P> {
    fn add_assign(&mut self, other: Self) {
        for (lhs, rhs) in self.coefficients.iter_mut().zip(other.coefficients.iter()) {
            *lhs += *rhs;
        }
    }
}

impl<P: ParameterSet + Copy> Add for RingElement<P> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        assert_eq!(
            self.coefficients.len(),
            other.coefficients.len(),
            "RingElements must be of the same length"
        );

        let mut result = [F::zero(); 256];
        for (i, item) in self.coefficients.iter().enumerate().take(256) {
            result[i] = *item + other.coefficients[i];
        }
        RingElement::new(result)
    }
}

impl<P: ParameterSet + Copy> Sub for RingElement<P> {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        assert_eq!(
            self.coefficients.len(),
            other.coefficients.len(),
            "RingElements must be of the same length"
        );

        let mut result = [F::zero(); 256];
        for (i, item) in self.coefficients.iter().enumerate().take(256) {
            result[i] = *item - other.coefficients[i];
        }
        RingElement::new(result)
    }
}

impl<P: ParameterSet + Copy + PartialEq> PartialEq for RingElement<P> {
    fn eq(&self, other: &Self) -> bool {
        if self.coefficients.len() != other.coefficients.len() {
            return false;
        }
        self.coefficients
            .iter()
            .zip(other.coefficients.iter())
            .all(|(a, b)| a == b)
    }
}
