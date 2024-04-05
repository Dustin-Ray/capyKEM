use crate::constants::ml_kem_constants::n;
use crate::math::field_element::FieldElement as F;

use core::fmt;
use core::iter::Sum;
use core::ops::AddAssign;
use core::ops::{Add, Sub};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

/// A polynomial is an element of the ring R. It is an array of 256 coefficients
/// which themselves are [F].
#[derive(Clone, Copy)]
pub struct RingElement {
    pub coefs: [F; n],
}

impl RingElement {
    // Create a new RingElement from a vector of FieldElements
    pub fn new(val: [F; n]) -> Self {
        RingElement { coefs: val }
    }

    pub fn zero() -> Self {
        [F::new(0); n].into()
    }

    // TODO: make eta const/generic
    pub fn sample_poly_cbd(s: &[u8], b: u8) -> RingElement {
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

impl fmt::Debug for RingElement {
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

impl From<[F; n]> for RingElement {
    fn from(val: [F; n]) -> Self {
        RingElement::new(val)
    }
}

impl AddAssign for RingElement {
    fn add_assign(&mut self, other: Self) {
        for (lhs, rhs) in self.coefs.iter_mut().zip(other.coefs.iter()) {
            *lhs += *rhs;
        }
    }
}

impl Add for RingElement {
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

impl Sum for RingElement {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(RingElement::zero(), |acc, elem| acc + elem)
    }
}

impl Sub for RingElement {
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

impl PartialEq for RingElement {
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
