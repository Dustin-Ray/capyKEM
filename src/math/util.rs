/// Safely truncate an unsigned integer value to shorter representation
pub trait Truncate<T> {
    fn truncate(self) -> T;
}

macro_rules! define_truncate {
    ($from:ident, $to:ident) => {
        impl Truncate<$to> for $from {
            fn truncate(self) -> $to {
                // This line is marked unsafe because the `unwrap_unchecked` call is UB when its
                // `self` argument is `Err`.  It never will be, because we explicitly zeroize the
                // high-order bits before converting.  We could have used `unwrap()`, but chose to
                // avoid the possibility of panic.
                unsafe { (self & $from::from($to::MAX)).try_into().unwrap_unchecked() }
            }
        }
    };
}

define_truncate!(u32, u16);
define_truncate!(u64, u32);
define_truncate!(usize, u8);
define_truncate!(u128, u16);
define_truncate!(u128, u8);

#[cfg(test)]
mod tests {

    extern crate alloc;
    use crate::constants::{K_MOD_ROOTS, K_NTT_ROOTS};
    use alloc::vec::Vec;
    const Q: u16 = 3329;

    fn bitreverse(mut i: u16) -> u16 {
        let mut ret = 0;
        for _ in 0..7 {
            let bit = i & 1;
            ret <<= 1;
            ret |= bit;
            i >>= 1;
        }
        ret
    }

    fn mod_pow(base: u32, exp: u32, modulus: u32) -> u16 {
        let mut result = 1;
        let mut base = base % modulus;
        let mut exp = exp;

        while exp > 0 {
            if exp % 2 == 1 {
                result = (result * base) % modulus;
            }
            exp >>= 1;
            base = (base * base) % modulus;
        }

        result as u16
    }

    #[test]
    fn create_and_test_kntt_roots() {
        let kntt_roots: Vec<u16> = (0..128)
            .map(|i| mod_pow(17, bitreverse(i).into(), Q.into()))
            .collect();
        assert_eq!(kntt_roots, K_NTT_ROOTS);
    }

    // #[test]
    // fn create_and_test_k_inverse_ntt_roots() {
    //     let k_inverse_ntt_roots: Vec<u16> = (0..128)
    //         .map(|i| {
    //             let exp = Q - 1 - bitreverse(i);
    //             mod_pow(17, exp.into(), Q.into())
    //         })
    //         .collect();

    //     assert_eq!(k_inverse_ntt_roots, K_INVERSE_NTT_ROOTS);
    // }

    #[test]
    fn create_and_test_k_mod_roots() {
        let k_mod_roots: Vec<u16> = (0..128)
            .map(|i| mod_pow(17, (2 * bitreverse(i) + 1).into(), Q.into()))
            .collect();

        assert_eq!(k_mod_roots, K_MOD_ROOTS);
    }
}
