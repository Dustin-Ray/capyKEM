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

#[cfg(test)]
mod tests {
    use crate::constants::{K_INVERSE_NTT_ROOTS, K_MOD_ROOTS, K_NTT_ROOTS};

    use super::*;

    const Q: u16 = 3329;
    #[test]
    fn create_and_test_kntt_roots() {
        let kntt_roots: Vec<u16> = (0..128)
            .map(|i| mod_pow(17, bitreverse(i).into(), Q.into()))
            .collect();
        assert_eq!(kntt_roots, K_NTT_ROOTS);
    }

    #[test]
    fn create_and_test_k_inverse_ntt_roots() {
        let k_inverse_ntt_roots: Vec<u16> = (0..128)
            .map(|i| {
                let exp = Q - 1 - bitreverse(i);
                mod_pow(17, exp.into(), Q.into())
            })
            .collect();

        assert_eq!(k_inverse_ntt_roots, K_INVERSE_NTT_ROOTS);
    }

    #[test]
    fn create_and_test_k_mod_roots() {
        let k_mod_roots: Vec<u16> = (0..128)
            .map(|i| mod_pow(17, (2 * bitreverse(i) + 1).into(), Q.into()))
            .collect();

        assert_eq!(k_mod_roots, K_MOD_ROOTS);
    }
}
