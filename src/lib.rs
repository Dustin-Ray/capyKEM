#![no_std]
extern crate alloc;
extern crate rand;

use alloc::vec::Vec;

use crate::constants::parameter_sets::ParameterSet;
use core::marker::PhantomData;

#[allow(non_upper_case_globals)]
pub mod constants;
#[allow(non_upper_case_globals)]
mod fips203;
#[allow(non_upper_case_globals)]
pub mod math;

/// A container for the 32-byte shared secret key
/// that will be exchanged via encapsulation/decapsulation.
/// It is defined over generic parameter set defined
/// by NIST of various security levels.
#[derive(Clone, Debug)]
pub struct Secret<P> {
    pub s: [u8; 32],
    _marker: PhantomData<P>,
}

impl<P: ParameterSet> Secret<P> {
    #[must_use]
    pub fn new(data: [u8; 32]) -> Self {
        Secret::<P> {
            s: data,
            _marker: PhantomData,
        }
    }
}

pub trait Encrypt {
    type Output;

    fn encrypt(&self, key: &[u8], rand: &[u8; 32]) -> Self::Output;
}

pub trait Decrypt {
    type Output;

    fn decrypt(&self, key: &[u8], ciphertext: &[u8]) -> Self::Output;
}

impl<P: ParameterSet + Copy> Encrypt for Secret<P> {
    type Output = Vec<u8>;

    fn encrypt(&self, ek_pke: &[u8], rand: &[u8; 32]) -> Self::Output {
        self.k_pke_encrypt(ek_pke, rand)
    }
}

impl<P: ParameterSet + Copy> Decrypt for Secret<P> {
    type Output = Vec<u8>;

    fn decrypt(&self, key: &[u8], ciphertext: &[u8]) -> Self::Output {
        self.k_pke_decrypt(key, ciphertext)
    }
}
