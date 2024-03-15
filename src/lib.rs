// #![no_std]
extern crate alloc;
use constants::parameter_sets::ParameterSet;
use core::marker::PhantomData;

#[allow(dead_code)]
#[allow(non_upper_case_globals)]
mod constants;
#[allow(dead_code)]
mod fips203;
#[allow(dead_code)]
pub mod math;

/// A container for the 32-byte shared secret key
/// that will be exchanged via encapsulation/decapsulation.
/// It is defined over generic parameter set defined
/// by NIST of various security levels.
#[derive(Debug)]
pub struct Secret<P> {
    pub m: [u8; 32],
    _marker: PhantomData<P>,
}

impl<P: ParameterSet> Secret<P> {
    #[must_use]
    pub fn new(data: [u8; 32]) -> Self {
        Secret::<P> {
            m: data,
            _marker: PhantomData,
        }
    }
}
