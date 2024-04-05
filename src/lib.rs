// #![no_std]
extern crate alloc;

use crate::constants::parameter_sets::ParameterSet;
use core::marker::PhantomData;

#[allow(non_upper_case_globals)]
mod constants;
#[allow(non_upper_case_globals)]
mod fips203;
#[allow(non_upper_case_globals)]
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
