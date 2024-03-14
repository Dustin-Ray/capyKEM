use core::marker::PhantomData;

use constants::parameter_sets::ParameterSet;

#[allow(dead_code)]
#[allow(non_upper_case_globals)]
mod constants;
#[allow(dead_code)]
mod fips203;
#[allow(dead_code)]
mod math;

pub struct Message<P> {
    pub m: [u8; 32],
    _marker: PhantomData<P>,
}

impl<P: ParameterSet> Message<P> {
    // Constructor function
    pub fn new(data: [u8; 32]) -> Self {
        Message::<P> {
            m: data,
            _marker: PhantomData,
        }
    }
}
