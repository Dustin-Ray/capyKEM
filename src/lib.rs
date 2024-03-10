#[allow(dead_code)]
#[allow(non_upper_case_globals)]
mod constants;
#[allow(dead_code)]
mod fips203;
#[allow(dead_code)]
mod math;

pub struct Message {
    pub data: Vec<u8>,
    pub k: u16,
    pub du: u16,
    pub dv: u16,
    pub eta_1: u16,
    pub eta_2: u16,
}
