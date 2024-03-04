pub mod ml_kem_constants {
    pub const Q: u16 = 3329;
    pub const N: u16 = 256;
    pub const ENCODE_SIZE_12: u16 = N * 12 / 8;
}

pub mod barrett_constants {
    pub const MULTIPLIER: u16 = 5039; // 4¹² / q,
    pub const SHIFT: u16 = 24; // log₂(4¹²)
}
