pub mod ml_kem_constants {
    pub const q: u16 = 3329;
    pub const n: u16 = 256;
}

pub mod barrett_constants {
    pub const MULTIPLIER: u16 = 5039; // 4¹² / q,
    pub const SHIFT: u16 = 24; // log₂(4¹²)
}
