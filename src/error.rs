use core::fmt;

/// Generic error type for ML-KEM operations
/// 
/// Error messages are intentionally generic to avoid leaking information
/// that could be used in timing or other side-channel attacks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KemError {
    /// Invalid input provided to a function
    InvalidInput,
    /// Decapsulation failed
    DecapsulationFailure,
    /// Encoding/decoding error
    EncodingError,
}

impl fmt::Display for KemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            KemError::InvalidInput => write!(f, "Invalid input"),
            KemError::DecapsulationFailure => write!(f, "Decapsulation failed"),
            KemError::EncodingError => write!(f, "Encoding error"),
        }
    }
}

// Note: no_std mode - no std::error::Error implementation

/// Type alias for Results using KemError
pub type Result<T> = core::result::Result<T, KemError>;

