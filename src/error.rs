//! Error type for the motif engine.

use std::fmt;

/// An error produced while parsing a request or running a scan.
///
/// The engine turns these into a JSON `{"error": "..."}` object at the
/// boundary, so callers of [`scan`](crate::scan) never see the type directly.
/// It is public so that code using the engine as a library (for example a
/// future scanner worker) can match on the cause.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EngineError {
    /// The request body was not valid JSON or did not match the expected shape.
    InvalidRequest(String),
    /// The sequence was empty after cleaning.
    EmptySequence,
    /// A motif pattern contained no bases.
    EmptyPattern,
    /// A motif pattern contained a character that is not a valid IUPAC code.
    InvalidIupacCode(char),
}

impl fmt::Display for EngineError {
    /// Write a short, human readable description of the error.
    ///
    /// # Arguments
    ///
    /// * `f` - The formatter to write into.
    ///
    /// # Returns
    ///
    /// `Ok(())` once the message has been written.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EngineError::InvalidRequest(msg) => write!(f, "invalid request: {msg}"),
            EngineError::EmptySequence => write!(f, "sequence is empty after cleaning"),
            EngineError::EmptyPattern => write!(f, "pattern has no bases"),
            EngineError::InvalidIupacCode(c) => write!(f, "'{c}' is not a valid IUPAC code"),
        }
    }
}

impl std::error::Error for EngineError {}
