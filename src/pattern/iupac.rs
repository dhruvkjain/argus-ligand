//! IUPAC nucleotide codes as bitmasks.
//!
//! Each base or ambiguity code is turned into a 4 bit mask over the set
//! {A, C, G, T}. Matching one position is then a single bitwise test, which
//! keeps the motif search loop fast. For example `N` (any base) is `1111` and
//! `W` (A or T) is `1001`.

use crate::error::EngineError;

/// A 4 bit set over {A, C, G, T}. Bit 0 is A, bit 1 is C, bit 2 is G, bit 3 is T.
pub(crate) type Mask = u8;

const A: Mask = 0b0001;
const C: Mask = 0b0010;
const G: Mask = 0b0100;
const T: Mask = 0b1000;

/// Map a single uppercase nucleotide byte to its IUPAC bitmask.
///
/// Accepts the four bases and all standard IUPAC ambiguity codes. This is also
/// used to decide whether a character is a valid nucleotide letter at all: a
/// return of `None` means it is not.
///
/// # Arguments
///
/// * `b` - One uppercase ASCII byte. `U` is not accepted here; callers map it
///   to `T` before calling.
///
/// # Returns
///
/// `Some(mask)` for a valid base or IUPAC code, or `None` for anything else.
pub(crate) fn base_mask(b: u8) -> Option<Mask> {
    Some(match b {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' => T,
        b'R' => A | G,
        b'Y' => C | T,
        b'S' => G | C,
        b'W' => A | T,
        b'K' => G | T,
        b'M' => A | C,
        b'B' => C | G | T,
        b'D' => A | G | T,
        b'H' => A | C | T,
        b'V' => A | C | G,
        b'N' => A | C | G | T,
        _ => return None,
    })
}

/// Compile a motif pattern string into one bitmask per position.
///
/// Input is uppercased and `U` is mapped to `T` before lookup. Whitespace
/// inside the pattern is ignored so patterns can be pasted with spaces.
///
/// # Arguments
///
/// * `pattern` - A motif in IUPAC codes, for example `TATAWAWR`.
///
/// # Returns
///
/// A vector of masks, one per pattern position, ready for [`crate::pattern::motif::find_motif`].
///
/// # Errors
///
/// Returns [`EngineError::EmptyPattern`] if the pattern has no bases.
/// Returns [`EngineError::InvalidIupacCode`] if a character is not a valid code.
pub(crate) fn compile_iupac(pattern: &str) -> Result<Vec<Mask>, EngineError> {
    let mut masks = Vec::with_capacity(pattern.len());
    for ch in pattern.bytes() {
        if ch.is_ascii_whitespace() {
            continue;
        }
        let up = ch.to_ascii_uppercase();
        let up = if up == b'U' { b'T' } else { up };
        match base_mask(up) {
            Some(m) => masks.push(m),
            None => return Err(EngineError::InvalidIupacCode(up as char)),
        }
    }
    if masks.is_empty() {
        return Err(EngineError::EmptyPattern);
    }
    Ok(masks)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compiles_valid_pattern() {
        let masks = compile_iupac("TATAWAWR").unwrap();
        assert_eq!(masks.len(), 8);
        assert_eq!(masks[4], A | T); // W
        assert_eq!(masks[7], A | G); // R
    }

    #[test]
    fn whitespace_is_ignored() {
        assert_eq!(compile_iupac("GA AT TC").unwrap().len(), 6);
    }

    #[test]
    fn rejects_invalid_code() {
        assert_eq!(
            compile_iupac("TAXQ"),
            Err(EngineError::InvalidIupacCode('X'))
        );
    }

    #[test]
    fn rejects_empty_pattern() {
        assert_eq!(compile_iupac("   "), Err(EngineError::EmptyPattern));
    }
}
