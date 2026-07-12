//! Sequence cleaning and strand operations.
//!
//! Everything downstream works on a clean byte slice of uppercase nucleotide
//! letters, so the cleaning step here is the single place that decides what
//! counts as valid input.

use crate::pattern::base_mask;

/// Which DNA strand a hit was found on.
///
/// The engine searches the forward strand directly and searches the minus
/// strand by scanning the reverse complement. This enum records which case a
/// result came from and provides the symbol used in the JSON output.
#[derive(Clone, Copy)]
pub(crate) enum Strand {
    /// The given forward strand.
    Forward,
    /// The reverse complement strand.
    Reverse,
}

impl Strand {
    /// Return the single character used for this strand in output: `'+'` or `'-'`.
    ///
    /// # Returns
    ///
    /// `'+'` for [`Strand::Forward`], `'-'` for [`Strand::Reverse`].
    pub(crate) fn symbol(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }
}

/// Clean raw input into a slice of uppercase nucleotide bytes.
///
/// Removes FASTA header lines (those starting with `>`) and comment lines
/// (starting with `;`), drops all whitespace, uppercases every letter, and maps
/// `U` to `T` so RNA input is accepted. Any remaining character that is not a
/// valid nucleotide or IUPAC code is dropped and counted in a warning.
///
/// # Arguments
///
/// * `raw` - The pasted sequence, either raw letters or FASTA text.
///
/// # Returns
///
/// A tuple of the cleaned bytes and a list of warnings. The warnings list is
/// empty when nothing was dropped.
pub(crate) fn clean_sequence(raw: &str) -> (Vec<u8>, Vec<String>) {
    let mut out = Vec::with_capacity(raw.len());
    let mut dropped = 0usize;
    for line in raw.lines() {
        if line.starts_with('>') || line.starts_with(';') {
            continue;
        }
        for ch in line.bytes() {
            if ch.is_ascii_whitespace() {
                continue;
            }
            let up = ch.to_ascii_uppercase();
            let up = if up == b'U' { b'T' } else { up };
            if base_mask(up).is_some() {
                out.push(up);
            } else {
                dropped += 1;
            }
        }
    }
    let mut warnings = Vec::new();
    if dropped > 0 {
        warnings.push(format!("dropped {dropped} non-nucleotide character(s)"));
    }
    (out, warnings)
}

/// Return the reverse complement of a sequence.
///
/// Complements each base using IUPAC rules, then reverses the order, which is
/// the standard way to read the opposite DNA strand.
///
/// # Arguments
///
/// * `seq` - Cleaned nucleotide bytes.
///
/// # Returns
///
/// A new vector holding the reverse complement.
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Return the IUPAC complement of a single base byte.
///
/// # Arguments
///
/// * `b` - One uppercase nucleotide or IUPAC code byte.
///
/// # Returns
///
/// The complement byte. An unrecognized byte is returned unchanged.
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_fasta_header_and_whitespace() {
        let (seq, warnings) = clean_sequence(">header\nACGT ACGT\nacgt");
        assert_eq!(seq, b"ACGTACGTACGT");
        assert!(warnings.is_empty());
    }

    #[test]
    fn maps_u_to_t_and_drops_junk() {
        let (seq, warnings) = clean_sequence("ACGU123");
        assert_eq!(seq, b"ACGT");
        assert_eq!(warnings.len(), 1);
    }

    #[test]
    fn reverse_complement_is_correct() {
        assert_eq!(reverse_complement(b"AACGTT"), b"AACGTT");
        assert_eq!(reverse_complement(b"GAATTC"), b"GAATTC");
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
    }
}
