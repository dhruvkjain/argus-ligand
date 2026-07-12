//! Open reading frame detection and protein translation.

use crate::sequence::Strand;
use crate::types::OrfMatch;

/// Find open reading frames in all three frames of a sequence.
///
/// An ORF here is the span from the first `ATG` start codon to the next stop
/// codon (`TAA`, `TAG`, or `TGA`) in the same frame. Frames do not nest: after
/// a stop, the search looks for the next start. This runs over the three frames
/// of the given strand only. To cover all six frames, call it once per strand.
///
/// Coordinates are mapped back to the forward strand, so a minus strand ORF is
/// reported in forward strand positions.
///
/// # Arguments
///
/// * `seq` - Sequence bytes to search, oriented on `strand`.
/// * `strand` - Which strand `seq` represents.
/// * `min_aa` - Minimum protein length in amino acids, excluding the stop.
///   Shorter ORFs are skipped.
/// * `fwd_len` - Length of the forward strand, used to map minus strand ORFs.
///
/// # Returns
///
/// A list of ORFs in forward strand coordinates. Order follows frame then
/// position, so it is not globally sorted by position.
pub(crate) fn find_orfs(
    seq: &[u8],
    strand: Strand,
    min_aa: usize,
    fwd_len: usize,
) -> Vec<OrfMatch> {
    let mut orfs = Vec::new();
    for frame in 0..3 {
        let mut open_at: Option<usize> = None;
        let mut i = frame;
        while i + 3 <= seq.len() {
            let codon = &seq[i..i + 3];
            if open_at.is_none() && codon == b"ATG" {
                open_at = Some(i);
            } else if let Some(start) = open_at {
                if is_stop(codon) {
                    let end = i + 3; // include the stop codon
                    let aa_length = (i - start) / 3; // residues before the stop
                    if aa_length >= min_aa {
                        orfs.push(build_orf(seq, start, end, strand, frame, aa_length, fwd_len));
                    }
                    open_at = None;
                }
            }
            i += 3;
        }
    }
    orfs
}

/// Build one [`OrfMatch`], translating the protein and mapping coordinates.
///
/// # Arguments
///
/// * `seq` - The strand oriented sequence the ORF was found in.
/// * `start` - Start index of the ORF in `seq`, at the start codon.
/// * `end` - End index of the ORF in `seq`, just past the stop codon.
/// * `strand` - Which strand `seq` represents.
/// * `frame` - Reading frame 0, 1, or 2.
/// * `aa_length` - Protein length in amino acids, excluding the stop.
/// * `fwd_len` - Forward strand length, used to map minus strand coordinates.
///
/// # Returns
///
/// A fully populated [`OrfMatch`] in forward strand coordinates.
fn build_orf(
    seq: &[u8],
    start: usize,
    end: usize,
    strand: Strand,
    frame: usize,
    aa_length: usize,
    fwd_len: usize,
) -> OrfMatch {
    let protein = translate(&seq[start..end - 3]);
    let (fwd_start, fwd_end) = match strand {
        Strand::Reverse => (fwd_len - end, fwd_len - start),
        Strand::Forward => (start, end),
    };
    OrfMatch {
        start: fwd_start,
        end: fwd_end,
        strand: strand.symbol(),
        frame,
        aa_length,
        protein,
    }
}

/// Test whether a codon is a stop codon.
///
/// # Arguments
///
/// * `codon` - Exactly three nucleotide bytes.
///
/// # Returns
///
/// `true` for `TAA`, `TAG`, or `TGA`, `false` otherwise.
fn is_stop(codon: &[u8]) -> bool {
    matches!(codon, b"TAA" | b"TAG" | b"TGA")
}

/// Translate whole codons to a protein string using the standard genetic code.
///
/// Reads the input three bytes at a time. A trailing partial codon is ignored.
/// A codon that contains an ambiguity code translates to `X`.
///
/// # Arguments
///
/// * `seq` - Nucleotide bytes, read from the start in frame.
///
/// # Returns
///
/// The protein as a string, one letter per amino acid.
fn translate(seq: &[u8]) -> String {
    let mut protein = String::with_capacity(seq.len() / 3);
    let mut i = 0;
    while i + 3 <= seq.len() {
        protein.push(codon_to_aa(&seq[i..i + 3]));
        i += 3;
    }
    protein
}

/// Map one codon to its amino acid letter under the standard genetic code.
///
/// # Arguments
///
/// * `c` - Exactly three nucleotide bytes.
///
/// # Returns
///
/// The one letter amino acid code, `*` for a stop, or `X` for a codon that
/// contains an ambiguity code.
fn codon_to_aa(c: &[u8]) -> char {
    match c {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        _ => 'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: &str) -> Vec<u8> {
        s.bytes().collect()
    }

    #[test]
    fn finds_basic_orf() {
        // ATG AAA AAA TAA translates to MKK, length 3.
        let seq = v("ATGAAAAAATAA");
        let orfs = find_orfs(&seq, Strand::Forward, 1, seq.len());
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].protein, "MKK");
        assert_eq!(orfs[0].aa_length, 3);
        assert_eq!((orfs[0].start, orfs[0].end), (0, 12));
    }

    #[test]
    fn respects_min_length() {
        let seq = v("ATGAAAAAATAA"); // protein length 3
        assert_eq!(find_orfs(&seq, Strand::Forward, 4, seq.len()).len(), 0);
        assert_eq!(find_orfs(&seq, Strand::Forward, 3, seq.len()).len(), 1);
    }

    #[test]
    fn translate_handles_ambiguity() {
        assert_eq!(translate(b"ATGNNN"), "MX");
    }
}
