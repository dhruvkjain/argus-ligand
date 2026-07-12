//! Motif search over compiled IUPAC masks.

use crate::pattern::iupac::{base_mask, Mask};
use crate::sequence::Strand;
use crate::types::MotifMatch;

/// Find every position where a compiled pattern matches a sequence.
///
/// `hay` is the sequence to search, already oriented on `strand`. A position
/// matches when, for every pattern position, the sequence base is allowed by
/// the pattern mask. Ambiguous sequence bases only match when their whole set
/// of possibilities fits inside the pattern mask, so they do not over match.
///
/// To search the minus strand, pass the reverse complement as `hay` and the
/// original forward length as `fwd_len`. The returned coordinates are mapped
/// back to the forward strand so all results share one coordinate system.
///
/// # Arguments
///
/// * `hay` - Sequence bytes to search, oriented on `strand`.
/// * `masks` - Compiled pattern, one mask per position, from [`crate::pattern::iupac::compile_iupac`].
/// * `strand` - Which strand `hay` represents.
/// * `fwd_len` - Length of the forward strand, used to map minus strand hits.
///
/// # Returns
///
/// A list of matches in forward strand coordinates, in ascending position
/// order. Empty when the pattern is longer than the sequence or nothing matches.
pub(crate) fn find_motif(
    hay: &[u8],
    masks: &[Mask],
    strand: Strand,
    fwd_len: usize,
) -> Vec<MotifMatch> {
    let mut hits = Vec::new();
    let plen = masks.len();
    if hay.len() < plen {
        return hits;
    }
    for i in 0..=(hay.len() - plen) {
        if !matches_at(hay, masks, i) {
            continue;
        }
        let matched = String::from_utf8_lossy(&hay[i..i + plen]).into_owned();
        let (start, end) = match strand {
            Strand::Reverse => (fwd_len - (i + plen), fwd_len - i),
            Strand::Forward => (i, i + plen),
        };
        hits.push(MotifMatch {
            start,
            end,
            strand: strand.symbol(),
            matched,
        });
    }
    hits
}

/// Strand symbol for a palindromic site that matches on both strands at once.
pub(crate) const PALINDROME: char = '±';

/// Search one or both strands and merge palindromic hits.
///
/// Runs [`find_motif`] on the forward strand, and if `both_strands` is set, on
/// the reverse complement too. A palindromic pattern (one that reads the same on
/// both strands, like the EcoRI site `GAATTC`) produces a forward hit and a
/// reverse hit at the exact same coordinates. Those describe one physical site,
/// so they are merged into a single match whose strand is [`PALINDROME`] (`±`)
/// instead of being counted twice.
///
/// # Arguments
///
/// * `fwd` - The forward strand sequence.
/// * `rc` - The reverse complement of `fwd`, precomputed by the caller.
/// * `masks` - Compiled pattern from [`crate::pattern::iupac::compile_iupac`].
/// * `both_strands` - Whether to also search the reverse complement.
///
/// # Returns
///
/// Matches in forward strand coordinates, sorted by position, with palindromic
/// sites collapsed to one `±` entry each.
pub(crate) fn find_motif_stranded(
    fwd: &[u8],
    rc: &[u8],
    masks: &[Mask],
    both_strands: bool,
) -> Vec<MotifMatch> {
    let mut hits = find_motif(fwd, masks, Strand::Forward, fwd.len());
    if both_strands {
        for m in find_motif(rc, masks, Strand::Reverse, fwd.len()) {
            match hits.iter_mut().find(|e| e.start == m.start && e.end == m.end) {
                // Same span on both strands means the site is palindromic.
                Some(existing) => existing.strand = PALINDROME,
                None => hits.push(m),
            }
        }
        hits.sort_by(|a, b| a.start.cmp(&b.start).then(a.strand.cmp(&b.strand)));
    }
    hits
}

/// Test whether the pattern matches `hay` starting at index `at`.
///
/// # Arguments
///
/// * `hay` - Sequence bytes to test against.
/// * `masks` - Compiled pattern.
/// * `at` - Start index in `hay`. The caller guarantees the pattern fits.
///
/// # Returns
///
/// `true` if every pattern position is satisfied, `false` otherwise.
fn matches_at(hay: &[u8], masks: &[Mask], at: usize) -> bool {
    for (j, &m) in masks.iter().enumerate() {
        match base_mask(hay[at + j]) {
            // The sequence base's possibility set must fit inside the pattern set.
            Some(sm) if sm & !m == 0 => {}
            _ => return false,
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pattern::iupac::compile_iupac;

    fn v(s: &str) -> Vec<u8> {
        s.bytes().collect()
    }

    #[test]
    fn literal_match_forward() {
        let masks = compile_iupac("GAATTC").unwrap();
        let hits = find_motif(&v("AAAGAATTCAAA"), &masks, Strand::Forward, 12);
        assert_eq!(hits.len(), 1);
        assert_eq!((hits[0].start, hits[0].end), (3, 9));
        assert_eq!(hits[0].strand, '+');
    }

    #[test]
    fn iupac_tata_box_matches_variants() {
        // TATAWAWR catches TATAAAA and TATATAA that a literal TATAAA misses.
        let masks = compile_iupac("TATAWAWR").unwrap();
        assert_eq!(find_motif(&v("CCTATAAAAGG"), &masks, Strand::Forward, 11).len(), 1);
        assert_eq!(find_motif(&v("CCTATATAAGG"), &masks, Strand::Forward, 11).len(), 1);
        assert_eq!(find_motif(&v("CCTATGGGCGG"), &masks, Strand::Forward, 11).len(), 0);
    }

    #[test]
    fn palindromic_site_counted_once() {
        // GAATTC is its own reverse complement: one site, marked palindromic.
        let seq = v("AAAGAATTCAAA");
        let rc = crate::sequence::reverse_complement(&seq);
        let masks = compile_iupac("GAATTC").unwrap();
        let hits = find_motif_stranded(&seq, &rc, &masks, true);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, PALINDROME);
    }

    #[test]
    fn non_palindromic_keeps_separate_strands() {
        // ATG at position 0 on the forward strand; CAT at 5..8 reads ATG on the
        // minus strand. Different spans, so both are kept.
        let seq = v("ATGCCCAT");
        let rc = crate::sequence::reverse_complement(&seq);
        let masks = compile_iupac("ATG").unwrap();
        let hits = find_motif_stranded(&seq, &rc, &masks, true);
        assert_eq!(hits.len(), 2);
        assert!(hits.iter().any(|h| h.strand == '+'));
        assert!(hits.iter().any(|h| h.strand == '-'));
    }

    #[test]
    fn minus_strand_coords_map_back() {
        // The EcoRI site is its own reverse complement, so a forward hit implies
        // a minus strand hit at the same coordinates.
        let seq = v("AAAGAATTCAAA");
        let masks = compile_iupac("GAATTC").unwrap();
        let rc = crate::sequence::reverse_complement(&seq);
        let hits = find_motif(&rc, &masks, Strand::Reverse, seq.len());
        assert_eq!(hits.len(), 1);
        assert_eq!((hits[0].start, hits[0].end), (3, 9));
        assert_eq!(hits[0].strand, '-');
    }
}
