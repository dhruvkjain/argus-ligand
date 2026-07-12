//! GC-rich region and CpG-island detection with a sliding window.
//!
//! GC content and CpG ratio are strand independent (complementing a strand
//! swaps C with G and leaves both counts and the palindromic CpG dinucleotide
//! unchanged), so this scan does not take a strand option.

use crate::types::GcMatch;

/// GC content and observed/expected CpG ratio of a window.
///
/// The CpG ratio uses the Gardiner-Garden and Frommer definition:
/// `observed CpG * length / (count C * count G)`.
///
/// # Arguments
///
/// * `win` - The window bytes to measure.
///
/// # Returns
///
/// A tuple `(gc_fraction, cpg_ratio)`. Both are 0.0 when the window is empty or
/// has no C and G bases.
fn window_stats(win: &[u8]) -> (f64, f64) {
    if win.is_empty() {
        return (0.0, 0.0);
    }
    let mut c = 0usize;
    let mut g = 0usize;
    let mut cpg = 0usize;
    for (i, &b) in win.iter().enumerate() {
        match b {
            b'C' => c += 1,
            b'G' => g += 1,
            _ => {}
        }
        if b == b'C' && win.get(i + 1) == Some(&b'G') {
            cpg += 1;
        }
    }
    let n = win.len() as f64;
    let gc = (c + g) as f64 / n;
    let expected = (c as f64 * g as f64) / n;
    let ratio = if expected > 0.0 {
        cpg as f64 / expected
    } else {
        0.0
    };
    (gc, ratio)
}

/// Find GC-rich regions and CpG islands with a sliding window.
///
/// Slides a window of size `min_len` across the sequence. A window qualifies
/// when its GC fraction is at least `min_gc` and its CpG ratio is at least
/// `min_cpg_ratio`. Overlapping qualifying windows are merged into regions, and
/// each region's stats are recomputed over its full span. Set `min_cpg_ratio`
/// to 0 to find any GC-rich region, or about 0.6 for classic CpG islands.
///
/// # Arguments
///
/// * `seq` - Cleaned nucleotide bytes.
/// * `min_len` - Window size and minimum region length.
/// * `min_gc` - Minimum GC fraction, 0.0 to 1.0.
/// * `min_cpg_ratio` - Minimum observed/expected CpG ratio.
///
/// # Returns
///
/// Merged regions in ascending position order. Empty when the sequence is
/// shorter than `min_len` or nothing qualifies.
pub(crate) fn find_gc_regions(
    seq: &[u8],
    min_len: usize,
    min_gc: f64,
    min_cpg_ratio: f64,
) -> Vec<GcMatch> {
    let win = min_len.max(1);
    if seq.len() < win {
        return Vec::new();
    }
    let mut regions: Vec<(usize, usize)> = Vec::new();
    for i in 0..=(seq.len() - win) {
        let (gc, ratio) = window_stats(&seq[i..i + win]);
        if gc >= min_gc && ratio >= min_cpg_ratio {
            let end = i + win;
            match regions.last_mut() {
                // Overlapping or adjacent window: extend the current region.
                Some(last) if i <= last.1 => last.1 = end,
                _ => regions.push((i, end)),
            }
        }
    }
    regions
        .into_iter()
        .map(|(start, end)| {
            let (gc, ratio) = window_stats(&seq[start..end]);
            GcMatch {
                start,
                end,
                length: end - start,
                gc_percent: gc * 100.0,
                cpg_ratio: ratio,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: &str) -> Vec<u8> {
        s.bytes().collect()
    }

    #[test]
    fn finds_a_gc_rich_region() {
        // A GC-rich stretch flanked by AT-rich sequence.
        let seq = v("AAAAAAAAAAGCGCGCGCGCGCGCGCAAAAAAAAAA");
        let regions = find_gc_regions(&seq, 10, 0.7, 0.0);
        assert_eq!(regions.len(), 1);
        assert!(regions[0].gc_percent > 70.0);
        // The region covers the GC core (positions 10..24).
        assert!(regions[0].start <= 10 && regions[0].end >= 24);
    }

    #[test]
    fn at_rich_sequence_has_no_gc_region() {
        let seq = v("ATATATATATATATATATAT");
        assert!(find_gc_regions(&seq, 10, 0.5, 0.0).is_empty());
    }

    #[test]
    fn cpg_threshold_excludes_low_cpg_region() {
        // 100% GC but no CpG dinucleotide (no C followed by G).
        let seq = v("GGGGGGGGGGCCCCCCCCCC");
        assert!(!find_gc_regions(&seq, 10, 0.5, 0.0).is_empty()); // GC-rich: yes
        assert!(find_gc_regions(&seq, 10, 0.5, 0.6).is_empty()); // CpG island: no
    }
}
