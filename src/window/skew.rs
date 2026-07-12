//! GC skew, used to locate bacterial replication origins.

/// Locate the replication origin and terminus from cumulative GC skew.
///
/// GC skew is `(G - C) / (G + C)` measured over a sliding window. The running
/// sum of the per window skew tends to reach its minimum near the replication
/// origin and its maximum near the terminus in many bacterial genomes.
///
/// # Arguments
///
/// * `seq` - Cleaned nucleotide bytes.
/// * `window` - Sliding window size. Clamped to the sequence length.
///
/// # Returns
///
/// `Some((origin, terminus, min_cumulative, max_cumulative))` where `origin` and
/// `terminus` are window start positions. `None` if the sequence is empty.
pub(crate) fn gc_skew_landmarks(seq: &[u8], window: usize) -> Option<(usize, usize, f64, f64)> {
    if seq.is_empty() {
        return None;
    }
    // Clamp so an oversized window falls back to the whole sequence.
    let win = window.clamp(1, seq.len());
    let mut cumulative = 0.0f64;
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;
    let mut min_pos = 0usize;
    let mut max_pos = 0usize;
    for i in 0..=(seq.len() - win) {
        let w = &seq[i..i + win];
        let g = w.iter().filter(|&&b| b == b'G').count() as f64;
        let c = w.iter().filter(|&&b| b == b'C').count() as f64;
        let skew = if g + c > 0.0 { (g - c) / (g + c) } else { 0.0 };
        cumulative += skew;
        if cumulative < min_val {
            min_val = cumulative;
            min_pos = i;
        }
        if cumulative > max_val {
            max_val = cumulative;
            max_pos = i;
        }
    }
    Some((min_pos, max_pos, min_val, max_val))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: &str) -> Vec<u8> {
        s.bytes().collect()
    }

    #[test]
    fn skew_returns_landmarks() {
        let seq = v("GGGGGGGGGGCCCCCCCCCC");
        let (origin, terminus, min_v, max_v) = gc_skew_landmarks(&seq, 5).unwrap();
        assert!(min_v <= max_v);
        assert!(origin <= seq.len() && terminus <= seq.len());
    }

    #[test]
    fn skew_clamps_oversized_window() {
        // A window larger than the sequence falls back to the whole sequence.
        assert!(gc_skew_landmarks(&v("ACGTGC"), 100).is_some());
    }

    #[test]
    fn skew_none_when_empty() {
        assert!(gc_skew_landmarks(&v(""), 10).is_none());
    }
}
