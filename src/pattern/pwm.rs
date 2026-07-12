//! Position weight matrix (PWM) scanning.
//!
//! A PWM scores how well each position of a window matches a set of example
//! sites, which suits fuzzy, probabilistic motifs like transcription factor
//! binding sites that no single fixed pattern captures. The matrix is built from
//! equal length example sites as a log-odds score against a uniform background,
//! then slid across the sequence. A window is reported when its total score
//! reaches a fraction of the best possible score.

use crate::error::EngineError;
use crate::types::PwmMatch;

/// A built position weight matrix: per position log-odds scores and the best
/// possible total score.
pub(crate) struct Pssm {
    /// One `[A, C, G, T]` score row per position.
    scores: Vec<[f64; 4]>,
    /// The maximum total score any window can reach.
    pub(crate) max_score: f64,
}

/// Map a base byte to its column index.
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Build a PWM from example sites.
///
/// All sites must be the same non-zero length and contain only A, C, G, T
/// (U is accepted as T). Frequencies use a Laplace pseudocount, and scores are
/// `log2(frequency / 0.25)` against a uniform background.
///
/// # Arguments
///
/// * `sites` - Equal length example binding sites.
///
/// # Returns
///
/// The built [`Pssm`].
///
/// # Errors
///
/// Returns [`EngineError::InvalidRequest`] if there are no sites, the sites are
/// not all the same non-zero length, or a site contains a non-nucleotide letter.
pub(crate) fn build_pssm(sites: &[String]) -> Result<Pssm, EngineError> {
    let cleaned: Vec<Vec<u8>> = sites
        .iter()
        .map(|s| {
            s.trim()
                .bytes()
                .filter(|b| !b.is_ascii_whitespace())
                .map(|b| {
                    let up = b.to_ascii_uppercase();
                    if up == b'U' {
                        b'T'
                    } else {
                        up
                    }
                })
                .collect()
        })
        .filter(|v: &Vec<u8>| !v.is_empty())
        .collect();

    if cleaned.is_empty() {
        return Err(EngineError::InvalidRequest(
            "PWM needs at least one example site".to_string(),
        ));
    }
    let len = cleaned[0].len();
    if cleaned.iter().any(|s| s.len() != len) {
        return Err(EngineError::InvalidRequest(
            "PWM example sites must all be the same length".to_string(),
        ));
    }

    let pseudo = 1.0f64;
    let denom = cleaned.len() as f64 + 4.0 * pseudo;
    let mut scores = Vec::with_capacity(len);
    let mut max_score = 0.0f64;
    for pos in 0..len {
        let mut counts = [0.0f64; 4];
        for site in &cleaned {
            match base_index(site[pos]) {
                Some(idx) => counts[idx] += 1.0,
                None => {
                    return Err(EngineError::InvalidRequest(format!(
                        "PWM site has a non-nucleotide letter '{}'",
                        site[pos] as char
                    )))
                }
            }
        }
        let mut row = [0.0f64; 4];
        let mut best = f64::NEG_INFINITY;
        for b in 0..4 {
            let freq = (counts[b] + pseudo) / denom;
            row[b] = (freq / 0.25).log2();
            best = best.max(row[b]);
        }
        scores.push(row);
        max_score += best;
    }

    Ok(Pssm { scores, max_score })
}

/// Score a window against the matrix.
///
/// Returns `None` if the window contains a non-ACGT base, which cannot be scored.
fn score_window(pssm: &Pssm, win: &[u8]) -> Option<f64> {
    let mut total = 0.0;
    for (pos, &b) in win.iter().enumerate() {
        total += pssm.scores[pos][base_index(b)?];
    }
    Some(total)
}

/// Scan a sequence with a PWM.
///
/// Reports every window whose score reaches `threshold` times the matrix's
/// maximum possible score. The reverse strand is searched by scanning the
/// reverse complement and mapping coordinates back to the forward strand.
///
/// # Arguments
///
/// * `seq` - The forward strand sequence.
/// * `rc` - The reverse complement of `seq`, precomputed by the caller.
/// * `pssm` - The built matrix.
/// * `threshold` - Fraction (0 to 1) of `pssm.max_score` a hit must reach.
/// * `both_strands` - Whether to also search the reverse complement.
///
/// # Returns
///
/// Hits in ascending forward strand position.
pub(crate) fn scan_pwm(
    seq: &[u8],
    rc: &[u8],
    pssm: &Pssm,
    threshold: f64,
    both_strands: bool,
) -> Vec<PwmMatch> {
    let len = pssm.scores.len();
    let cutoff = threshold * pssm.max_score;
    let mut hits = Vec::new();
    if len == 0 || seq.len() < len {
        return hits;
    }

    for i in 0..=(seq.len() - len) {
        if let Some(score) = score_window(pssm, &seq[i..i + len]) {
            if score >= cutoff {
                hits.push(PwmMatch {
                    start: i,
                    end: i + len,
                    strand: '+',
                    score,
                });
            }
        }
    }
    if both_strands {
        let fwd_len = seq.len();
        for i in 0..=(rc.len() - len) {
            if let Some(score) = score_window(pssm, &rc[i..i + len]) {
                if score >= cutoff {
                    hits.push(PwmMatch {
                        start: fwd_len - (i + len),
                        end: fwd_len - i,
                        strand: '-',
                        score,
                    });
                }
            }
        }
    }
    hits.sort_by(|a, b| a.start.cmp(&b.start));
    hits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::reverse_complement;

    fn sites(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn builds_and_scores_consensus() {
        let pssm = build_pssm(&sites(&["TATAAA", "TATAAA", "TATATA"])).unwrap();
        assert!(pssm.max_score > 0.0);
        let seq: Vec<u8> = "GGGTATAAACCC".bytes().collect();
        let rc = reverse_complement(&seq);
        let hits = scan_pwm(&seq, &rc, &pssm, 0.8, false);
        assert!(hits.iter().any(|h| h.start == 3 && h.strand == '+'));
    }

    #[test]
    fn rejects_unequal_lengths() {
        assert!(build_pssm(&sites(&["TATA", "TAT"])).is_err());
    }

    #[test]
    fn rejects_empty() {
        assert!(build_pssm(&[]).is_err());
    }

    #[test]
    fn rejects_bad_base() {
        assert!(build_pssm(&sites(&["TAXA", "TATA"])).is_err());
    }
}
