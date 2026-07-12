//! The scan orchestrator.
//!
//! This module ties the other modules together. It parses a request, cleans the
//! sequence once, computes the reverse complement once, then runs each scan and
//! collects the results.

use crate::error::EngineError;
use crate::orf::find_orfs;
use crate::pattern::{
    build_pssm, compile_iupac, find_motif_stranded, find_restriction_sites, scan_pwm,
};
use crate::sequence::{clean_sequence, reverse_complement, Strand};
use crate::types::{ScanRequest, ScanResponse, ScanResult, ScanSpec};
use crate::window::{find_gc_regions, gc_skew_landmarks};

/// Run every scan in a JSON request and return the results as JSON.
///
/// This is the single entry point to the engine. The input is a JSON
/// `ScanRequest` and the output is a JSON `ScanResponse`. The function never
/// panics: if the input cannot be handled it returns a JSON object with an
/// `error` field instead.
///
/// # Arguments
///
/// * `input` - A JSON string with a `sequence` and a list of `scans`. See the
///   crate root for the request shape.
///
/// # Returns
///
/// A JSON string. On success it is a `ScanResponse`. On failure it is
/// `{"error": "..."}` describing what went wrong.
///
/// # Example
///
/// ```rust
/// use argus_ligand::scan;
///
/// let request = r#"{"sequence":"GAATTC","scans":[
///     {"type":"motif","pattern":"GAATTC","both_strands":false}]}"#;
/// let out = scan(request);
/// assert!(out.contains("\"count\":1"));
/// ```
pub fn scan(input: &str) -> String {
    match run(input) {
        Ok(resp) => serde_json::to_string(&resp)
            .unwrap_or_else(|e| format!("{{\"error\":\"serialize: {e}\"}}")),
        Err(e) => serde_json::json!({ "error": e.to_string() }).to_string(),
    }
}

/// Parse and run a request into a typed response.
///
/// Cleaning and reverse complement happen once and are shared across all scans.
/// A bad motif pattern does not fail the whole request: that scan becomes a
/// [`ScanResult::Error`] entry and the other scans still run.
///
/// # Arguments
///
/// * `input` - The raw JSON request string.
///
/// # Returns
///
/// A [`ScanResponse`] holding one result per scan.
///
/// # Errors
///
/// Returns [`EngineError::InvalidRequest`] if the JSON does not parse or match
/// the expected shape. Returns [`EngineError::EmptySequence`] if the sequence is
/// empty after cleaning.
fn run(input: &str) -> Result<ScanResponse, EngineError> {
    let req: ScanRequest =
        serde_json::from_str(input).map_err(|e| EngineError::InvalidRequest(e.to_string()))?;
    run_request(req)
}

/// Run an already parsed request.
///
/// This is the shared core used by both the JSON string entry point [`scan`] and
/// the AI path, which builds a [`ScanRequest`] from typed scans rather than from
/// JSON text. Cleaning and reverse complement happen once and are shared across
/// all scans.
///
/// # Arguments
///
/// * `req` - The parsed request holding a sequence and a list of scans.
///
/// # Returns
///
/// A [`ScanResponse`] holding one result per scan.
///
/// # Errors
///
/// Returns [`EngineError::EmptySequence`] if the sequence is empty after cleaning.
pub(crate) fn run_request(req: ScanRequest) -> Result<ScanResponse, EngineError> {
    let (seq, mut warnings) = clean_sequence(&req.sequence);
    if seq.is_empty() {
        return Err(EngineError::EmptySequence);
    }

    let rc = reverse_complement(&seq);
    let mut results = Vec::with_capacity(req.scans.len());

    for spec in req.scans {
        results.push(run_one(spec, &seq, &rc));
    }

    warnings.shrink_to_fit();
    Ok(ScanResponse {
        seq_length: seq.len(),
        results,
        warnings,
    })
}

/// Run a single scan against the forward sequence and its reverse complement.
///
/// # Arguments
///
/// * `spec` - The scan to run.
/// * `seq` - The cleaned forward sequence.
/// * `rc` - The reverse complement of `seq`, precomputed by the caller.
///
/// # Returns
///
/// The [`ScanResult`] for this scan. An invalid motif pattern yields a
/// [`ScanResult::Error`] rather than failing the request.
fn run_one(spec: ScanSpec, seq: &[u8], rc: &[u8]) -> ScanResult {
    match spec {
        ScanSpec::Motif {
            pattern,
            both_strands,
        } => match compile_iupac(&pattern) {
            Ok(masks) => {
                let matches = find_motif_stranded(seq, rc, &masks, both_strands);
                ScanResult::Motif {
                    count: matches.len(),
                    pattern,
                    matches,
                }
            }
            Err(e) => ScanResult::Error {
                message: format!("motif '{pattern}': {e}"),
            },
        },
        ScanSpec::Orf {
            min_aa,
            both_strands,
        } => {
            let mut matches = find_orfs(seq, Strand::Forward, min_aa, seq.len());
            if both_strands {
                matches.extend(find_orfs(rc, Strand::Reverse, min_aa, seq.len()));
            }
            ScanResult::Orf {
                count: matches.len(),
                matches,
            }
        }
        ScanSpec::Gc {
            min_len,
            min_gc,
            min_cpg_ratio,
        } => {
            let matches = find_gc_regions(seq, min_len, min_gc, min_cpg_ratio);
            ScanResult::Gc {
                count: matches.len(),
                matches,
            }
        }
        ScanSpec::GcSkew { window } => match gc_skew_landmarks(seq, window) {
            Some((origin, terminus, min_skew, max_skew)) => ScanResult::GcSkew {
                window: window.clamp(1, seq.len()),
                origin,
                terminus,
                min_skew,
                max_skew,
            },
            None => ScanResult::Error {
                message: "sequence is shorter than the GC skew window".to_string(),
            },
        },
        ScanSpec::Restriction {
            enzymes,
            both_strands,
        } => {
            let matches = find_restriction_sites(seq, rc, &enzymes, both_strands);
            ScanResult::Restriction {
                count: matches.len(),
                matches,
            }
        }
        ScanSpec::Pwm {
            sites,
            threshold,
            both_strands,
        } => match build_pssm(&sites) {
            Ok(pssm) => {
                let matches = scan_pwm(seq, rc, &pssm, threshold, both_strands);
                ScanResult::Pwm {
                    count: matches.len(),
                    max_score: pssm.max_score,
                    matches,
                }
            }
            Err(e) => ScanResult::Error {
                message: format!("pwm: {e}"),
            },
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn end_to_end_json() {
        let req = r#"{"sequence":">demo\nATGAAAAAATAAGAATTC","scans":[
            {"type":"motif","pattern":"GAATTC","both_strands":false},
            {"type":"orf","min_aa":1,"both_strands":false}]}"#;
        let out = scan(req);
        assert!(out.contains("\"seq_length\":18"));
        assert!(out.contains("MKK"));
        assert!(out.contains("\"count\":1"));
    }

    #[test]
    fn invalid_json_returns_error_object() {
        let out = scan("not json");
        assert!(out.contains("\"error\""));
        assert!(out.contains("invalid request"));
    }

    #[test]
    fn empty_sequence_returns_error_object() {
        let out = scan(r#"{"sequence":">only header","scans":[]}"#);
        assert!(out.contains("sequence is empty"));
    }

    #[test]
    fn invalid_pattern_does_not_fail_whole_request() {
        let req = r#"{"sequence":"ACGT","scans":[
            {"type":"motif","pattern":"TAXQ","both_strands":false}]}"#;
        let out = scan(req);
        assert!(out.contains("not a valid IUPAC code"));
    }
}
