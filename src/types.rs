//! Request and response types for the `/scan` API.
//!
//! These map one to one to the JSON shapes documented in the crate root. Serde
//! handles the conversion, so the engine works with typed values rather than
//! raw JSON.

use serde::{Deserialize, Serialize};

/// A full scan request: one sequence and a list of scans to run on it.
#[derive(Deserialize)]
pub(crate) struct ScanRequest {
    /// Raw nucleotides or FASTA text. Headers and whitespace are stripped later.
    pub(crate) sequence: String,
    /// The scans to run, in order.
    pub(crate) scans: Vec<ScanSpec>,
}

/// One scan to run, tagged by its `type` field in JSON.
///
/// Deserialized from a `/scan` request, and also serialized back in an `/ask`
/// response so the caller can see how the AI interpreted their prompt.
#[derive(Deserialize, Serialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub(crate) enum ScanSpec {
    /// A literal or IUPAC motif search.
    Motif {
        /// The pattern in IUPAC codes, for example `TATAWAWR`.
        pattern: String,
        /// Whether to also search the reverse complement strand.
        #[serde(default = "default_true")]
        both_strands: bool,
    },
    /// An open reading frame search.
    Orf {
        /// Minimum protein length in amino acids, excluding the stop codon.
        #[serde(default = "default_min_aa")]
        min_aa: usize,
        /// Whether to also search the reverse complement strand.
        #[serde(default = "default_true")]
        both_strands: bool,
    },
    /// A GC-rich region or CpG-island search. Strand independent.
    Gc {
        /// Minimum region length in nucleotides.
        #[serde(default = "default_min_len")]
        min_len: usize,
        /// Minimum GC fraction, from 0.0 to 1.0.
        #[serde(default = "default_min_gc")]
        min_gc: f64,
        /// Minimum observed/expected CpG ratio. 0 accepts any GC-rich region;
        /// about 0.6 is the classic CpG-island threshold.
        #[serde(default = "default_min_cpg_ratio")]
        min_cpg_ratio: f64,
    },
    /// GC skew across a sliding window, used to locate replication origins.
    #[serde(rename = "gc_skew")]
    GcSkew {
        /// Sliding window size in nucleotides.
        #[serde(default = "default_skew_window")]
        window: usize,
    },
    /// Restriction enzyme recognition-site mapping.
    Restriction {
        /// Enzyme names to search. Empty means all known enzymes.
        #[serde(default)]
        enzymes: Vec<String>,
        /// Whether to also search the reverse complement strand.
        #[serde(default = "default_true")]
        both_strands: bool,
    },
    /// A position weight matrix scan, built from example sites.
    Pwm {
        /// Equal-length example binding sites the matrix is built from.
        sites: Vec<String>,
        /// Score threshold as a fraction (0 to 1) of the maximum possible score.
        #[serde(default = "default_pwm_threshold")]
        threshold: f64,
        /// Whether to also search the reverse complement strand.
        #[serde(default = "default_true")]
        both_strands: bool,
    },
}

/// Default for `both_strands`: search both strands unless told otherwise.
fn default_true() -> bool {
    true
}

/// Default minimum protein length for an ORF scan.
fn default_min_aa() -> usize {
    30
}

/// Default minimum region length for a GC scan (the classic CpG-island size).
fn default_min_len() -> usize {
    200
}

/// Default minimum GC fraction for a GC scan.
fn default_min_gc() -> f64 {
    0.5
}

/// Default minimum CpG ratio for a GC scan. 0 means "any GC-rich region".
fn default_min_cpg_ratio() -> f64 {
    0.0
}

/// Default sliding window size for a GC skew scan.
fn default_skew_window() -> usize {
    100
}

/// Default PWM score threshold as a fraction of the maximum possible score.
fn default_pwm_threshold() -> f64 {
    0.8
}

/// The full response: the cleaned sequence length and one result per scan.
#[derive(Serialize)]
pub(crate) struct ScanResponse {
    /// Length of the sequence after cleaning, in nucleotides.
    pub(crate) seq_length: usize,
    /// One entry per scan in the request, in the same order.
    pub(crate) results: Vec<ScanResult>,
    /// Non fatal notes, for example dropped characters. Omitted from JSON when empty.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub(crate) warnings: Vec<String>,
}

/// The result of one scan, tagged by `type` in JSON.
#[derive(Serialize)]
#[serde(tag = "type", rename_all = "lowercase")]
pub(crate) enum ScanResult {
    /// Motif search results.
    Motif {
        /// The pattern that was searched.
        pattern: String,
        /// Number of matches found.
        count: usize,
        /// The matches, in ascending forward strand position.
        matches: Vec<MotifMatch>,
    },
    /// ORF search results.
    Orf {
        /// Number of ORFs found.
        count: usize,
        /// The ORFs found.
        matches: Vec<OrfMatch>,
    },
    /// GC-rich region / CpG-island results.
    Gc {
        /// Number of regions found.
        count: usize,
        /// The regions found.
        matches: Vec<GcMatch>,
    },
    /// GC skew landmark results.
    #[serde(rename = "gc_skew")]
    GcSkew {
        /// The window size used.
        window: usize,
        /// Position of minimum cumulative skew (putative replication origin).
        origin: usize,
        /// Position of maximum cumulative skew (putative terminus).
        terminus: usize,
        /// The minimum cumulative skew value.
        min_skew: f64,
        /// The maximum cumulative skew value.
        max_skew: f64,
    },
    /// Restriction site results.
    Restriction {
        /// Number of sites found.
        count: usize,
        /// The sites found.
        matches: Vec<RestrictionMatch>,
    },
    /// Position weight matrix results.
    Pwm {
        /// Number of hits found.
        count: usize,
        /// The maximum possible score of the built matrix, for context.
        max_score: f64,
        /// The hits found.
        matches: Vec<PwmMatch>,
    },
    /// A per scan error, for example an invalid pattern. Other scans still run.
    Error {
        /// What went wrong with this scan.
        message: String,
    },
}

/// One motif match, in forward strand coordinates.
#[derive(Serialize)]
pub(crate) struct MotifMatch {
    /// Start position, 0 based, inclusive.
    pub(crate) start: usize,
    /// End position, 0 based, exclusive.
    pub(crate) end: usize,
    /// Strand the match was found on: `'+'` or `'-'`.
    pub(crate) strand: char,
    /// The matched substring as it reads on `strand`.
    pub(crate) matched: String,
}

/// One GC-rich region or CpG island.
#[derive(Serialize)]
pub(crate) struct GcMatch {
    /// Start position, 0 based, inclusive.
    pub(crate) start: usize,
    /// End position, 0 based, exclusive.
    pub(crate) end: usize,
    /// Region length in nucleotides.
    pub(crate) length: usize,
    /// GC content of the region as a percentage, 0 to 100.
    pub(crate) gc_percent: f64,
    /// Observed/expected CpG ratio of the region.
    pub(crate) cpg_ratio: f64,
}

/// One restriction enzyme site.
#[derive(Serialize)]
pub(crate) struct RestrictionMatch {
    /// Start position, 0 based, inclusive.
    pub(crate) start: usize,
    /// End position, 0 based, exclusive.
    pub(crate) end: usize,
    /// Strand: `'+'`, `'-'`, or `'±'` for a palindromic site.
    pub(crate) strand: char,
    /// The enzyme name, for example `EcoRI`.
    pub(crate) enzyme: String,
    /// The recognition site, for example `GAATTC`.
    pub(crate) site: String,
}

/// One position weight matrix hit, in forward strand coordinates.
#[derive(Serialize)]
pub(crate) struct PwmMatch {
    /// Start position, 0 based, inclusive.
    pub(crate) start: usize,
    /// End position, 0 based, exclusive.
    pub(crate) end: usize,
    /// Strand the hit was found on: `'+'` or `'-'`.
    pub(crate) strand: char,
    /// The match score under the matrix.
    pub(crate) score: f64,
}

/// One open reading frame, in forward strand coordinates.
#[derive(Serialize)]
pub(crate) struct OrfMatch {
    /// Start position of the ORF including the start codon, 0 based, inclusive.
    pub(crate) start: usize,
    /// End position of the ORF including the stop codon, 0 based, exclusive.
    pub(crate) end: usize,
    /// Strand the ORF was found on: `'+'` or `'-'`.
    pub(crate) strand: char,
    /// Reading frame on that strand: 0, 1, or 2.
    pub(crate) frame: usize,
    /// Protein length in amino acids, excluding the stop.
    pub(crate) aa_length: usize,
    /// The translated protein, one letter per amino acid.
    pub(crate) protein: String,
}
