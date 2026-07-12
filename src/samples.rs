//! Bundled real sample sequences.
//!
//! These are genuine records downloaded from NCBI GenBank, embedded into the
//! binary at compile time. They give the UI real data to scan instead of only
//! the tiny synthetic demos. The engine strips the FASTA header on its own, so
//! the raw file content is served straight to the sequence box.

// These are consumed only by the wasm-only Worker routes, so they look dead on a
// host build. Allow that off-wasm only; the wasm target still checks normally.
#![cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]

use serde::Serialize;

/// One bundled sample: a stable id, a human label, and its FASTA text.
pub(crate) struct Sample {
    /// Stable id used in the `/samples/{id}` route.
    pub(crate) id: &'static str,
    /// Human readable label shown in the picker.
    pub(crate) label: &'static str,
    /// The raw FASTA text, header included.
    pub(crate) fasta: &'static str,
}

/// All bundled samples, in display order.
pub(crate) const SAMPLES: &[Sample] = &[
    Sample {
        id: "puc19",
        label: "pUC19 cloning vector (2,686 bp)",
        fasta: include_str!("../public/samples/puc19.fasta"),
    },
    Sample {
        id: "pbr322",
        label: "pBR322 plasmid (4,361 bp)",
        fasta: include_str!("../public/samples/pbr322.fasta"),
    },
    Sample {
        id: "lambda",
        label: "Phage lambda genome (48,502 bp)",
        fasta: include_str!("../public/samples/lambda.fasta"),
    },
];

/// A sample entry without its sequence, for the picker list.
#[derive(Serialize)]
pub(crate) struct SampleInfo {
    /// Stable id used in the `/samples/{id}` route.
    pub(crate) id: &'static str,
    /// Human readable label.
    pub(crate) label: &'static str,
}

/// List the bundled samples as id and label pairs.
///
/// # Returns
///
/// A vector of [`SampleInfo`], one per sample, in display order.
pub(crate) fn list() -> Vec<SampleInfo> {
    SAMPLES
        .iter()
        .map(|s| SampleInfo {
            id: s.id,
            label: s.label,
        })
        .collect()
}

/// Look up a sample's FASTA text by id.
///
/// # Arguments
///
/// * `id` - The sample id from the picker.
///
/// # Returns
///
/// `Some(fasta)` if a sample with that id exists, `None` otherwise.
pub(crate) fn fasta(id: &str) -> Option<&'static str> {
    SAMPLES.iter().find(|s| s.id == id).map(|s| s.fasta)
}
