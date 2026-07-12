//! argus-ligand: a DNA motif scanner as a 100% Rust Cloudflare Worker.
//!
//! The crate has two layers. The engine layer is pure Rust and does the DNA
//! work: cleaning input, matching motifs, and finding open reading frames. The
//! Worker layer ([`worker_app`], compiled only for `wasm32`) serves the UI and
//! exposes the engine over HTTP. The single engine entry point is [`scan`].
//!
//! Because the `worker` dependency is limited to the `wasm32` target, running
//! `cargo test` on a normal machine compiles only the engine layer, which keeps
//! the tests fast and free of the Worker runtime.
//!
//! # Request shape
//!
//! `POST /scan` takes a JSON object like this:
//!
//! ```json
//! {
//!   "sequence": ">optional_fasta_header\nACGTACGT...",
//!   "scans": [
//!     { "type": "motif", "pattern": "TATAWAWR", "both_strands": true },
//!     { "type": "orf", "min_aa": 30, "both_strands": true }
//!   ]
//! }
//! ```
//!
//! `sequence` may be raw nucleotides or FASTA. Motif patterns use IUPAC codes.
//! Coordinates in the response are 0 based and half open, `[start, end)`.
//!
//! # Example
//!
//! ```rust
//! use argus_ligand::scan;
//!
//! let request = r#"{"sequence":"ATGAAATAA","scans":[
//!     {"type":"orf","min_aa":1,"both_strands":false}]}"#;
//! let out = scan(request);
//! assert!(out.contains("\"protein\":\"MK\""));
//! ```

mod ai;
mod engine;
mod error;
mod orf;
mod pattern;
mod samples;
mod sequence;
mod types;
mod window;

#[cfg(target_arch = "wasm32")]
mod worker_app;

pub use engine::scan;
pub use error::EngineError;
