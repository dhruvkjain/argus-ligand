//! Pattern-matching scans.
//!
//! Everything here finds where a sequence matches a pattern: exact and IUPAC
//! motifs, restriction enzyme sites (a motif search over a named table), and
//! position weight matrices for scored, probabilistic motifs. The `iupac`
//! submodule holds the shared code -> bitmask primitives.

pub(crate) mod iupac;
pub(crate) mod motif;
pub(crate) mod pwm;
pub(crate) mod restriction;

pub(crate) use iupac::{base_mask, compile_iupac};
pub(crate) use motif::find_motif_stranded;
pub(crate) use pwm::{build_pssm, scan_pwm};
pub(crate) use restriction::find_restriction_sites;
