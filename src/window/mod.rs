//! Sliding-window composition scans.
//!
//! These measure base composition over a window rather than matching a pattern:
//! GC-rich regions and CpG islands (`gc`), and GC skew for locating replication
//! origins (`skew`).

mod gc;
mod skew;

pub(crate) use gc::find_gc_regions;
pub(crate) use skew::gc_skew_landmarks;
