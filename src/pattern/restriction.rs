//! Restriction enzyme site mapping.
//!
//! Each enzyme has a recognition sequence in IUPAC codes. Finding its sites is
//! just a motif search, so this module reuses [`crate::pattern::motif::find_motif_stranded`]
//! and tags each hit with the enzyme name. Most recognition sites are
//! palindromic, so they come back with the `±` strand from the motif engine.

use crate::pattern::iupac::compile_iupac;
use crate::pattern::motif::find_motif_stranded;
use crate::types::RestrictionMatch;

/// Built-in table of common enzymes: name and recognition sequence.
pub(crate) const ENZYMES: &[(&str, &str)] = &[
    ("EcoRI", "GAATTC"),
    ("BamHI", "GGATCC"),
    ("HindIII", "AAGCTT"),
    ("NotI", "GCGGCCGC"),
    ("XhoI", "CTCGAG"),
    ("PstI", "CTGCAG"),
    ("SmaI", "CCCGGG"),
    ("KpnI", "GGTACC"),
    ("SacI", "GAGCTC"),
    ("SalI", "GTCGAC"),
    ("XbaI", "TCTAGA"),
    ("SphI", "GCATGC"),
    ("NcoI", "CCATGG"),
    ("NdeI", "CATATG"),
    ("EcoRV", "GATATC"),
    ("HaeIII", "GGCC"),
    ("AluI", "AGCT"),
    ("TaqI", "TCGA"),
    ("NheI", "GCTAGC"),
    ("BglII", "AGATCT"),
    ("ClaI", "ATCGAT"),
    ("MluI", "ACGCGT"),
    ("ApaI", "GGGCCC"),
    ("DraI", "TTTAAA"),
    ("HpaI", "GTTAAC"),
    // CCTNAGG. Its site spans beta-globin codons 5-6; the sickle cell mutation
    // (GAG -> GTG) destroys it, which is the classic RFLP diagnostic.
    ("MstII", "CCTNAGG"),
];

/// Map restriction sites in a sequence.
///
/// Searches each enzyme in the built-in table, or only the named subset when
/// `enzymes` is non-empty. Enzyme names are matched case insensitively; unknown
/// names are ignored.
///
/// # Arguments
///
/// * `seq` - The forward strand sequence.
/// * `rc` - The reverse complement of `seq`, precomputed by the caller.
/// * `enzymes` - Enzyme names to include. Empty means all known enzymes.
/// * `both_strands` - Whether to also search the reverse complement.
///
/// # Returns
///
/// All sites found, sorted by start position.
pub(crate) fn find_restriction_sites(
    seq: &[u8],
    rc: &[u8],
    enzymes: &[String],
    both_strands: bool,
) -> Vec<RestrictionMatch> {
    let mut hits = Vec::new();
    for &(name, site) in ENZYMES {
        if !enzymes.is_empty() && !enzymes.iter().any(|e| e.eq_ignore_ascii_case(name)) {
            continue;
        }
        if let Ok(masks) = compile_iupac(site) {
            for m in find_motif_stranded(seq, rc, &masks, both_strands) {
                hits.push(RestrictionMatch {
                    start: m.start,
                    end: m.end,
                    strand: m.strand,
                    enzyme: name.to_string(),
                    site: site.to_string(),
                });
            }
        }
    }
    hits.sort_by_key(|h| h.start);
    hits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::reverse_complement;

    fn v(s: &str) -> Vec<u8> {
        s.bytes().collect()
    }

    #[test]
    fn finds_ecori_site() {
        let seq = v("AAAGAATTCAAA");
        let rc = reverse_complement(&seq);
        let hits = find_restriction_sites(&seq, &rc, &[], true);
        assert!(hits.iter().any(|h| h.enzyme == "EcoRI" && h.start == 3));
    }

    #[test]
    fn filters_by_enzyme_name() {
        let seq = v("AAAGAATTCAAAGGATCCAAA"); // EcoRI + BamHI
        let rc = reverse_complement(&seq);
        let only_bam = find_restriction_sites(&seq, &rc, &["BamHI".to_string()], true);
        assert!(only_bam.iter().all(|h| h.enzyme == "BamHI"));
        assert_eq!(only_bam.len(), 1);
    }

    #[test]
    fn unknown_enzyme_yields_nothing() {
        let seq = v("AAAGAATTCAAA");
        let rc = reverse_complement(&seq);
        assert!(find_restriction_sites(&seq, &rc, &["Nonsense".to_string()], true).is_empty());
    }

    #[test]
    fn mstii_distinguishes_sickle_cell() {
        // Beta-globin exon 1. The healthy allele carries an MstII site
        // (CCTGAGG); the sickle mutation GAG -> GTG turns it into CCTGTGG,
        // which MstII no longer cuts. This is the classic RFLP diagnostic.
        let healthy = v("ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCC");
        let sickle = v("ATGGTGCACCTGACTCCTGTGGAGAAGTCTGCC");
        let only_mstii = ["MstII".to_string()];
        let h = find_restriction_sites(&healthy, &reverse_complement(&healthy), &only_mstii, true);
        let s = find_restriction_sites(&sickle, &reverse_complement(&sickle), &only_mstii, true);
        assert_eq!(h.len(), 1); // cuts the healthy allele
        assert_eq!(s.len(), 0); // sickle mutation destroys the site
    }
}
