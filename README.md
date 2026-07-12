# argus-ligand (finds features/motifs in DNA)
(i know is ligand used in checmical compounds but using rdkit on cloudflare workers :sigh: )
 
You can paste a sequence and either build a scan by hand or ask for what you want in plain English, and Workers AI turns that request into a scan.
> One Cloudflare Worker, One Workers AI (Llama 3.1 8B, JSON mode), One KV Cache, written in Rust (via [`workers-rs`](https://github.com/cloudflare/workers-rs)). The same Rust program serves the web page and does the DNA scanning.

> This is educational only. It is not a medical device, it works on example
> teaching sequences, and nothing here should inform any real health decision.
> Real clinical testing compares a person's DNA against a reference genome and
> databases of known variants, which is a different kind of pipeline.


</br>
</br>

## Example: seeing sickle cell anemia

Sickle cell anemia is a disease which is identified with only one letter. 
In the beta-globin gene the codon `GAG` (glutamate) becomes `GTG`
(valine): a single A to T change. 

Two of the scans here show that change on an example sequence.

Healthy allele:

```
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCC
```

Sickle allele (the single base flipped, `GAG` to `GTG`):

```
ATGGTGCACCTGACTCCTGTGGAGAAGTCTGCC
```

- **Way 1, motif search :** 
Search the motif `CCTGAGG`. It is present in the healthy allele and absent in the sickle allele. Searching `CCTGTGG` does the opposite. So one motif tells the two apart.

- **Way 2, restriction site (how it is really diagnosed) :** 
The healthy sequence `CCTGAGG` is a cut site for the enzyme MstII, whose recognition sequence is
`CCTNAGG`. The sickle mutation turns it into `CCTGTGG`, which MstII no longer
cuts. So a restriction scan for MstII finds one site in the healthy allele and
none in the sickle allele. "The enzyme stops cutting" is a genuine historical
diagnostic, called an RFLP (restriction fragment length polymorphism).

> Other conditions map onto the same ideas: many traits and disorders are single
> base changes a motif can spot, while repeat expansion disorders (for example
> Huntington's disease, caused by too many `CAG` repeats) are about counting a
> repeated unit.

</br>
</br>


## A tini-tiny biology (just asked AI to explained in simpler terms)

DNA is a long string over four letters called bases: A, C, G, and T. It is double stranded. Each base pairs with a complement (A with T, C with G), so the second strand is the reverse complement of the first. 

Soooo, a feature can sit on either strand, which is why most scans can look at both.

> A **motif** is a short recurring pattern in the sequence that carries meaning, for example a spot where a protein binds. 

> Motifs are written in **IUPAC codes**, which add letters that stand for "any of these bases". 

> For example W means A or T, R means A or G, and N means any base. So the pattern `TATAWAWR` describes a whole family of related sequences rather than one fixed string. That is why a real promoter search uses `TATAWAWR` and not just `TATAAA`: the plain string would miss the common variants `TATATAA`, `TATAAAA`, `TATAAAG`, and so on.

The scanner supports these analyses:

- **Motif (literal or IUPAC).** Where a pattern occurs, on either strand. The
  classic example is the TATA box, a signal found near the start of many genes.
- **ORF (open reading frame).** A stretch from a start codon (ATG) to a stop
  codon in the same reading frame. ORFs are candidate genes; translating one
  gives the protein it would code for.
- **GC content and CpG islands.** Regions unusually rich in G and C, measured
  with a sliding window. CpG islands (stretches with many CG dinucleotides) often
  sit at gene promoters.
- **GC skew.** The running imbalance between G and C along the sequence. Where
  the cumulative skew turns around, it marks the likely replication origin in
  many bacterial genomes.
- **Restriction sites.** Short specific sequences that a restriction enzyme cuts,
  for example EcoRI cuts `GAATTC`. Mapping these is a staple of cloning.
- **PWM (position weight matrix).** A scored, probabilistic motif built from
  example sites. It suits fuzzy patterns like transcription factor binding sites,
  where no single fixed string fits every real occurrence.

A site that reads the same on both strands (a palindrome, like the EcoRI site
`GAATTC`) is reported once with strand `±` instead of being double counted.

</br>
</br>


## Architecture
<p align="center">
<img width="711" height="701" alt="image" src="https://github.com/user-attachments/assets/fdcd5846-e80d-4080-91ad-b9e7e23bd4c4" />
</p>

</br>
</br>

## Project layout

```
📁 argus-ligand/
│
├── 📁 src/
│   ├── lib.rs            crate root: module wiring and re-exports
│   ├── error.rs          EngineError type
│   ├── types.rs          request and response JSON types (serde)
│   ├── sequence.rs       cleaning, reverse complement, Strand enum
│   ├── orf.rs            ORF finding and protein translation
│   ├── engine.rs         the scan orchestrator
│   ├── ai.rs             plain-English decomposition via Workers AI
│   ├── samples.rs        bundled real genomes, served from /samples
│   ├── worker_app.rs     HTTP handler + rate limit (wasm32 only)
│   │
│   ├── 📁 pattern/       pattern-matching scans
│   │   ├── mod.rs
│   │   ├── iupac.rs      IUPAC code to bitmask
│   │   ├── motif.rs      motif search over compiled masks
│   │   ├── restriction.rs  enzyme table + site mapping
│   │   └── pwm.rs        position weight matrix scoring
│   │
│   └── 📁 window/        sliding-window composition scans
│       ├── mod.rs
│       ├── gc.rs         GC-rich regions and CpG islands
│       └── skew.rs       GC skew (replication origin)
│
└── 📁 public/
    ├── index.html        UI markup (bento layout)
    ├── index.css         UI styles (terminal-brutalist theme)
    ├── index.js          UI logic
    │
    └── 📁 samples/
        ├── puc19.fasta   pUC19 cloning vector (2,686 bp)
        ├── pbr322.fasta  pBR322 plasmid (4,361 bp)
        └── lambda.fasta  phage lambda genome (48,502 bp)
```

</br>
</br>

## Where does the data come from?

A sequence enters three ways:

1. **You paste it** into the box (raw nucleotides or FASTA; `>` header and
   whitespace lines are stripped automatically).
2. **A built-in real genome** from the picker. Three genuine records are bundled
   from NCBI GenBank and embedded in the binary: pUC19, pBR322, and phage lambda.
3. **A tiny synthetic demo** (the TATA-box and ORF buttons) for showing a
   specific pattern quickly.

</br>
</br>

## Develop

```bash
# One-time toolchain
rustup target add wasm32-unknown-unknown
cargo install worker-build

# Test the pure Rust engine (native, no WASM, instant)
cargo test

# Run locally
npx wrangler dev

# Deploy to the internet (asks you to log in to Cloudflare)
npx wrangler deploy
```

</br>

---
<sub>Built with 💻 & ☕️ | © 2025 [Dhruv Jain](https://github.com/dhruvkjain)</sub>
