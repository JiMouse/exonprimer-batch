# Design method

How each amplicon is chosen, and what is checked before it is called OK.
[← back to the README](../README.md)

## What it does

For each gene in `genes.txt`:

1. Gets the exon structure of the **MANE Select** transcript (RefSeq NM reported) from the UCSC
   `mane` track — one request per chromosome for the whole genome. Ensembl REST is asked per gene
   only for what the track cannot answer: a transcript you name, a gene with no MANE transcript
   (Ensembl canonical is then used), or a row the track could not provide.
2. Takes each **coding exon** as a target (by default the coding part only; `--target full` for
   whole exons). Exons are numbered by their rank in the transcript.
3. Applies ExonPrimer's rules:
   - primer 3′ end at least `--min-dist` bp (default 50) from the exon/intron boundary,
     within a `--primer-region` window (default 150 bp) beyond that;
   - **PCR products are 200–500 bp** by default (`--product MIN MAX`);
   - exons separated by a short intron are **combined** into one amplicon if they still fit in
     one product — e.g. `E2-3` (turn off with `--no-combine`);
   - exons too long for one product get **several primer pairs** on overlapping fragments
     (`--overlap`, default 50 bp) — e.g. `E9a`, `E9b`, `E9c`. The largest target per pair is
     derived from the product range: 500 − 2 × (50 bp min distance + 20 nt primer + 20 bp margin)
     = 320 bp (override with `--max-target`);
   - repeats are soft-masked and primers with a 3′ end in a repeat are rejected
     (the equivalent of ExonPrimer's mispriming-library option).
4. Runs Primer3 (Tm 57–63 °C, opt 60; 18–27 nt, opt 20; GC 30–70 %; GC clamp), with the two
   primers of a pair matched to within **2 °C** (`--max-tm-diff`, and `--wt-tm-diff` for how hard
   Primer3 prefers a matched pair — see *Tm* below). If nothing is found it retries once with
   relaxed settings and flags the result.

### Safeguards (as in EasyExonPrimer, Wu & Munroe 2006)

**1. Common SNPs under primers.** Common variants from UCSC `dbSnp155Common` (dbSNP 155,
1000 Genomes MAF ≥ 1 %; raise the cut-off with `--snp-maf 0.05`) are masked before Primer3.

**2. Repeats.** Ensembl soft-masking (RepeatMasker, TRF, Dust) plus UCSC `rmsk` are masked.

Design runs in tiers and stops at the first that works:

| Tier | Accepted | Status |
|---|---|---|
| `clean` | no common SNP and no repeat anywhere under either primer | OK |
| `repeat-5p` | a repeat touches only the 5′ part of a primer (3′ end in unique sequence) | Check |
| `snp-5p` | a common SNP lies under a primer, but never in its last 5 bases (`--snp-3prime`) | Check |

`--strict` keeps only the `clean` tier (exons without a clean pair are reported as failed).
Every SNP left under a final primer is listed with rsID, MAF and distance from the 3′ end.

**3. High-homology sequences.**
- Targets overlapping a segmental duplication (UCSC `genomicSuperDups`, ≥ 90 % identity, ≥ 1 kb —
  pseudogenes and paralogs such as *PMS2*, *SMN1*, *CYP21A2*, *GBA1*) are flagged with the
  identity and the duplicated locus.
- With `--ispcr` and `--genome-2bit`, every candidate pair (best + extras) is run through
  UCSC **isPcr** genome-wide (same settings as the UCSC In-Silico PCR web tool). Specific pairs
  are preferred; any off-target product is listed. Hits on alt/fix haplotype contigs are ignored.
  The GitHub workflow downloads isPcr and the genome automatically.

An amplicon is **OK** only if it is `clean`, designed with normal settings, outside segmental
duplications and (when checked) specific. Otherwise it is **Check**, with the reasons listed.

Set the parameters to match the ones you use on the ExonPrimer form (`--product`, `--min-dist`,
`--primer-region`, `--max-target`, `--overlap`, `--tm`, `--size`, `--gc`, `--gc-clamp`).

