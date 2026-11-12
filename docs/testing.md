# Tests

What each suite pins down. [← back to the README](../README.md)

`node tests/test_hgvs.mjs` tests the variant resolver. It extracts the code from
`site/index.html` between the `HGVS-RESOLVER` markers, so it exercises what is actually served,
against synthetic genes on both strands: CDS mapping, intron offsets in the right direction per
strand, reference-base read-back, and that UTR, out-of-range, wrong-chromosome and deep-intronic
lookups are refused rather than answered wrongly.

`node tests/test_page.mjs` tests the page around it: the whole `<script>` runs against a stubbed
DOM and a stubbed `fetch` serving synthetic genes **and a stubbed UCSC API**. It covers which gene
a variant resolves to — by symbol, Ensembl transcript, RefSeq, or coordinate — the
transcript-safety rules (a version difference warns, a different accession is offered for
conversion and never renumbered), the UCSC route end to end (0-based half-open rows become
1-based exons, a foreign `c.` lands on the right genomic base, a shared exon is recognised, a
gene outside the panel is reported as such), ΔTm display and re-ordering including the guard
that never promotes a pair with an off-target product, that a gene symbol is never routed as
a variant, and the oligo order rows: seven tab-separated columns in the order of the synthesis
order sheet, `GENE-E<exon>-F`/`-R` names, the sequence length, an empty PCR-positive-control
cell and the three fixed parameters.

`python tests/test_offline.py` builds a synthetic genome with genes on both strands and fake
Ensembl/UCSC/isPcr responses. The fakes serve the UCSC `mane` track and the batched
`POST /sequence/region` too, so the fast paths are the ones under test: that a model built from
the MANE track is **identical** to one built from Ensembl (exons, transcript, CDS bounds, region),
that a gene missing from the track falls back to Ensembl rather than being skipped, that one POST
covers many genes and the sequences match the single-region ones, and that a dead batch endpoint
falls back instead of losing data. It also pins both published Tm values and that exactly two are
published, the ΔTm ceiling, and that the salt conditions actually reach Primer3 — including a
check of `tm_ucsc` against the figures UCSC in-silico PCR itself prints for a known pair, and
one that the magnesium is really zeroed rather than merely the tables switched, so a primer3-py
upgrade or a tidy-up that moves either fails loudly instead of silently shifting a column. Beyond that it checks
primer coordinates against the genome, minimum distance, coverage of every coding base, SNP
avoidance and the 3′-end rule, the MAF cut-off, repeat avoidance and the 5′-only fallback,
segmental-duplication alerts, and that an off-target in-silico PCR product demotes a pair.

