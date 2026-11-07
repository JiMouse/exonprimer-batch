# Variant search

What the search box accepts, and the three things it refuses to guess.
[← back to the README](../README.md)

## Find the primers for a variant

The search box takes a variant as well as a gene symbol, and answers with the amplicon that
sequences it:

```
NM_000110.4:c.1129G>A      POLG:c.2243G>C       POLG:p.Trp748Ser
c.2243-12T>C               CFTR:c.1521_1523del  chr15:89327201
```

It resolves `c.`, `n.`, `g.` and `p.` descriptions, intron offsets (`c.1129+5`), ranges
(`c.1129_1130del`) and plain coordinates. A gene symbol always wins over variant parsing, so a
symbol is never misread. Links are shareable: `…/#v=POLG:c.2243G%3EC`.

This runs entirely in the browser — the published JSON already carries each transcript's CDS
structure, strand and amplicon spans, so no server or API is involved.

Three things it deliberately does **not** guess:

- **A different transcript.** The panel is built on MANE Select, so `c.` coordinates from
  another transcript are never silently renumbered against ours. A matching accession with a
  different version is used but flagged; a *different* accession is offered for conversion
  instead (below).
- **Deep intronic variants.** Primers reach about `--primer-region` bp into the flanks, so a
  variant further in is reported as covered by no amplicon, never rounded to the nearest exon.
- **UTR variants.** `c.-45` and `c.*120` lie outside the coding target and are refused.

When the variant states a reference allele it is checked against the amplicon sequence, which
catches a wrong genome build or a typo. Primers are matched on the **PCR product** span, not the
exon, so an intronic variant inside the product is correctly reported as sequenceable.

### A transcript the panel does not use

`NM_001316329.2:c.500A>G` is numbered on a *PRKCG* transcript that is not MANE Select, so the
published data cannot map it. The page says so and offers one button. On click it fetches that
accession's exon structure from the **UCSC API** (`api.genome.ucsc.edu`, which sends
`Access-Control-Allow-Origin: *`), converts the position in the browser with the same arithmetic
used for our own transcripts, and shows the amplicon. Still no server of ours, and nothing is
requested until you click.

Because an amplicon is keyed on the **genomic position**, the transcript question disappears once
the position is known — and alternative transcripts very often share the exon outright. For
*PRKCG*, `NM_001316329.2` and the MANE `NM_002739.5` have **17 of 19 exons at identical
coordinates**, so the pair to order is literally the same one. The page says which transcript did
the numbering, whether the exon is shared, and warns when it is not. If UCSC reports a gene the
panel does not cover (`NM_006088.6` → *TUBB4B*), you get that message instead of a
mis-mapped exon.

A gene or accession that is simply absent says so — *check the spelling, otherwise it is not in
our database* — with the closest symbols in the panel as links, rather than claiming the input
was unreadable HGVS.

The ΔTm checkbox above the search box is described in
[Tm](tm.md#on-the-site-alternatives-re-ordered-by-δtm).
