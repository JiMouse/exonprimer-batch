# exonprimer-batch

### → **[jimouse.github.io/exonprimer-batch](https://jimouse.github.io/exonprimer-batch)**

Intronic PCR/Sanger primers for **every coding exon of a list of genes**, designed the
ExonPrimer way, exported to **Excel**, and published as a **GitHub Pages site** where you pick a
gene and an exon — or paste a variant — to get the precomputed primers. Search a gene (`MTHFR`),
a variant (`POLG:c.2243G>C`) or a position (`chr15:89327201`), and download the whole panel as
one Excel workbook.

The gene list covers **6,437 genes** (GRCh38, MANE Select) — the original panel merged with the
**OMIM-morbid set**, every gene with an OMIM disease. **6,370** of them reach the site; the 67
that do not are listed in the Problems sheet — 59 have no protein-coding transcript at all
(ncRNAs and pseudogenes such as *ATXN8OS*, *CYP21A1P*, *DUX4L1*, which have no coding exon to
target) and 8 are symbols Ensembl cannot resolve.

> ExonPrimer itself (IHG, Helmholtz Munich) is a web form that does not allow automated
> access, so this project reproduces its method locally with Primer3 instead of querying it.

## Documentation

| | |
|---|---|
| **[Design method](docs/method.md)** | ExonPrimer's rules, combining and splitting exons, and the SNP / repeat / segmental-duplication / in-silico-PCR safeguards behind every OK |
| **[Tm](docs/tm.md)** | which of the two Tm columns to anneal to, why UCSC prints a different number, and how the two primers of a pair are matched |
| **[Variant search](docs/variant-search.md)** | what the search box accepts, and the three things it refuses to guess |
| **[Building and publishing](docs/build-and-deploy.md)** | local options, what a full build costs, the sharded GitHub workflow and the data cache |
| **[Tests](docs/testing.md)** | what each suite pins down |

## In short

For each gene in `genes.txt`:

- exon structure from the **MANE Select** transcript (UCSC `mane` track, Ensembl only for what
  the track cannot answer);
- each **coding exon** is a target, with primer 3′ ends at least 50 bp into the flanking intron;
- products **200–500 bp** — exons separated by a short intron are merged into one amplicon
  (`E2-3`), exons too long for one product are split into overlapping fragments (`E9a`, `E9b`);
- Primer3 at Tm 57–63 °C (opt 60), 18–27 nt, GC 30–70 %, the two primers matched to within 2 °C;
- common SNPs (UCSC `dbSnp155Common`, MAF ≥ 1 %) and repeats are masked **before** design;
  segmental duplications and genome-wide isPcr products are flagged.

An amplicon is **OK** only if it is clean of SNPs and repeats, designed with normal settings,
outside segmental duplications and (when checked) specific. Otherwise it is **Check**, with the
reasons listed. Every rule above has a flag that changes it — see
[Design method](docs/method.md).

## Run locally

```bash
pip install -r requirements.txt
python exonprimer_batch.py genes.txt -o exon_primers.xlsx            # Excel only
python exonprimer_batch.py genes.txt -o site/exon_primers.xlsx --json-dir site/data
python -m http.server -d site 8000                                   # preview the site
python exonprimer_batch.py --help                                    # all options
```

Ensembl and UCSC responses are cached in `.cache/`, so reruns with new parameters are fast.
GRCh37, genome-wide in-silico PCR and the GitHub workflow are in
[Building and publishing](docs/build-and-deploy.md).

## Excel output

| Sheet | Content |
|---|---|
| Primers | One row per amplicon: status (OK / CHECK / FAIL), names (`GENE_E3F/R`), sequences, length, **both Tm values** (PCR conditions, and UCSC's formula), **ΔTm F–R**, GC, product size, distances to the exon, coordinates, design tier, SNPs under F and R, repeat bases, segmental duplication, in-silico PCR result, warnings, UCSC link |
| Alternatives | The best pair plus 2 alternatives per amplicon, each with SNPs and in-silico PCR result |
| Amplicons | Amplicon sequences, exons UPPER / introns lower |
| Genes | Transcript, RefSeq, coding exons, amplicons per gene |
| Problems | Genes not found, exons without primers |
| Parameters | All settings used |

## Which Tm to use at the bench

**`Tm PCR`.** It is `F Tm PCR (°C)` in the workbook and the plain `Tm °C` column on the site:
the value the design was actually selected on, under the conditions in the tube (50 mM Na⁺,
1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo). The usual rule of thumb — anneal a few degrees below the
lower primer's Tm, or touch down to it — applies to this column and to no other number on the
page. On the current build 97.9 % of best pairs have both primers inside the 57–63 °C window, so
one programme — 58–60 °C, or a 65 → 58 °C touchdown — runs essentially the whole panel.

`Tm UCSC` is there for one purpose: so that a primer pasted into UCSC in-silico PCR gives a
number you recognise. It reads 5–8 °C higher on a 20-mer because UCSC keeps the older Breslauer
tables — a difference of convention, not a disagreement about the primer. Never set a
temperature from it. [Tm](docs/tm.md) has the arithmetic, the ΔTm caveat, and why a third column
was published and then removed.

## Triaging the Check rows

A large panel produces thousands of **Check** amplicons, and the reasons carry very different
risk: an off-target product may amplify a paralog, while a repeat touching a primer's 5′ end is
usually harmless. In sheet order they look alike, so most of the reading goes on the benign
majority.

```bash
python scripts/triage_checks.py site/data                              # a local build
python scripts/triage_checks.py https://<user>.github.io/<repo>/data   # the published site
```

It counts the amplicons by their most serious reason — off-target, segmental duplication, no
product, SNP under a primer, relaxed settings, 5′ repeat — and then lists the ones where **no
alternative pair avoids the problem**. Those are the real redesign queue: everywhere else, the
Alternatives sheet already holds a pair that sidesteps the issue.

## Before ordering primers

- Start with `scripts/triage_checks.py`, worst reason first, rather than reading the Primers
  sheet in order.
- Review every **Check** row: SNPs near the 3′ end, repeats, segmental duplications, off-target
  products.
- For genes with known pseudogenes, a specific in-silico PCR result does not rule out
  mismatched priming on the paralog — consider long-range PCR or paralog-specific design.
- GC-rich first exons may need relaxed settings or manual design (flagged in the output).

## Development

Developed with AI assistance (Claude Code). The method, parameters and clinical interpretation
are the author's; every primer in the output should be reviewed before ordering, as above.
