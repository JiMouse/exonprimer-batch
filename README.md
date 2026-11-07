# exonprimer-batch

### → **[jimouse.github.io/exonprimer-batch](https://jimouse.github.io/exonprimer-batch)**

Intronic PCR/Sanger primers for **every coding exon of a list of genes**, designed the
ExonPrimer way, exported to **Excel**, and published as a **GitHub Pages site** where you pick a
gene and an exon — or paste a variant — to get the precomputed primers.

The gene list covers **6,437 genes** (GRCh38, MANE Select) — the original panel merged with the
**OMIM-morbid set**, every gene with an OMIM disease. Each primer pair is checked against common
SNPs, repeats, segmental duplications and genome-wide in-silico PCR.
Search a gene (`MTHFR`), a variant (`POLG:c.2243G>C`) or a position (`chr15:89327201`), and
download the whole panel as one Excel workbook.

> ExonPrimer itself (IHG, Helmholtz Munich) is a web form that does not allow automated
> access, so this project reproduces its method locally with Primer3 instead of querying it.

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

## Tm, and why two tools disagree

The Tm columns are Primer3's nearest-neighbour value (SantaLucia 1998). What matters is the
conditions, which are now stated rather than inherited from a library default:

| | default | flag |
|---|---|---|
| monovalent cations (Na⁺) | 50 mM | `--salt` |
| divalent cations (Mg²⁺) | 1.5 mM | `--salt` |
| total dNTP | 0.6 mM | `--salt` |
| annealing oligo | 50 nM | `--salt` |

These were previously left to `primer3-py`'s defaults, which is why the reported Tm is a
**PCR-condition** Tm. The Mg²⁺ term is worth about **6 °C** on a 20-mer: `PRKCG_E5F`
`ACTGTGGCTCCCTCCTCTAC` reads **60.3 °C** here and **54.4 °C** with Mg²⁺ and dNTP at zero, from
the very same tables. Most simple web calculators model Na⁺ alone, so they will read about 6 °C
lower — that is a difference of convention, not a disagreement about the primer. Both values are
now reported side by side (`F Tm PCR` / `F Tm Na+`) so no one has to guess.

The nearest-neighbour **ΔTm** between the two primers of a pair is stable across conventions.
Simple GC rules are not: for the pair above they give a ΔTm of 2.0 °C where both
nearest-neighbour models give 0.4 °C.

### Matching the two primers of a pair

`--max-tm-diff` (default **2 °C**) is the ceiling, and `--wt-tm-diff` (default **1**) is how hard
Primer3 prefers a matched pair. The weight matters more than the ceiling: **Primer3's own default
for it is 0**, so without it a matched pair is merely *allowed*, never *preferred* — which is why
10 % of the previous build's best pairs exceeded 2 °C.

Both levers change which pair wins, so they change primers. Measured on 179 amplicons across 10
genes, against the previous settings (ceiling 3, weight 0):

| settings | primers unchanged | amplicons with ΔTm > 2 °C | median ΔTm |
|---|---|---|---|
| ceiling 3, weight 0 (before) | — | 18 | 0.5 |
| ceiling 2, weight 0 | 91 % | 7 | 0.5 |
| ceiling 2, weight 1 (now) | 39 % | **0** | 0.2 |

The count of OK amplicons is the same in all three, so nothing is lost either way — it is a
choice between matched pairs and leaving a published panel alone. `--wt-tm-diff 0` keeps 91 % of
the primers and still removes most of the wide pairs, if that trade suits you better.

## Excel output

| Sheet | Content |
|---|---|
| Primers | One row per amplicon: status (OK / CHECK / FAIL), names (`GENE_E3F/R`), sequences, length, **both Tm values** (PCR conditions and Na⁺ only), **ΔTm F–R**, GC, product size, distances to the exon, coordinates, design tier, SNPs under F and R, repeat bases, segmental duplication, in-silico PCR result, warnings, UCSC link |
| Alternatives | The best pair plus 2 alternatives per amplicon, each with SNPs and in-silico PCR result |
| Amplicons | Amplicon sequences, exons UPPER / introns lower |
| Genes | Transcript, RefSeq, coding exons, amplicons per gene |
| Problems | Genes not found, exons without primers |
| Parameters | All settings used |

## Run locally

```bash
pip install -r requirements.txt
python exonprimer_batch.py genes.txt -o exon_primers.xlsx            # Excel only
python exonprimer_batch.py genes.txt -o site/exon_primers.xlsx --json-dir site/data
python -m http.server -d site 8000                                   # preview the site
python exonprimer_batch.py --help                                    # all options
```

Ensembl and UCSC responses are cached in `.cache/`, so reruns with new parameters are fast.
UCSC is queried at most once per second (`--ucsc-delay`), which is the rate UCSC asks for.
`--assembly GRCh37` uses the GRCh37 Ensembl server and hg19 annotations (no MANE there,
canonical is used).

Genome-wide in-silico PCR locally. isPcr lives in the `blat/` subdirectory of the UCSC
binaries; replace `linux.x86_64` with `macOSX.arm64` (Apple Silicon) or `macOSX.x86_64` (Intel):

```bash
mkdir -p tools
curl -o tools/isPcr https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/isPcr && chmod +x tools/isPcr
curl -o tools/hg38.2bit https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit   # ~800 MB
python exonprimer_batch.py genes.txt --ispcr tools/isPcr --genome-2bit tools/hg38.2bit
```

isPcr indexes the whole genome in memory (several GB of RAM). If it fails, the run continues
and the output says "in-silico PCR not run". On macOS the downloaded binary is unsigned, so
Gatekeeper may need `xattr -d com.apple.quarantine tools/isPcr`.

### How long a build takes

Designing is not the cost — a fully cached run does 5 genes in 0.86 s at 97 % CPU, so the 6,437
genes are about **18 minutes of computation**. Everything else is waiting on public APIs, which
is why the pipeline is shaped around cutting the *number* of requests rather than their size.
Measured per gene, in series:

| | before | now |
|---|---|---|
| exon structure | 1 Ensembl `/lookup/symbol` per gene, 6.4 s | UCSC `mane`, **24 requests for the whole genome** |
| soft-masked sequence | 1 Ensembl GET per gene, 14.1 s | `POST /sequence/region`, 50 regions at a time, **0.88 s** |
| SNPs / repeats / segdups | 3 UCSC requests per gene, 2.1 s + 3 s of mandated delay | unchanged |
| Primer3 and output | 0.17 s | unchanged |

Ensembl's cost is per-request overhead, not payload: the same endpoint took 19.4 s for a 25 KB
sequence and 6.7 s for a 218 KB one. Batching is the whole fix — 8,942 Ensembl requests become
about 115, which is also the courteous thing to do to a public server.

The UCSC annotation is deliberately left alone. UCSC's API documentation asks for *"a maximum of
one hit every second"* and runs a `botDelay` system that can impose *"a restriction of access"*;
if it throttled us, `annotate()` would carry on with the annotation missing, every amplicon of
that gene would be flagged `snps not checked` and demoted to **CHECK**, and the build would still
publish. Trading an hour for that risk is a bad deal. Fetching the tracks as bigBed files instead
is UCSC's own recommendation for bulk use, but querying `dbSnp155Common.bb` (1.96 GB) by remote
range costs 2.0 s per gene — no better than the API — so it would only pay off by downloading and
caching the file, which is a separate piece of work.

In the sharded workflow (4 concurrent jobs) that is roughly **2 h 45 for 6,437 genes**, against
about 46 h before.

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

### Primer pairs matched on Tm

The reported Tm is Primer3's nearest-neighbour value (SantaLucia 1998) under **PCR conditions**:
50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM primer. A calculator that models Na⁺ alone reads
about **6 °C lower for the same oligo** — so check which convention a tool uses before concluding
that two tools disagree. The nearest-neighbour ΔTm between the two primers of a pair is stable
across conventions; simple GC rules are not.

Every amplicon view shows **ΔTm** for the pair, and the checkbox above the search box re-orders
the alternatives already published so that a pair matched to within 2 °C comes first. It changes
nothing in the data — Primer3's own rank 1 stays listed, and a pair is never promoted if it would
trade away genome-wide specificity. On the current build, 10 % of best pairs exceed 2 °C and
two thirds of those already have a matched alternative.

## Publish the website on GitHub

1. Create a repository and push this folder to its `main` branch.
2. **Settings → Pages → Build and deployment → Source: GitHub Actions.**
3. Edit `genes.txt` and push — or run the workflow by hand from the **Actions** tab
   (you can choose the assembly and pass extra options there).
4. The site appears at `https://<user>.github.io/<repo>/` — for this repository,
   <https://jimouse.github.io/exonprimer-batch>. Links like `…/#MTHFR/5` open a gene and exon
   directly, and `…/#v=POLG:c.2243G%3EC` opens a variant. The full Excel is downloadable from
   the site and from each workflow run's artifacts.

The workflow runs the offline tests, designs primers, and deploys `site/` — nothing generated is
committed to the repository.

### Large gene lists

A few thousand genes do not fit in one GitHub job (6 h limit), so the build is **sharded**:

| Job | What it does |
|---|---|
| `plan` | runs the offline tests, strips comments from `genes.txt` and `split`s it into chunks |
| `design` | a matrix job per chunk, each designing its own genes and uploading `part.xlsx` + `data/` |
| `merge` | `scripts/merge_parts.py` stitches the chunks into one `site/` and one workbook |
| `deploy` | publishes the merged `site/` to Pages |

Tune it from the **Actions** tab: `chunk_size` (genes per job, default 150) and `max_parallel`
(concurrent jobs, default 4). **Keep `max_parallel` low** — every job queries Ensembl and UCSC,
and parallel jobs multiply that load on public servers.

### Changing the website without rebuilding the data

The built data is also published as an asset on the `data-cache` prerelease, named after a
fingerprint of the inputs that determine it — `genes.txt`, `exonprimer_batch.py`,
`scripts/merge_parts.py`, the assembly, the isPcr flag and any extra options. The website is
deliberately **not** part of that fingerprint.

So editing `site/index.html` leaves the fingerprint unchanged: `plan` finds the published asset,
the 30 `design` jobs are skipped, `merge` restores the data and the site redeploys in a couple of
minutes. Editing `genes.txt` changes the fingerprint and the data is rebuilt as usual. Tick
**force_rebuild** in the Actions tab to rebuild regardless.

The three newest fingerprints are kept and older assets pruned. A partial build is never
published, so a missing chunk can't be baked in and restored by later website changes. The
`data-cache` release is not a software release — if you archive releases to Zenodo, exclude it.

Each chunk caches its own Ensembl/UCSC responses, so a re-run only pays for genes that are
genuinely new. `fail-fast` is off and `merge` runs even if some chunks fail: the site is built
from whatever succeeded, and the run log warns how many chunks are missing. Re-running the
workflow fills the gaps from cache.

The merged workbook is written with openpyxl's write-only writer, so memory stays flat
regardless of how many amplicons there are.

Expect retired gene symbols to fail: Ensembl does not follow HGNC aliases, so `GBA` or `AARS`
resolve to nothing and land in the Problems sheet. Use current symbols.

## Tests

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
that never promotes a pair with an off-target product, and that a gene symbol is never routed as
a variant.

`python tests/test_offline.py` builds a synthetic genome with genes on both strands and fake
Ensembl/UCSC/isPcr responses. The fakes serve the UCSC `mane` track and the batched
`POST /sequence/region` too, so the fast paths are the ones under test: that a model built from
the MANE track is **identical** to one built from Ensembl (exons, transcript, CDS bounds, region),
that a gene missing from the track falls back to Ensembl rather than being skipped, that one POST
covers many genes and the sequences match the single-region ones, and that a dead batch endpoint
falls back instead of losing data. It also pins both Tm values, the ΔTm ceiling, and that the salt
conditions actually reach Primer3. Beyond that it checks primer coordinates against the genome,
minimum distance,
coverage of every coding base, SNP avoidance and the 3′-end rule, the MAF cut-off, repeat
avoidance and the 5′-only fallback, segmental-duplication alerts, and that an off-target
in-silico PCR product demotes a pair.

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
