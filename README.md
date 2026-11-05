# exonprimer-batch

### → **[jimouse.github.io/exonprimer-batch](https://jimouse.github.io/exonprimer-batch)**

Intronic PCR/Sanger primers for **every coding exon of a list of genes**, designed the
ExonPrimer way, exported to **Excel**, and published as a **GitHub Pages site** where you pick a
gene and an exon — or paste a variant — to get the precomputed primers.

The live site covers **4,471 genes / 60,677 amplicons** (GRCh38, MANE Select), each primer pair
checked against common SNPs, repeats, segmental duplications and genome-wide in-silico PCR.
Search a gene (`MTHFR`), a variant (`POLG:c.2243G>C`) or a position (`chr15:89327201`), and
download the whole panel as one Excel workbook.

> ExonPrimer itself (IHG, Helmholtz Munich) is a web form that does not allow automated
> access, so this project reproduces its method locally with Primer3 instead of querying it.

## What it does

For each gene in `genes.txt`:

1. Gets the exon structure from Ensembl REST — the **MANE Select** transcript (RefSeq NM reported),
   falling back to the Ensembl canonical transcript, or a transcript you name.
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
4. Runs Primer3 (Tm 57–63 °C, opt 60; 18–27 nt, opt 20; GC 30–70 %; GC clamp). If nothing is
   found it retries once with relaxed settings and flags the result.

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

## Excel output

| Sheet | Content |
|---|---|
| Primers | One row per amplicon: status (OK / CHECK / FAIL), names (`GENE_E3F/R`), sequences, length, Tm, GC, product size, distances to the exon, coordinates, design tier, SNPs under F and R, repeat bases, segmental duplication, in-silico PCR result, warnings, UCSC link |
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
UCSC is queried at most once per second (`--ucsc-delay`).
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

- **A different transcript.** The panel is built on MANE Select. A matching accession with a
  different version is used but flagged; a different accession is refused rather than
  renumbered, because `c.` coordinates would not agree.
- **Deep intronic variants.** Primers reach about `--primer-region` bp into the flanks, so a
  variant further in is reported as covered by no amplicon, never rounded to the nearest exon.
- **UTR variants.** `c.-45` and `c.*120` lie outside the coding target and are refused.

When the variant states a reference allele it is checked against the amplicon sequence, which
catches a wrong genome build or a typo. Primers are matched on the **PCR product** span, not the
exon, so an intronic variant inside the product is correctly reported as sequenceable.

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
DOM and a stubbed `fetch` serving synthetic genes. It covers which gene a variant resolves to —
by symbol, Ensembl transcript, RefSeq, or coordinate — the transcript-safety rules (a version
difference warns, a different accession is refused), what the result panel renders, and that a
gene symbol is never routed as a variant.

`python tests/test_offline.py` builds a synthetic genome with genes on both strands and fake
Ensembl/UCSC/isPcr responses. It checks primer coordinates against the genome, minimum distance,
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
