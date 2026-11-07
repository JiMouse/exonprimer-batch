# Building and publishing

Local options beyond the quickstart, what a full build costs, and how the sharded
GitHub workflow publishes the site. [← back to the README](../README.md)

## Local options

Ensembl and UCSC responses are cached in `.cache/`, so reruns with new parameters are fast.
UCSC is queried at most once per second (`--ucsc-delay`), which is the rate UCSC asks for.
`--assembly GRCh37` uses the GRCh37 Ensembl server and hg19 annotations (no MANE there,
canonical is used).

## Genome-wide in-silico PCR locally

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

## How long a build takes

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

## Publish the website on GitHub

1. Create a repository and push this folder to its `main` branch.
2. **Settings → Pages → Build and deployment → Source: GitHub Actions.**
3. Edit `genes.txt` and push — or run the workflow by hand from the **Actions** tab
   (you can choose the assembly and pass extra options there).
4. The site appears at `https://<user>.github.io/<repo>/` — for this repository,
   <https://jimouse.github.io/exonprimer-batch>. Links like `…/#MTHFR/5` open a gene and exon
   directly, and `…/#v=POLG:c.2243G%3EC` opens a variant. The full Excel is downloadable from
   the site and from each workflow run's artifacts.

The workflow runs the offline tests, designs primers, deploys `site/` and commits the built data
back to `main` as `Built data: N files changed`, so `git log -p site/data/genes/POLG.json` is the
history of that gene's primers. The workbook is not committed — it is a build artifact and a
download on the site.

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
the 43 `design` jobs are skipped, `merge` restores the data and the site redeploys in a couple of
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

