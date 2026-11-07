#!/usr/bin/env python3
"""
exonprimer_batch.py - ExonPrimer-style intronic PCR primer design for every
coding exon of a list of genes, with Excel and JSON (website) output.

Method (ExonPrimer, IHG Helmholtz Munich; safeguards as in EasyExonPrimer,
Wu & Munroe, Appl Bioinformatics 2006):
  * Exon structure: Ensembl REST, MANE Select transcript (falls back to the
    Ensembl canonical transcript, or a transcript you name in the gene list).
  * Each coding exon is a PCR target. Primers sit in the flanking sequence,
    at least --min-dist bp from the exon/intron boundary, within a
    --primer-region bp window beyond that.
  * PCR products are 200-500 bp by default (--product). Exons separated by
    a small intron are combined when they fit; exons too long for one
    product are split into overlapping fragments, each with its own pair.
  * Common SNPs (dbSNP 155 "Common", 1000 Genomes MAF >= 1%, UCSC) and
    repeats (RepeatMasker / Ensembl soft-masking) are masked before Primer3,
    in tiers: first no SNP and no repeat anywhere under the primers; then a
    repeat may touch the 5' part of a primer; then a SNP may sit outside the
    3'-terminal --snp-3prime bases. Every remaining SNP/repeat is reported.
  * High homology: targets in segmental duplications (UCSC genomicSuperDups)
    are flagged, and, when isPcr and the genome .2bit are given, every
    candidate pair is checked genome-wide by in-silico PCR; specific pairs
    are preferred and off-target products are reported.

Usage:
  python exonprimer_batch.py genes.txt -o primers.xlsx --json-dir site/data \
      [--ispcr tools/isPcr --genome-2bit tools/hg38.2bit]

genes.txt: one HGNC symbol per line; optionally a transcript ID after it
(e.g. "BRCA1  ENST00000357654") to override the MANE Select choice.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlencode

import primer3
import requests

__version__ = "2.1.0"

SERVERS = {
    "GRCh38": "https://rest.ensembl.org",
    "GRCh37": "https://grch37.rest.ensembl.org",
}
UCSC_DB = {"GRCh38": "hg38", "GRCh37": "hg19"}
UCSC_API = "https://api.genome.ucsc.edu"
SNP_TRACK, RMSK_TRACK, SEGDUP_TRACK = "dbSnp155Common", "rmsk", "genomicSuperDups"
MANE_TRACK = "mane"
# Chromosomes the MANE track is fetched for, one request each.
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
TIERS = {
    "clean": "no common SNP and no repeat under either primer",
    "repeat-5p": "a repeat overlaps the 5' part of a primer (3' end outside repeats)",
    "snp-5p": "a common SNP lies under a primer, outside its 3'-terminal bases",
}


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def tm_na(seq: str) -> float:
    """Nearest-neighbour Tm with monovalent salt only (no Mg2+, no dNTP)."""
    return round(primer3.calc_tm(seq.upper(), mv_conc=50.0, dv_conc=0.0,
                                 dntp_conc=0.0, dna_conc=50.0), 1)


def ucsc_chrom(chrom: str) -> str:
    return "chrM" if chrom in ("MT", "M") else f"chr{chrom}"


class _Http:
    def __init__(self, cache_dir: str | None, delay: float, name: str):
        self.session = requests.Session()
        self.session.headers["User-Agent"] = f"exonprimer-batch/{__version__}"
        self.cache_dir = Path(cache_dir) / name if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.delay = delay

    def post_json(self, path_or_url: str, body: dict, key: str, params=None):
        """POST with the same cache, retry and backoff as get_json."""
        url = path_or_url if path_or_url.startswith("http") else getattr(self, "base", "") + path_or_url
        return self._request("POST", url, key, params, json.dumps(body))

    def get_json(self, url: str, key: str, params=None, headers=None):
        return self._request("GET", url, key, params, None, headers)

    def _request(self, method: str, url: str, key: str, params=None, body=None, headers=None):
        cache = None
        if self.cache_dir:
            cache = self.cache_dir / (re.sub(r"[^A-Za-z0-9._-]+", "_", key)[:200] + ".json")
            if cache.exists():
                return json.loads(cache.read_text())
        for attempt in range(8):
            try:
                hdr = dict(headers or {})
                if body is not None:
                    hdr.setdefault("Content-Type", "application/json")
                    hdr.setdefault("Accept", "application/json")
                r = self.session.request(method, url, params=params, headers=hdr,
                                         data=body, timeout=600 if body else 120)
            except requests.RequestException:
                time.sleep(2 ** attempt)
                continue
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 2 ** attempt)))
                continue
            if r.status_code in (500, 502, 503, 504):
                time.sleep(2 ** attempt)
                continue
            if r.status_code in (400, 404):
                try:
                    msg = r.json().get("error", r.text)
                except ValueError:
                    msg = r.text
                raise LookupError(msg)
            r.raise_for_status()
            data = r.json()
            if cache:
                cache.write_text(json.dumps(data))
            time.sleep(self.delay)
            return data
        raise RuntimeError(f"request failed repeatedly: {url}")


# --------------------------------------------------------------------------
# Ensembl (exon structure + soft-masked sequence)
# --------------------------------------------------------------------------
class Ensembl(_Http):
    def __init__(self, assembly: str = "GRCh38", cache_dir: str | None = None):
        super().__init__(cache_dir, 0.07, "ensembl")  # < 15 requests/s
        self.base = SERVERS[assembly]

    def _get(self, path, params=None):
        return self.get_json(self.base + path, path + "_" + urlencode(params or {}), params,
                             {"Content-Type": "application/json"})

    def gene(self, symbol):
        return self._get(f"/lookup/symbol/homo_sapiens/{symbol}", {"expand": 1, "mane": 1})

    def transcript(self, transcript_id):
        return self._get(f"/lookup/id/{transcript_id}", {"expand": 1, "mane": 1})

    def sequence(self, chrom, start, end, strand):
        return self._get(f"/sequence/region/human/{chrom}:{start}..{end}:{strand}", {"mask": "soft"})["seq"]

    def sequences(self, regions: list[str]) -> dict[str, str]:
        """Soft-masked sequence for several regions in one POST.

        Keyed by the `query` string Ensembl echoes for each result, so a
        reordered or short response can never silently misalign sequences
        with genes. At most 50 regions per request, which is Ensembl's limit.
        """
        if not regions:
            return {}
        data = self.post_json("/sequence/region/human", {"regions": regions},
                              "seqbatch_" + hashlib.sha1("|".join(regions).encode()).hexdigest(),
                              {"mask": "soft"})
        out = {}
        for item in data:
            q = item.get("query")
            if q and item.get("seq"):
                out[q] = item["seq"]
        return out


# --------------------------------------------------------------------------
# UCSC (common SNPs, RepeatMasker, segmental duplications)
# --------------------------------------------------------------------------
class UCSC(_Http):
    def __init__(self, assembly: str = "GRCh38", cache_dir: str | None = None, delay: float = 1.0):
        super().__init__(cache_dir, delay, "ucsc")
        self.db = UCSC_DB[assembly]

    def track(self, track: str, chrom: str, start: int, end: int) -> list[dict]:
        """Items of a track overlapping chrom:start-end (1-based inclusive)."""
        q = f"genome={self.db};track={track};chrom={chrom};start={start - 1};end={end};maxItemsOutput=-1"
        data = self.get_json(f"{UCSC_API}/getData/track?{q}", f"{self.db}_{track}_{chrom}_{start}_{end}")
        items = data.get(track, [])
        if isinstance(items, dict):  # some tracks return {chrom: [items]}
            items = [x for v in items.values() for x in v]
        return items

    def mane(self, chrom: str) -> list[dict]:
        """Every MANE transcript on a chromosome, in one request."""
        q = f"genome={self.db};track={MANE_TRACK};chrom={chrom};maxItemsOutput=-1"
        data = self.get_json(f"{UCSC_API}/getData/track?{q}", f"{self.db}_{MANE_TRACK}_{chrom}")
        items = data.get(MANE_TRACK, [])
        if isinstance(items, dict):
            items = [x for v in items.values() for x in v]
        return items

    def snps(self, chrom, start, end, min_maf: float) -> list[dict]:
        out = []
        for it in self.track(SNP_TRACK, chrom, start, end):
            maf = max_maf(it.get("minorAlleleFreq"))
            if maf is not None and maf < min_maf:
                continue
            s0, e0 = int(it["chromStart"]), int(it["chromEnd"])
            out.append({"rs": it.get("name", ""), "start": s0 + 1, "end": max(e0, s0 + 1),
                        "maf": None if maf is None else round(maf, 4),
                        "class": it.get("class", ""), "ref": it.get("ref", ""),
                        "alts": it.get("alts", "")})
        return out

    def repeats(self, chrom, start, end) -> list[dict]:
        return [{"start": int(it["genoStart"]) + 1, "end": int(it["genoEnd"]),
                 "name": it.get("repName", ""), "class": it.get("repClass", "")}
                for it in self.track(RMSK_TRACK, chrom, start, end)]

    def segdups(self, chrom, start, end) -> list[dict]:
        return [{"start": int(it["chromStart"]) + 1, "end": int(it["chromEnd"]),
                 "other": f'{it.get("otherChrom", "")}:{int(it.get("otherStart", 0)) + 1}-{it.get("otherEnd", "")}',
                 "identity": round(100 * float(it.get("fracMatch", 0)), 1)}
                for it in self.track(SEGDUP_TRACK, chrom, start, end)]


def max_maf(v) -> float | None:
    """Highest minor-allele frequency across dbSNP frequency sources."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        vals = [float(v)]
    else:
        parts = v.split(",") if isinstance(v, str) else v
        vals = []
        for x in parts:
            try:
                vals.append(float(x))
            except (TypeError, ValueError):
                pass
    vals = [x for x in vals if 0 <= x <= 1]
    return max(vals) if vals else None


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------
@dataclass
class Params:
    min_dist: int = 50
    primer_region: int = 150
    product_min: int = 200      # PCR product length range (bp)
    product_max: int = 500
    max_target: int = 0         # 0 = derived from product_max (see __post_init__)
    overlap: int = 50
    combine: bool = True
    target: str = "coding"
    tm_min: float = 57.0
    tm_opt: float = 60.0
    tm_max: float = 63.0
    size_min: int = 18
    size_opt: int = 20
    size_max: int = 27
    gc_min: float = 30.0
    gc_max: float = 70.0
    gc_clamp: int = 1
    max_tm_diff: float = 2.0    # |Tm(F) - Tm(R)| ceiling for a pair
    wt_tm_diff: float = 1.0     # how hard Primer3 prefers a matched pair
    # Thermodynamics. These were previously left to primer3-py's defaults,
    # which is why the reported Tm is a PCR-condition Tm: the 1.5 mM Mg2+ is
    # worth about 6 degrees on a 20-mer, so a calculator that models Na+
    # alone reads that much lower for the very same oligo.
    salt_mono: float = 50.0     # monovalent cations, mM
    salt_div: float = 1.5       # divalent cations (Mg2+), mM
    dntp: float = 0.6           # total dNTP, mM
    dna_conc: float = 50.0      # annealing oligo concentration, nM
    max_poly_x: int = 4
    num_return: int = 3
    relax: bool = True
    # safeguards
    snp_maf: float = 0.01       # SNPs at or above this MAF are avoided
    snp_3prime: int = 5         # never accept a SNP in this many 3'-terminal bases
    fallback: bool = True       # allow tiers repeat-5p / snp-5p when "clean" fails
    ispcr_max_size: int = 4000  # in-silico PCR: max product size considered

    def __post_init__(self):
        # Largest target that still fits in product_max with primers at
        # min_dist plus a 20-bp margin on each side for primer placement.
        if not 0 < self.product_min <= self.product_max:
            raise ValueError("product range must be MIN <= MAX (bp)")
        if self.max_target <= 0:
            self.max_target = self.product_max - 2 * (self.min_dist + self.size_opt + 20)
        if self.max_target < 30 or self.max_target + 2 * (self.min_dist + self.size_min) > self.product_max:
            raise ValueError(f"max_target {self.max_target} bp cannot fit in a {self.product_max}-bp product "
                             f"with primers {self.min_dist} bp from the exon")
        if self.overlap >= self.max_target:
            raise ValueError("overlap must be smaller than max_target")

    def primer3_global(self, relaxed: bool = False) -> dict:
        g = {
            "PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0,
            "PRIMER_OPT_SIZE": self.size_opt, "PRIMER_MIN_SIZE": self.size_min, "PRIMER_MAX_SIZE": self.size_max,
            "PRIMER_OPT_TM": self.tm_opt, "PRIMER_MIN_TM": self.tm_min, "PRIMER_MAX_TM": self.tm_max,
            "PRIMER_MIN_GC": self.gc_min, "PRIMER_MAX_GC": self.gc_max, "PRIMER_GC_CLAMP": self.gc_clamp,
            "PRIMER_PAIR_MAX_DIFF_TM": self.max_tm_diff, "PRIMER_MAX_POLY_X": self.max_poly_x,
            "PRIMER_MAX_NS_ACCEPTED": 0, "PRIMER_EXPLAIN_FLAG": 1,
            # Primer3's own default for this weight is 0, so without it a
            # matched pair is merely allowed, never preferred.
            "PRIMER_PAIR_WT_DIFF_TM": self.wt_tm_diff,
            "PRIMER_SALT_MONOVALENT": self.salt_mono, "PRIMER_SALT_DIVALENT": self.salt_div,
            "PRIMER_DNTP_CONC": self.dntp, "PRIMER_DNA_CONC": self.dna_conc,
        }
        if relaxed:
            g.update({
                "PRIMER_MIN_TM": self.tm_min - 3, "PRIMER_MAX_TM": self.tm_max + 3,
                "PRIMER_MIN_SIZE": max(16, self.size_min - 2), "PRIMER_MAX_SIZE": self.size_max + 3,
                "PRIMER_MIN_GC": max(15.0, self.gc_min - 10), "PRIMER_MAX_GC": min(85.0, self.gc_max + 10),
                "PRIMER_GC_CLAMP": 0, "PRIMER_PAIR_MAX_DIFF_TM": self.max_tm_diff + 1,
                "PRIMER_MAX_POLY_X": self.max_poly_x + 1,
            })
        return g


@dataclass
class Exon:
    number: int
    start: int
    end: int
    cds_start: int | None
    cds_end: int | None
    ensembl_id: str = ""

    @property
    def coding(self) -> bool:
        return self.cds_start is not None


@dataclass
class Target:
    exons: list[int]
    lo: int
    hi: int
    label: str = ""
    fragment: str = ""


@dataclass
class GeneModel:
    gene: str
    gene_id: str
    transcript_id: str
    transcript_version: str
    refseq: str
    selection: str
    chrom: str
    strand: int
    exons: list[Exon] = field(default_factory=list)
    region_start: int = 0
    region_end: int = 0
    seq: str = ""                                   # transcript-oriented, soft-masked
    snps: list[dict] = field(default_factory=list)  # genomic, from UCSC
    repeats: list[dict] = field(default_factory=list)
    segdups: list[dict] = field(default_factory=list)
    annotated: dict = field(default_factory=dict)   # which annotations were loaded
    rep: bytearray = field(default_factory=bytearray)      # local repeat mask
    snp_at: dict = field(default_factory=dict)             # local pos -> [snp]

    def to_local(self, g: int) -> int:
        return g - self.region_start if self.strand == 1 else self.region_end - g

    def to_genomic(self, i: int) -> int:
        return self.region_start + i if self.strand == 1 else self.region_end - i

    def local_span(self, a: int, b: int) -> tuple[int, int]:
        x, y = self.to_local(a), self.to_local(b)
        return (x, y) if x <= y else (y, x)

    def clip_span(self, a: int, b: int) -> tuple[int, int] | None:
        a, b = max(a, self.region_start), min(b, self.region_end)
        return self.local_span(a, b) if a <= b else None

    def build_masks(self):
        self.rep = bytearray(1 if c.islower() else 0 for c in self.seq)
        for r in self.repeats:
            sp = self.clip_span(r["start"], r["end"])
            if sp:
                self.rep[sp[0]:sp[1] + 1] = b"\x01" * (sp[1] - sp[0] + 1)
        self.snp_at = {}
        for s in self.snps:
            sp = self.clip_span(s["start"], s["end"])
            if sp:
                for i in range(sp[0], sp[1] + 1):
                    self.snp_at.setdefault(i, []).append(s)


# --------------------------------------------------------------------------
# Build the gene model
# --------------------------------------------------------------------------
def pick_transcript(gene: dict) -> tuple[dict, str]:
    txs = gene.get("Transcript", [])
    for t in txs:
        for m in t.get("MANE", []) or []:
            if m.get("type") == "MANE_Select":
                return t, "MANE Select"
    for t in txs:
        if t.get("is_canonical") and t.get("Translation"):
            return t, "Ensembl canonical"
    coding = [t for t in txs if t.get("Translation")]
    if coding:
        return max(coding, key=lambda t: t["Translation"]["length"]), "longest CDS"
    raise LookupError("no protein-coding transcript")


class ManeTable:
    """MANE Select transcripts for the whole genome, from the UCSC `mane` track.

    One request per chromosome replaces one Ensembl /lookup/symbol per gene -
    24 requests instead of thousands, which is the difference between a build
    that takes hours and one that takes minutes. The track carries everything
    build_model needs: the versioned ENST, the RefSeq NM, the CDS bounds and
    the exon blocks.

    Checked against the 4,471 genes of the previous build: exons, transcript
    and CDS bounds all identical for 4,449 (99.51 %). The rest fall back to
    Ensembl, which is also where genes with no MANE transcript at all go, so
    the transcript actually used never changes silently.
    """

    def __init__(self, ucsc, log=print):
        self.by_symbol: dict[str, dict] = {}
        if ucsc is None:
            return
        for chrom in CHROMS:
            try:
                rows = ucsc.mane(f"chr{chrom}")
            except Exception as ex:
                log(f"warning: MANE track for chr{chrom} not loaded ({ex}); those genes use Ensembl")
                continue
            for r in rows:
                if r.get("maneStat") != "MANE Select":
                    continue
                sym = r.get("geneName2") or ""
                if sym:
                    self.by_symbol.setdefault(sym.upper(), r)
        log(f"MANE table: {len(self.by_symbol)} genes from {len(CHROMS)} UCSC requests")

    def model(self, symbol: str, pad: int) -> GeneModel | None:
        r = self.by_symbol.get(symbol.upper())
        if not r:
            return None
        strand = -1 if r.get("strand") == "-" else 1
        cds_lo, cds_hi = r["thickStart"] + 1, r["thickEnd"]       # 1-based inclusive
        if cds_hi < cds_lo:
            return None                                          # no CDS: let Ensembl decide
        st = r["chromStart"]
        sizes = [int(v) for v in str(r["blockSizes"]).rstrip(",").split(",") if v != ""]
        offs = [int(v) for v in str(r["chromStarts"]).rstrip(",").split(",") if v != ""]
        if not sizes or len(sizes) != len(offs):
            return None
        blocks = sorted(((st + o + 1, st + o + sz) for o, sz in zip(offs, sizes)),
                        reverse=(strand == -1))
        exons = []
        for i, (lo, hi) in enumerate(blocks, 1):
            cs, ce = max(lo, cds_lo), min(hi, cds_hi)
            exons.append(Exon(i, lo, hi, cs if cs <= ce else None, ce if cs <= ce else None))
        if not any(e.coding for e in exons):
            return None
        m = GeneModel(gene=symbol, gene_id=str(r.get("geneName", "")), transcript_id=str(r["name"]).split(".")[0],
                      transcript_version=str(r["name"]), refseq=str(r.get("ncbiId", "")),
                      selection="MANE Select", chrom=str(r["chrom"]).replace("chr", ""),
                      strand=strand, exons=exons)
        m.region_start = max(1, min(e.start for e in exons) - pad)
        m.region_end = max(e.end for e in exons) + pad
        return m


def fetch_sequences(ens, models: list[GeneModel], batch: int = 50, log=print):
    """Soft-masked sequence for many models in POST batches.

    Ensembl's per-request overhead dominates a GET per gene (14 s measured,
    independent of payload size), while POST /sequence/region takes 50 regions
    at once for 0.88 s per gene - the same server, the same mask=soft, and the
    sequences are byte-identical to the individual GETs. Results are matched
    by the `query` field they echo back, never by position.

    A model whose sequence does not come back keeps seq == "" and is reported
    by the caller, so a short read can never be mistaken for a short gene.
    """
    todo = [m for m in models if not m.seq]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        regions = [f"{m.chrom}:{m.region_start}..{m.region_end}:{m.strand}" for m in chunk]
        try:
            got = ens.sequences(regions)
        except Exception as ex:
            log(f"warning: batch sequence fetch failed ({ex}); falling back to one request per gene")
            got = {}
            for m, q in zip(chunk, regions):
                try:
                    got[q] = ens.sequence(m.chrom, m.region_start, m.region_end, m.strand)
                except Exception as ex2:
                    log(f"    warning: {m.gene} sequence not loaded ({ex2})")
        for m, q in zip(chunk, regions):
            seq = got.get(q)
            if not seq:
                continue
            m.seq = seq
            if len(m.seq) != m.region_end - m.region_start + 1:
                m.region_end = m.region_start + len(m.seq) - 1
        log(f"  sequences {min(i + batch, len(todo))}/{len(todo)}")


def build_model(ens, symbol: str, transcript_id: str | None, pad: int,
                with_seq: bool = True) -> GeneModel:
    """Exon structure from Ensembl. With with_seq=False the sequence is left
    for fetch_sequences() to collect in batches."""
    gene = ens.gene(symbol)
    if transcript_id:
        tx = ens.transcript(transcript_id.split(".")[0])
        selection = "user"
        if not tx.get("Translation"):
            raise LookupError(f"{transcript_id} is not protein-coding")
    else:
        tx, selection = pick_transcript(gene)
    refseq = next((m.get("refseq_match", "") for m in tx.get("MANE", []) or []
                   if m.get("type") == "MANE_Select"), "")
    strand = tx["strand"]
    tl = tx["Translation"]
    exons = []
    for i, e in enumerate(sorted(tx["Exon"], key=lambda e: e["start"], reverse=(strand == -1)), 1):
        cs, ce = max(e["start"], tl["start"]), min(e["end"], tl["end"])
        exons.append(Exon(i, e["start"], e["end"], cs if cs <= ce else None, ce if cs <= ce else None, e.get("id", "")))
    m = GeneModel(gene=gene.get("display_name", symbol), gene_id=gene["id"], transcript_id=tx["id"],
                  transcript_version=f'{tx["id"]}.{tx.get("version", "")}'.rstrip("."), refseq=refseq,
                  selection=selection, chrom=str(tx["seq_region_name"]), strand=strand, exons=exons)
    m.region_start = max(1, min(e.start for e in exons) - pad)
    m.region_end = max(e.end for e in exons) + pad
    if with_seq:
        m.seq = ens.sequence(m.chrom, m.region_start, m.region_end, strand)
        if len(m.seq) != m.region_end - m.region_start + 1:
            m.region_end = m.region_start + len(m.seq) - 1
    return m


def annotate(m: GeneModel, ucsc, p: Params, log=print):
    """Load common SNPs, repeats and segmental duplications for the gene region."""
    chrom = ucsc_chrom(m.chrom)
    for kind, fn in (("snps", lambda: ucsc.snps(chrom, m.region_start, m.region_end, p.snp_maf)),
                     ("repeats", lambda: ucsc.repeats(chrom, m.region_start, m.region_end)),
                     ("segdups", lambda: ucsc.segdups(chrom, m.region_start, m.region_end))):
        if ucsc is None:
            m.annotated[kind] = False
            continue
        try:
            setattr(m, kind, fn())
            m.annotated[kind] = True
        except Exception as ex:
            m.annotated[kind] = False
            log(f"    warning: {m.gene} {kind} not loaded from UCSC ({ex})")
    m.build_masks()


# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------
def make_targets(m: GeneModel, p: Params) -> list[Target]:
    base = []
    for e in m.exons:
        if e.coding:
            a, b = (e.cds_start, e.cds_end) if p.target == "coding" else (e.start, e.end)
            base.append(Target([e.number], *m.local_span(a, b)))
    merged: list[Target] = []
    for t in base:
        if p.combine and merged and (t.hi - merged[-1].lo + 1) <= p.max_target:
            merged[-1] = Target(merged[-1].exons + t.exons, merged[-1].lo, t.hi)
        else:
            merged.append(t)
    out = []
    for t in merged:
        label = str(t.exons[0]) if len(t.exons) == 1 else f"{t.exons[0]}-{t.exons[-1]}"
        length = t.hi - t.lo + 1
        if length <= p.max_target:
            out.append(Target(t.exons, t.lo, t.hi, label))
            continue
        n = math.ceil((length - p.overlap) / (p.max_target - p.overlap))
        frag = math.ceil((length + (n - 1) * p.overlap) / n)
        for k in range(n):
            lo = t.lo + k * (frag - p.overlap)
            hi = min(t.hi, lo + frag - 1) if k < n - 1 else t.hi
            letter = chr(ord("a") + k) if k < 26 else f"_{k + 1}"
            out.append(Target(t.exons, lo, hi, f"{label}{letter}", f"{k + 1}/{n}"))
    return out


# --------------------------------------------------------------------------
# Primer design for one target (tiered masking)
# --------------------------------------------------------------------------
def masked_template(m: GeneModel, w0: int, w1: int, tier: str) -> str:
    """Template for Primer3. 'N' = forbidden under a primer; lowercase =
    repeat (Primer3 rejects primers whose 3' end is lowercase)."""
    out = []
    for i in range(w0, w1 + 1):
        c = m.seq[i].upper()
        rep, snp = m.rep[i], i in m.snp_at
        if snp and tier in ("clean", "repeat-5p"):
            c = "N"
        elif rep:
            c = "N" if tier == "clean" else c.lower()
        out.append(c)
    return "".join(out)


def design_target(m: GeneModel, t: Target, p: Params) -> dict:
    w0 = max(0, t.lo - p.min_dist - p.primer_region)
    w1 = min(len(m.seq) - 1, t.hi + p.min_dist + p.primer_region)
    left_len = (t.lo - p.min_dist) - w0
    right_start = (t.hi + p.min_dist + 1) - w0
    right_len = (w1 + 1) - (t.hi + p.min_dist + 1)
    result = {"pairs": [], "relaxed": False, "tier": "", "explain": ""}
    if left_len < p.size_min or right_len < p.size_min:
        result["explain"] = "not enough flanking sequence"
        return result
    ok = [[0, left_len, right_start, right_len]]
    min_prod = max(p.product_min, (t.hi - t.lo + 1) + 2 * p.min_dist)
    max_prod = min(p.product_max, w1 - w0 + 1)
    if min_prod > max_prod:
        result["explain"] = f"target too long for a {p.product_max}-bp product"
        return result
    n_cand = max(6, 2 * p.num_return)  # extra candidates for the specificity check
    tiers = list(TIERS) if p.fallback else ["clean"]
    for relaxed in ([False, True] if p.relax else [False]):
        for tier in tiers:
            g = p.primer3_global(relaxed)
            g["PRIMER_PRODUCT_SIZE_RANGE"] = [[min_prod, max_prod]]
            g["PRIMER_LOWERCASE_MASKING"] = 1
            g["PRIMER_NUM_RETURN"] = n_cand * (4 if tier == "snp-5p" else 1)
            res = primer3.design_primers(
                {"SEQUENCE_ID": f"{m.gene}_E{t.label}",
                 "SEQUENCE_TEMPLATE": masked_template(m, w0, w1, tier),
                 "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": ok}, g)
            pairs = [pair_record(m, t, res, i, w0) for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0))]
            if tier == "snp-5p":  # no SNP in the 3'-terminal bases
                pairs = [b for b in pairs
                         if all(s["pos3"] > p.snp_3prime for s in b["f_snps"] + b["r_snps"])]
            if pairs:
                pairs = pairs[:n_cand]
                for k, b in enumerate(pairs, 1):
                    b["rank"] = k
                result.update(pairs=pairs, relaxed=relaxed, tier=tier)
                return result
            if tier == "clean" and not relaxed:
                result["explain"] = "; ".join(
                    f"{k.split('_')[1].lower()}: {res.get(k, '')}"
                    for k in ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN"))
    return result


def primer_features(m: GeneModel, lo: int, hi: int, three_prime_hi: bool) -> tuple[list, int, list]:
    """SNPs (with distance from the 3' end, 1 = terminal base), repeat-masked
    bases and repeat names under a primer occupying local [lo, hi]."""
    snps, seen = [], set()
    for i in range(lo, hi + 1):
        for s in m.snp_at.get(i, []):
            if s["rs"] in seen:
                continue
            seen.add(s["rs"])
            snps.append({"rs": s["rs"], "maf": s["maf"], "pos3": (hi - i + 1) if three_prime_hi else (i - lo + 1)})
    rep_bases = sum(m.rep[lo:hi + 1])
    names = []
    if rep_bases:
        g1, g2 = sorted((m.to_genomic(lo), m.to_genomic(hi)))
        names = sorted({f'{r["name"]} ({r["class"]})' for r in m.repeats if r["start"] <= g2 and r["end"] >= g1})
    return sorted(snps, key=lambda s: s["pos3"]), rep_bases, names


def pair_record(m: GeneModel, t: Target, res: dict, i: int, w0: int) -> dict:
    ls, llen = res[f"PRIMER_LEFT_{i}"]
    rs, rlen = res[f"PRIMER_RIGHT_{i}"]  # rs = 3'-most base of the right primer
    a_lo, a_hi = w0 + ls, w0 + rs
    g1, g2 = sorted((m.to_genomic(a_lo), m.to_genomic(a_hi)))
    exon_mask = bytearray(a_hi - a_lo + 1)
    for e in m.exons:
        lo, hi = m.local_span(e.start, e.end)
        for j in range(max(lo, a_lo), min(hi, a_hi) + 1):
            exon_mask[j - a_lo] = 1
    amp = "".join(c.upper() if exon_mask[j] else c.lower() for j, c in enumerate(m.seq[a_lo:a_hi + 1]))
    f_lo, f_hi = a_lo, a_lo + llen - 1
    r_lo, r_hi = a_hi - rlen + 1, a_hi
    f_snps, f_rep, f_rnames = primer_features(m, f_lo, f_hi, True)
    r_snps, r_rep, r_rnames = primer_features(m, r_lo, r_hi, False)
    return {
        "rank": i + 1,
        "f_seq": m.seq[f_lo:f_hi + 1].upper(),
        "r_seq": rc(m.seq[r_lo:r_hi + 1]).upper(),
        "f_tm": round(res[f"PRIMER_LEFT_{i}_TM"], 1), "r_tm": round(res[f"PRIMER_RIGHT_{i}_TM"], 1),
        # The same nearest-neighbour tables with the divalent cations and
        # dNTPs removed - the convention most simple web calculators use.
        # Reported alongside so two tools can be compared without guessing
        # which conditions each assumed.
        "f_tm_na": tm_na(m.seq[f_lo:f_hi + 1]), "r_tm_na": tm_na(rc(m.seq[r_lo:r_hi + 1])),
        "d_tm": round(abs(round(res[f"PRIMER_LEFT_{i}_TM"], 1) - round(res[f"PRIMER_RIGHT_{i}_TM"], 1)), 1),
        "f_gc": round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 1), "r_gc": round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 1),
        "f_len": llen, "r_len": rlen,
        "product": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
        "penalty": round(res[f"PRIMER_PAIR_{i}_PENALTY"], 2),
        "f_dist": t.lo - f_hi - 1,
        "r_dist": r_lo - t.hi - 1,
        "chrom": m.chrom, "start": g1, "end": g2,
        "amplicon": amp,
        "f_snps": f_snps, "r_snps": r_snps,
        "f_repeat": f_rep, "r_repeat": r_rep, "repeat_names": sorted(set(f_rnames + r_rnames)),
        "ispcr": {"checked": False},
    }


# --------------------------------------------------------------------------
# Whole-gene driver
# --------------------------------------------------------------------------
def design_gene(m: GeneModel, p: Params, assembly: str) -> dict:
    amplicons, failures = [], []
    for t in make_targets(m, p):
        r = design_target(m, t, p)
        tg1, tg2 = sorted((m.to_genomic(t.lo), m.to_genomic(t.hi)))
        sd = [s for s in m.segdups if s["start"] <= tg2 and s["end"] >= tg1]
        rec = {"label": t.label, "exons": t.exons, "fragment": t.fragment,
               "target_start": tg1, "target_end": tg2, "target_len": t.hi - t.lo + 1,
               "relaxed": r["relaxed"], "tier": r["tier"], "candidates": r["pairs"],
               "segdups": sorted({(s["other"], s["identity"]) for s in sd}, key=lambda x: -x[1])}
        rec["segdups"] = [{"other": o, "identity": i} for o, i in rec["segdups"]]
        if r["pairs"]:
            amplicons.append(rec)
        else:
            rec["reason"] = r["explain"]
            failures.append(rec)
    g = {"gene": m.gene, "gene_id": m.gene_id, "transcript": m.transcript_version, "refseq": m.refseq,
         "selection": m.selection, "assembly": assembly, "chrom": m.chrom, "strand": m.strand,
         "exon_count": len(m.exons), "annotated": m.annotated,
         "exons": [{"n": e.number, "id": e.ensembl_id, "start": e.start, "end": e.end,
                    "cds_start": e.cds_start, "cds_end": e.cds_end, "coding": e.coding,
                    "len": e.end - e.start + 1} for e in m.exons],
         "amplicons": amplicons, "failures": failures}
    finalize_gene(g, p)
    return g


def finalize_gene(g: dict, p: Params):
    """Rank pairs (specific first), keep num_return, write status/warnings."""
    for a in g["amplicons"]:
        cands = sorted(a["candidates"], key=lambda b: (bool(b["ispcr"].get("off_target")), b["rank"]))
        a["pairs"] = [dict(b, rank=k) for k, b in enumerate(cands[:p.num_return], 1)]
        a["warnings"] = amplicon_warnings(g, a)
        a["status"] = "OK" if not a["warnings"] else "CHECK"


def amplicon_warnings(g: dict, a: dict) -> list[str]:
    b, w = a["pairs"][0], []
    for side in ("f", "r"):
        for s in b[f"{side}_snps"]:
            maf = f'MAF {s["maf"]:.2f}' if s["maf"] is not None else "common"
            w.append(f'{side.upper()} primer: SNP {s["rs"]} ({maf}) {s["pos3"]} nt from 3\' end')
        if b[f"{side}_repeat"]:
            w.append(f'{side.upper()} primer: {b[f"{side}_repeat"]} nt in repeat (5\' part)')
    if b["repeat_names"]:
        w.append("Repeats under primers: " + ", ".join(b["repeat_names"]))
    for s in a["segdups"][:3]:
        w.append(f'Target in segmental duplication: {s["identity"]}% identity with {s["other"]}')
    ip = b["ispcr"]
    if ip.get("checked"):
        if ip.get("off_target"):
            w.append(f'In-silico PCR: {len(ip["off_target"])} off-target product(s): ' + "; ".join(ip["off_target"][:3]))
        if not ip.get("on_target"):
            w.append("In-silico PCR: expected product not found")
    if a["relaxed"]:
        w.append("Relaxed Primer3 settings")
    for kind in ("snps", "repeats", "segdups"):
        if not g["annotated"].get(kind):
            w.append(f"{kind} not checked (UCSC annotation unavailable)")
    return w


# --------------------------------------------------------------------------
# Genome-wide specificity: UCSC isPcr
# --------------------------------------------------------------------------
def run_ispcr(results: list[dict], ispcr: str, genome: str, p: Params, log=print) -> bool:
    queries, where = [], {}
    for gi, g in enumerate(results):
        for ai, a in enumerate(g["amplicons"]):
            for pi, b in enumerate(a["candidates"]):
                pid = f"q{gi}_{ai}_{pi}"
                queries.append(f"{pid}\t{b['f_seq']}\t{b['r_seq']}")
                where[pid] = b
    if not queries:
        return True
    with tempfile.TemporaryDirectory() as td:
        q, out = Path(td) / "pairs.txt", Path(td) / "out.bed"
        q.write_text("\n".join(queries) + "\n")
        cmd = [ispcr, genome, str(q), str(out), "-out=bed", f"-maxSize={p.ispcr_max_size}",
               "-minPerfect=15", "-minGood=15"]
        log(f"Running in-silico PCR on {len(queries)} primer pairs ...")
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except (OSError, subprocess.CalledProcessError) as ex:
            log(f"warning: isPcr failed, specificity not checked ({getattr(ex, 'stderr', '') or ex})")
            return False
        apply_ispcr(where, out.read_text().splitlines(), results)
    return True


def apply_ispcr(where: dict, bed_lines: list[str], results: list[dict]):
    hits: dict[str, list] = {}
    for line in bed_lines:
        f = line.split("\t")
        if len(f) >= 4 and not line.startswith(("#", "track")):
            hits.setdefault(f[3], []).append((f[0], int(f[1]) + 1, int(f[2])))
    for pid, b in where.items():
        chrom = ucsc_chrom(b["chrom"])
        on, off = False, []
        for c, s, e in hits.get(pid, []):
            if c == chrom and s <= b["end"] and e >= b["start"]:
                on = True
            elif "_alt" not in c and "_fix" not in c:  # alt haplotypes duplicate the locus
                off.append(f"{c}:{s}-{e} ({e - s + 1} bp)")
        b["ispcr"] = {"checked": True, "on_target": on, "off_target": off}


# --------------------------------------------------------------------------
# Output helpers
# --------------------------------------------------------------------------
def ucsc_pcr_url(f: str, r: str, assembly: str) -> str:
    return "https://genome.ucsc.edu/cgi-bin/hgPcr?" + urlencode(
        {"db": UCSC_DB[assembly], "wp_target": "genome", "wp_f": f, "wp_r": r,
         "wp_size": 4000, "wp_perfect": 15, "wp_good": 15, "Submit": "submit"})


def snp_text(snps: list[dict]) -> str:
    return "; ".join(f'{s["rs"]} ({s["pos3"]} nt from 3\'' + (f', MAF {s["maf"]:.2f})' if s["maf"] is not None else ")")
                     for s in snps)


def ispcr_text(ip: dict) -> str:
    if not ip.get("checked"):
        return "not run"
    if ip.get("off_target"):
        return f'{1 + len(ip["off_target"]) if ip.get("on_target") else len(ip["off_target"])} products: ' + "; ".join(ip["off_target"])
    return "specific (1 product)" if ip.get("on_target") else "no product found"


def read_gene_list(path: str) -> list[tuple[str, str | None]]:
    out, seen = [], set()
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = re.split(r"[\s,;]+", line)
        if parts[0].upper() not in seen:
            seen.add(parts[0].upper())
            out.append((parts[0], parts[1] if len(parts) > 1 else None))
    return out


# --------------------------------------------------------------------------
# Excel output
# --------------------------------------------------------------------------
def write_excel(results: list[dict], errors: list[dict], p: Params, path: str, assembly: str, ispcr_ran: bool):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    font, mono = Font(name="Arial", size=10), Font(name="Consolas", size=10)
    bold = Font(name="Arial", size=10, bold=True, color="FFFFFF")
    link = Font(name="Arial", size=10, color="0563C1", underline="single")
    head_fill = PatternFill("solid", start_color="1F4E5A")
    band = PatternFill("solid", start_color="EEF4F5")
    status_fill = {"OK": PatternFill("solid", start_color="D5EFD9"),
                   "CHECK": PatternFill("solid", start_color="FDEBC8"),
                   "FAIL": PatternFill("solid", start_color="F7D3D0")}
    thin = Border(bottom=Side(style="thin", color="C9D6D9"))

    def header(ws, cols, widths, freeze="C2"):
        ws.append(cols)
        for i, w in enumerate(widths, 1):
            c = ws.cell(row=1, column=i)
            c.font, c.fill = bold, head_fill
            c.alignment = Alignment(vertical="center", wrap_text=True)
            ws.column_dimensions[get_column_letter(i)].width = w
        ws.row_dimensions[1].height = 30
        ws.freeze_panes = freeze

    def style_rows(ws, mono_cols=(), link_cols=(), status_col=None, wrap_cols=()):
        prev, shade = None, False
        for row in ws.iter_rows(min_row=2):
            if row[0].value != prev:
                shade, prev = not shade, row[0].value
            for c in row:
                c.font = mono if c.column in mono_cols else font
                c.border = thin
                if shade:
                    c.fill = band
                if c.column in wrap_cols:
                    c.alignment = Alignment(wrap_text=True, vertical="top")
            for col in link_cols:
                c = row[col - 1]
                if c.value:
                    c.hyperlink, c.value, c.font = c.value, "in-silico PCR", link
            if status_col and row[status_col - 1].value in status_fill:
                row[status_col - 1].fill = status_fill[row[status_col - 1].value]

    ws = wb.active
    ws.title = "Primers"
    cols = ["Gene", "Exon", "Status", "Fragment", "Transcript", "RefSeq (MANE)",
            "Forward name", "Forward primer (5'→3')", "F len", "F Tm PCR (°C)", "F Tm Na+ (°C)", "F GC %",
            "Reverse name", "Reverse primer (5'→3')", "R len", "R Tm PCR (°C)", "R Tm Na+ (°C)", "R GC %",
            "ΔTm F-R (°C)",
            "Product (bp)", "F 3' to exon (bp)", "Exon to R 3' (bp)",
            "Chr", f"Amplicon start ({assembly})", f"Amplicon end ({assembly})", "Target (bp)",
            "Design tier", "Common SNPs under F", "Common SNPs under R", "Repeat bases F/R",
            "Segmental duplication", "In-silico PCR", "Warnings", "UCSC check"]
    header(ws, cols, [10, 7, 8, 8, 19, 15, 16, 30, 6, 11, 11, 7, 16, 30, 6, 11, 11, 7, 11,
                      9, 9, 9, 5, 13, 13, 8,
                      11, 26, 26, 9, 30, 30, 60, 13])
    for g in results:
        for a in g["amplicons"]:
            b = a["pairs"][0]
            name = f'{g["gene"]}_E{a["label"]}'
            ws.append([g["gene"], a["label"], a["status"], a["fragment"], g["transcript"], g["refseq"],
                       name + "F", b["f_seq"], b["f_len"], b["f_tm"], b["f_tm_na"], b["f_gc"],
                       name + "R", b["r_seq"], b["r_len"], b["r_tm"], b["r_tm_na"], b["r_gc"],
                       b["d_tm"],
                       b["product"], b["f_dist"], b["r_dist"], b["chrom"], b["start"], b["end"], a["target_len"],
                       a["tier"] + (" (relaxed)" if a["relaxed"] else ""),
                       snp_text(b["f_snps"]), snp_text(b["r_snps"]),
                       f'{b["f_repeat"]}/{b["r_repeat"]}' if (b["f_repeat"] or b["r_repeat"]) else "",
                       "; ".join(f'{s["identity"]}% {s["other"]}' for s in a["segdups"][:3]),
                       ispcr_text(b["ispcr"]), "\n".join(a["warnings"]),
                       ucsc_pcr_url(b["f_seq"], b["r_seq"], assembly)])
        for f in g["failures"]:
            ws.append([g["gene"], f["label"], "FAIL", f["fragment"], g["transcript"], g["refseq"]]
                      + [""] * (len(cols) - 8) + ["No primer pair found. " + f.get("reason", ""), ""])
    style_rows(ws, mono_cols=(8, 14), link_cols=(len(cols),), status_col=3)
    ws.auto_filter.ref = ws.dimensions

    ws2 = wb.create_sheet("Alternatives")
    header(ws2, ["Gene", "Exon", "Rank", "Forward primer (5'→3')", "F Tm PCR (°C)", "F Tm Na+ (°C)",
                 "Reverse primer (5'→3')", "R Tm PCR (°C)", "R Tm Na+ (°C)", "ΔTm F-R (°C)",
                 "Product (bp)", "Penalty", "SNPs F", "SNPs R", "Repeat bases F/R",
                 "In-silico PCR", f"Amplicon start ({assembly})", f"Amplicon end ({assembly})"],
           [10, 7, 6, 30, 11, 11, 30, 11, 11, 11, 9, 8, 24, 24, 9, 28, 13, 13])
    for g in results:
        for a in g["amplicons"]:
            for b in a["pairs"]:
                ws2.append([g["gene"], a["label"], b["rank"], b["f_seq"], b["f_tm"], b["f_tm_na"],
                            b["r_seq"], b["r_tm"], b["r_tm_na"], b["d_tm"],
                            b["product"], b["penalty"], snp_text(b["f_snps"]), snp_text(b["r_snps"]),
                            f'{b["f_repeat"]}/{b["r_repeat"]}', ispcr_text(b["ispcr"]), b["start"], b["end"]])
    style_rows(ws2, mono_cols=(4, 7))
    ws2.auto_filter.ref = ws2.dimensions

    ws3 = wb.create_sheet("Amplicons")
    header(ws3, ["Gene", "Exon", "Product (bp)", "Amplicon sequence (exon UPPER, intron lower)"], [10, 7, 9, 120])
    for g in results:
        for a in g["amplicons"]:
            ws3.append([g["gene"], a["label"], a["pairs"][0]["product"], a["pairs"][0]["amplicon"]])
    style_rows(ws3, mono_cols=(4,))

    ws4 = wb.create_sheet("Genes")
    header(ws4, ["Gene", "Transcript", "RefSeq (MANE)", "Selection", "Chr", "Strand", "Exons (total)",
                 "Coding exons", "Amplicons", "OK", "To check", "Failed", "Common SNPs in region"],
           [10, 19, 15, 16, 5, 7, 9, 9, 9, 7, 8, 7, 10])
    for g in results:
        st = [a["status"] for a in g["amplicons"]]
        ws4.append([g["gene"], g["transcript"], g["refseq"], g["selection"], g["chrom"],
                    "+" if g["strand"] == 1 else "-", g["exon_count"], sum(e["coding"] for e in g["exons"]),
                    len(g["amplicons"]), st.count("OK"), st.count("CHECK"), len(g["failures"]),
                    g.get("n_snps", "") if g["annotated"].get("snps") else "not checked"])
    style_rows(ws4)

    ws5 = wb.create_sheet("Problems")
    header(ws5, ["Gene", "Exon", "Target (bp)", "Problem"], [10, 8, 9, 110])
    for g in results:
        for f in g["failures"]:
            ws5.append([g["gene"], f["label"], f["target_len"], "No primer pair found. " + f.get("reason", "")])
    for e in errors:
        ws5.append([e["gene"], "", "", e["error"]])
    style_rows(ws5)

    ws6 = wb.create_sheet("Parameters")
    header(ws6, ["Setting", "Value", "Meaning"], [26, 18, 95], freeze="A2")
    meaning = {
        "min_dist": "Minimal distance between primer 3' end and exon/intron boundary (bp)",
        "primer_region": "Width of the window, beyond min_dist, where primers may be placed (bp)",
        "product_min": "Minimal PCR product length (bp)",
        "product_max": "Maximal PCR product length (bp)",
        "max_target": "Targets larger than this are split into overlapping fragments, each with its own primer pair "
                      "(bp; derived from product_max - 2 x (min_dist + primer length + 20))",
        "overlap": "Overlap between fragments of a split exon (bp)",
        "combine": "Combine neighbouring exons into one amplicon when they fit in max_target",
        "target": "coding = coding part of each exon is the target; full = whole exon",
        "relax": "If no pair is found, retry with relaxed Tm/size/GC settings",
        "max_tm_diff": "Largest |Tm(F) - Tm(R)| accepted for a pair (degrees C)",
        "wt_tm_diff": "Primer3 penalty weight on that difference; its own default is 0, "
                      "so without this a matched pair is allowed but never preferred",
        "salt_mono": "Monovalent cations for the Tm calculation (mM)",
        "salt_div": "Divalent cations, Mg2+, for the Tm calculation (mM) - worth about 6 degrees "
                    "on a 20-mer, which is why a Na+-only calculator reads that much lower",
        "dntp": "Total dNTP concentration for the Tm calculation (mM)",
        "dna_conc": "Annealing oligo concentration for the Tm calculation (nM)",
        "snp_maf": f"SNPs from UCSC {SNP_TRACK} (1000 Genomes MAF >= 1%) at or above this MAF are avoided",
        "snp_3prime": "A SNP is never accepted within this many bases of a primer's 3' end",
        "fallback": "Allow tiers repeat-5p and snp-5p when no clean pair exists (always reported)",
        "ispcr_max_size": "In-silico PCR: largest product considered (bp)",
    }
    ws6.append(["Assembly", assembly, "Genome build used for coordinates (Ensembl REST)"])
    ws6.append(["Transcript choice", "MANE Select", "Falls back to Ensembl canonical; overridable in the gene list"])
    ws6.append(["Repeats", "masked", f"Ensembl soft-masking (RepeatMasker, TRF, Dust) + UCSC {RMSK_TRACK}"])
    ws6.append(["Homology", "flagged", f"UCSC {SEGDUP_TRACK} overlap; in-silico PCR " +
                ("run on all candidate pairs (isPcr)" if ispcr_ran else "NOT run - use the UCSC links")])
    for k, v in vars(p).items():
        ws6.append([k, str(v), meaning.get(k, "Primer3 setting")])
    for k, v in TIERS.items():
        ws6.append([f"Tier: {k}", "", v])
    ws6.append(["Status OK", "", "Clean tier, normal settings, no segmental duplication, no off-target product"])
    ws6.append(["Generated", time.strftime("%Y-%m-%d %H:%M"),
                f"exonprimer_batch.py v{__version__}, primer3-py {primer3.__version__}"])
    style_rows(ws6)

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)


def write_json(results, errors, p: Params, out_dir: str, assembly: str, ispcr_ran: bool):
    d = Path(out_dir)
    (d / "genes").mkdir(parents=True, exist_ok=True)
    index = []
    for g in results:
        slim = dict(g, amplicons=[{k: v for k, v in a.items() if k != "candidates"} for a in g["amplicons"]],
                    failures=[{k: v for k, v in f.items() if k != "candidates"} for f in g["failures"]])
        (d / "genes" / f'{g["gene"]}.json').write_text(json.dumps(slim, separators=(",", ":")))
        st =[a["status"] for a in g["amplicons"]]
        index.append({"gene": g["gene"], "transcript": g["transcript"], "refseq": g["refseq"],
                      "chrom": g["chrom"], "amplicons": len(st), "ok": st.count("OK"),
                      "check": st.count("CHECK"), "failures": len(g["failures"])})
    meta = {"assembly": assembly, "ucsc_db": UCSC_DB[assembly], "generated": time.strftime("%Y-%m-%d"),
            "version": __version__, "params": vars(p), "ispcr": ispcr_ran, "tiers": TIERS,
            "genes": index, "errors": errors}
    (d / "index.json").write_text(json.dumps(meta, indent=1))


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main(argv=None, ensembl=None, ucsc=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genes", help="text file with one gene symbol per line (optional transcript ID after it)")
    ap.add_argument("-o", "--out", default="exon_primers.xlsx", help="Excel output file")
    ap.add_argument("--json-dir", help="also write per-gene JSON for the website here")
    ap.add_argument("--assembly", choices=list(SERVERS), default="GRCh38")
    ap.add_argument("--cache", default=".cache", help="cache directory for web responses ('' to disable)")
    d = Params()
    g1 = ap.add_argument_group("ExonPrimer layout")
    g1.add_argument("--min-dist", type=int, default=d.min_dist)
    g1.add_argument("--primer-region", type=int, default=d.primer_region)
    g1.add_argument("--product", type=int, nargs=2, metavar=("MIN", "MAX"), default=[d.product_min, d.product_max],
                    help="PCR product length range in bp (default 200 500)")
    g1.add_argument("--max-target", type=int, default=0,
                    help="split targets above this size (default: derived from the product range)")
    g1.add_argument("--overlap", type=int, default=d.overlap)
    g1.add_argument("--no-combine", action="store_true", help="one amplicon per exon, never combine")
    g1.add_argument("--target", choices=["coding", "full"], default=d.target)
    g2 = ap.add_argument_group("Primer3")
    g2.add_argument("--tm", type=float, nargs=3, metavar=("MIN", "OPT", "MAX"), default=[d.tm_min, d.tm_opt, d.tm_max])
    g2.add_argument("--size", type=int, nargs=3, metavar=("MIN", "OPT", "MAX"), default=[d.size_min, d.size_opt, d.size_max])
    g2.add_argument("--gc", type=float, nargs=2, metavar=("MIN", "MAX"), default=[d.gc_min, d.gc_max])
    g2.add_argument("--gc-clamp", type=int, default=d.gc_clamp)
    g2.add_argument("--alternatives", type=int, default=d.num_return, help="pairs to keep per amplicon")
    g2.add_argument("--no-relax", action="store_true")
    g2.add_argument("--max-tm-diff", type=float, default=d.max_tm_diff,
                    help="largest |Tm(F) - Tm(R)| accepted for a pair (default 2)")
    g2.add_argument("--wt-tm-diff", type=float, default=d.wt_tm_diff,
                    help="Primer3 weight on the pair Tm difference (default 1; Primer3's own default is 0)")
    g2.add_argument("--salt", type=float, nargs=4, metavar=("MONO_MM", "DIV_MM", "DNTP_MM", "DNA_NM"),
                    default=[d.salt_mono, d.salt_div, d.dntp, d.dna_conc],
                    help="Tm conditions: monovalent mM, divalent mM, dNTP mM, oligo nM "
                         "(default 50 1.5 0.6 50 - PCR conditions)")
    g3 = ap.add_argument_group("SNPs, repeats, homology")
    g3.add_argument("--snp-maf", type=float, default=d.snp_maf, help="avoid SNPs with MAF >= this (default 0.01)")
    g3.add_argument("--snp-3prime", type=int, default=d.snp_3prime, help="no SNP allowed in this many 3' bases")
    g3.add_argument("--strict", action="store_true", help="only accept primers with no SNP and no repeat at all")
    g3.add_argument("--no-ucsc", action="store_true", help="skip UCSC SNP/repeat/segdup annotation")
    g3.add_argument("--ucsc-delay", type=float, default=1.0,
                    help="seconds between UCSC API requests (UCSC asks for <= 1 hit/s)")
    g3.add_argument("--no-mane-track", action="store_true",
                    help="never use the UCSC MANE track; ask Ensembl per gene (slow)")
    g3.add_argument("--seq-batch", type=int, default=50,
                    help="regions per Ensembl POST sequence request (max 50)")
    g3.add_argument("--ispcr", help="path to UCSC isPcr executable (genome-wide specificity check)")
    g3.add_argument("--genome-2bit", help="genome .2bit for isPcr (hg38.2bit / hg19.2bit)")
    g3.add_argument("--ispcr-max-size", type=int, default=d.ispcr_max_size)
    a = ap.parse_args(argv)

    try:
        p = Params(min_dist=a.min_dist, primer_region=a.primer_region, max_target=a.max_target,
                   overlap=a.overlap, product_min=a.product[0], product_max=a.product[1],
                   combine=not a.no_combine, target=a.target, tm_min=a.tm[0], tm_opt=a.tm[1], tm_max=a.tm[2],
                   size_min=a.size[0], size_opt=a.size[1], size_max=a.size[2], gc_min=a.gc[0], gc_max=a.gc[1],
                   gc_clamp=a.gc_clamp, num_return=a.alternatives, relax=not a.no_relax, snp_maf=a.snp_maf,
                   snp_3prime=a.snp_3prime, fallback=not a.strict, ispcr_max_size=a.ispcr_max_size,
                   max_tm_diff=a.max_tm_diff, wt_tm_diff=a.wt_tm_diff,
                   salt_mono=a.salt[0], salt_div=a.salt[1], dntp=a.salt[2], dna_conc=a.salt[3])
    except ValueError as ex:
        ap.error(str(ex))
    if bool(a.ispcr) != bool(a.genome_2bit):
        ap.error("--ispcr and --genome-2bit go together")

    def log(msg):
        print(msg, file=sys.stderr, flush=True)

    ens = ensembl or Ensembl(a.assembly, a.cache or None)
    if ucsc is None and not a.no_ucsc:
        ucsc = UCSC(a.assembly, a.cache or None, a.ucsc_delay)
    pad = p.min_dist + p.primer_region + 100
    genes = read_gene_list(a.genes)
    results, errors = [], []

    # ---- phase 1: exon structures ---------------------------------------
    # The UCSC MANE track covers almost every gene in one request per
    # chromosome; Ensembl is asked only for the rest - genes with an explicit
    # transcript, genes with no MANE transcript, and anything the track could
    # not provide. Ensembl remains the authority, it is just no longer asked
    # the same question thousands of times.
    mane = ManeTable(None if (a.no_ucsc or a.no_mane_track) else ucsc, log)
    models, from_mane = [], 0
    for i, (sym, tx) in enumerate(genes, 1):
        try:
            m = None if tx else mane.model(sym, pad)
            if m is not None:
                from_mane += 1
            else:
                m = build_model(ens, sym, tx, pad, with_seq=False)
            models.append((sym, m))
        except Exception as ex:
            errors.append({"gene": sym, "error": f"{type(ex).__name__}: {ex}"})
            log(f"[{i}/{len(genes)}] {sym:<10} ERROR {ex}")
    log(f"Structures: {len(models)}/{len(genes)} ({from_mane} from the MANE track, "
        f"{len(models) - from_mane} from Ensembl)")

    # ---- phase 2: sequences, in batches ---------------------------------
    log(f"Fetching soft-masked sequence for {len(models)} genes ...")
    fetch_sequences(ens, [m for _, m in models], a.seq_batch, log)
    short = [(sym, m) for sym, m in models if not m.seq]
    for sym, _ in short:
        errors.append({"gene": sym, "error": "LookupError: no sequence returned by Ensembl"})
        log(f"{sym:<10} ERROR no sequence returned")
    models = [(sym, m) for sym, m in models if m.seq]

    # ---- phase 3: annotate and design -----------------------------------
    for i, (sym, m) in enumerate(models, 1):
        try:
            annotate(m, None if a.no_ucsc else ucsc, p, log)
            g = design_gene(m, p, a.assembly)
            g["n_snps"] = len(m.snps)
            results.append(g)
            st = [x["status"] for x in g["amplicons"]]
            log(f"[{i}/{len(models)}] {g['gene']:<10} {g['transcript']:<20} {len(st)} amplicons "
                f"({st.count('OK')} OK, {st.count('CHECK')} to check), {len(g['failures'])} failed")
        except Exception as ex:
            errors.append({"gene": sym, "error": f"{type(ex).__name__}: {ex}"})
            log(f"[{i}/{len(models)}] {sym:<10} ERROR {ex}")

    ispcr_ran = False
    if a.ispcr:
        ispcr_ran = run_ispcr(results, a.ispcr, a.genome_2bit, p, log)
        if ispcr_ran:
            for g in results:
                finalize_gene(g, p)

    write_excel(results, errors, p, a.out, a.assembly, ispcr_ran)
    if a.json_dir:
        write_json(results, errors, p, a.json_dir, a.assembly, ispcr_ran)
    amps = [x for g in results for x in g["amplicons"]]
    log(f"Done: {len(results)} genes, {len(amps)} amplicons "
        f"({sum(x['status'] == 'OK' for x in amps)} OK), "
        f"{sum(len(g['failures']) for g in results)} failed targets, {len(errors)} gene errors -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
