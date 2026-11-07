"""Offline tests: synthetic genome, fake Ensembl / UCSC responses, fake isPcr.

Checks
  * primers match the genome at the reported coordinates on both strands,
    keep min_dist from the exon, and all coding bases are covered;
  * exon combining and splitting;
  * common SNPs are avoided, or tolerated only outside the 3'-terminal bases
    and reported;
  * repeats are avoided, or tolerated only in the 5' part and reported;
  * segmental duplications are flagged;
  * in-silico PCR off-target products demote a pair and are reported.

Run:  python tests/test_offline.py   (or python -m pytest tests/)
"""
import json
import os
import random
import stat
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import exonprimer_batch as epb  # noqa: E402

rc = epb.rc
random.seed(7)
CHROM = "".join(random.choice("ACGT") for _ in range(60000))
CHROM = CHROM[:9000] + CHROM[9000:9400].lower() + CHROM[9400:]  # a soft-masked repeat

GENES = {
    "PLUSGENE": dict(strand=1, exons=[(5000, 5300), (8000, 8120), (8300, 8420),
                                      (12000, 13400), (20000, 20600)], tl=(5150, 20250)),
    "MINUSGENE": dict(strand=-1, exons=[(30000, 30500), (33000, 33090),
                                        (36000, 37200), (40000, 40250)], tl=(30300, 40100)),
}


class FakeEnsembl:
    def __init__(self, *a, **k):
        pass

    def gene(self, sym):
        if sym not in GENES:
            raise LookupError(f"No valid lookup found for symbol {sym}")
        g = GENES[sym]
        return {"id": f"ENSG_{sym}", "display_name": sym, "Transcript": [{
            "id": f"ENST_{sym}", "version": 1, "strand": g["strand"], "seq_region_name": "1",
            "is_canonical": 1, "MANE": [{"type": "MANE_Select", "refseq_match": f"NM_{sym}.1"}],
            "Exon": [{"id": f"ENSE{i}", "start": s, "end": e} for i, (s, e) in enumerate(g["exons"])],
            "Translation": {"start": g["tl"][0], "end": g["tl"][1], "length": 999}}]}

    def sequence(self, chrom, start, end, strand):
        self.gets = getattr(self, "gets", 0) + 1
        s = CHROM[start - 1:end]
        return s if strand == 1 else rc(s)

    def sequences(self, regions):
        """The POST batch, keyed by the query string as Ensembl keys it."""
        self.posts = getattr(self, "posts", 0) + 1
        out = {}
        for q in regions:
            chrom, rest = q.split(":")
            span, strand = rest.rsplit(":", 1)
            a, b = (int(x) for x in span.split(".."))
            out[q] = self.sequence(chrom, a, b, int(strand))
        self.gets -= len(regions)          # charge the batch, not the pieces
        return out


class FakeUCSC:
    def __init__(self, snps=(), repeats=(), segdups=()):
        self._s, self._r, self._d = list(snps), list(repeats), list(segdups)

    def snps(self, chrom, start, end, min_maf):
        return [s for s in self._s if s["end"] >= start and s["start"] <= end and (s["maf"] or 1) >= min_maf]

    def repeats(self, chrom, start, end):
        return self._r

    def segdups(self, chrom, start, end):
        return self._d

    def mane(self, chrom):
        """The UCSC MANE track, in the shape the real API serves it:
        0-based half-open blocks, CDS in thickStart/thickEnd."""
        if chrom != "chr1":
            return []
        rows = []
        for sym, g in GENES.items():
            ex = sorted(g["exons"])
            st = ex[0][0] - 1
            rows.append({
                "chrom": "chr1", "chromStart": st, "chromEnd": ex[-1][1],
                "name": f"ENST_{sym}.1", "geneName": f"ENSG_{sym}", "geneName2": sym,
                "ncbiId": f"NM_{sym}.1", "maneStat": "MANE Select",
                "strand": "+" if g["strand"] == 1 else "-",
                "thickStart": g["tl"][0] - 1, "thickEnd": g["tl"][1],
                "blockSizes": ",".join(str(e - s + 1) for s, e in ex) + ",",
                "chromStarts": ",".join(str(s - 1 - st) for s, _ in ex) + ",",
            })
        return rows


def snp(pos, maf=0.2, rs=None):
    return {"rs": rs or f"rs{pos}", "start": pos, "end": pos, "maf": maf, "class": "snv"}


def run(sym, p=None, ucsc=None):
    p = p or epb.Params()
    m = epb.build_model(FakeEnsembl(), sym, None, p.min_dist + p.primer_region + 100)
    epb.annotate(m, ucsc if ucsc is not None else FakeUCSC(), p, log=lambda *_: None)
    return m, epb.design_gene(m, p, "GRCh38")


def amp(g, label):
    return next(a for a in g["amplicons"] if a["label"] == label)


def genomic(a, b):
    return CHROM[a - 1:b].upper()


def check_gene(sym, p=None, ucsc=None):
    p = p or epb.Params()
    m, g = run(sym, p, ucsc)
    assert g["refseq"] == f"NM_{sym}.1" and g["selection"] == "MANE Select"
    assert not g["failures"], g["failures"]
    covered = set()
    for a in g["amplicons"]:
        for b in a["candidates"]:
            seq = genomic(b["start"], b["end"])
            if m.strand == -1:
                seq = rc(seq)
            assert len(seq) == b["product"]
            assert p.product_min <= b["product"] <= p.product_max, (sym, a["label"], b["product"])
            assert seq.startswith(b["f_seq"]) and seq.endswith(rc(b["r_seq"])), (sym, a["label"])
            assert b["amplicon"].upper() == seq
            assert p.min_dist <= b["f_dist"] <= p.min_dist + p.primer_region
            assert b["r_dist"] >= p.min_dist
            assert b["start"] <= a["target_start"] and b["end"] >= a["target_end"]
        covered.update(range(a["target_start"], a["target_end"] + 1))
    for e in g["exons"]:
        if e["coding"]:
            assert set(range(e["cds_start"], e["cds_end"] + 1)) <= covered, e
    return m, g


# ---------------------------------------------------------------- layout
def test_plus_strand():
    _, g = check_gene("PLUSGENE")
    labels = [a["label"] for a in g["amplicons"]]
    # exons 2-3 span 421 bp: too long for one <=500-bp product -> one pair each
    assert labels[0] == "1" and "2" in labels and "3" in labels and "2-3" not in labels
    # 1401-bp exon -> 6 overlapping fragments, each with its own primer pair
    assert [l for l in labels if l.startswith("4")] == ["4a", "4b", "4c", "4d", "4e", "4f"]
    assert all(a["status"] == "OK" and a["tier"] == "clean" for a in g["amplicons"])


def test_minus_strand():
    _, g = check_gene("MINUSGENE")
    assert g["exons"][0]["start"] == 40000
    labels = [a["label"] for a in g["amplicons"]]
    assert labels[0] == "1" and labels[-1] == "4"


def test_no_combine_and_full_exon():
    _, g = check_gene("PLUSGENE", epb.Params(combine=False, target="full"))
    labels = [a["label"] for a in g["amplicons"]]
    assert "2" in labels and "3" in labels and "2-3" not in labels


def test_default_product_range_and_derived_target():
    p = epb.Params()
    assert (p.product_min, p.product_max, p.max_target) == (200, 500, 320)
    # small exon (exon 2 alone, 121 bp) still gives >= 200 bp
    _, g = check_gene("PLUSGENE", epb.Params(combine=False))
    assert amp(g, "2")["pairs"][0]["product"] >= 200


def test_custom_product_range():
    p = epb.Params(product_min=300, product_max=700)
    assert p.max_target == 520
    _, g = check_gene("PLUSGENE", p)
    assert [l for l in (a["label"] for a in g["amplicons"]) if l.startswith("4")] == ["4a", "4b", "4c"]
    assert "2-3" in [a["label"] for a in g["amplicons"]]  # now fits one product


def test_impossible_settings_rejected():
    for kw in (dict(product_max=200, max_target=300), dict(product_min=600, product_max=500)):
        try:
            epb.Params(**kw)
        except ValueError:
            continue
        raise AssertionError(kw)


# ---------------------------------------------------------------- SNPs
def test_snp_under_best_primer_is_avoided():
    for sym in GENES:
        _, g0 = run(sym)
        a0 = g0["amplicons"][-1]
        b0 = a0["pairs"][0]
        # put common SNPs under both best primers (3' end of F, middle of R)
        m0, _ = run(sym)
        f_hi = m0.to_genomic(m0.to_local(b0["start"] if m0.strand == 1 else b0["end"]) + b0["f_len"] - 1)
        r_mid = (b0["end"] - b0["r_len"] // 2) if m0.strand == 1 else (b0["start"] + b0["r_len"] // 2)
        uc = FakeUCSC(snps=[snp(f_hi), snp(r_mid)])
        _, g = check_gene(sym, ucsc=uc)
        a = amp(g, a0["label"])
        b = a["pairs"][0]
        assert a["tier"] == "clean" and a["status"] == "OK"
        assert not b["f_snps"] and not b["r_snps"]
        assert (b["f_seq"], b["r_seq"]) != (b0["f_seq"], b0["r_seq"])


def test_dense_snps_fall_back_to_5prime_only():
    m0, g0 = run("PLUSGENE")
    a0 = amp(g0, "5")
    # SNP every 8 bp across the whole left primer window of exon 5 (+ strand)
    w_start = a0["target_start"] - 50 - 150
    w_end = a0["target_start"] - 50
    uc = FakeUCSC(snps=[snp(x, 0.05) for x in range(w_start, w_end, 8)])
    _, g = check_gene("PLUSGENE", ucsc=uc)
    a = amp(g, "5")
    assert a["tier"] == "snp-5p" and a["status"] == "CHECK"
    for b in a["candidates"]:
        assert b["f_snps"] and all(s["pos3"] > 5 for s in b["f_snps"] + b["r_snps"])
    assert any("SNP rs" in w for w in a["warnings"])


def test_strict_mode_refuses_snps():
    m0, g0 = run("PLUSGENE")
    a0 = amp(g0, "5")
    w_start, w_end = a0["target_start"] - 200, a0["target_start"] - 50
    uc = FakeUCSC(snps=[snp(x, 0.05) for x in range(w_start, w_end, 8)])
    _, g = run("PLUSGENE", epb.Params(fallback=False), uc)
    assert any(f["label"] == "5" for f in g["failures"])


def test_maf_threshold():
    m0, g0 = run("PLUSGENE")
    a0 = amp(g0, "5")
    w_start, w_end = a0["target_start"] - 200, a0["target_start"] - 50
    rare = FakeUCSC(snps=[snp(x, 0.02) for x in range(w_start, w_end, 8)])
    _, g = run("PLUSGENE", epb.Params(snp_maf=0.05), rare)  # 2% SNPs ignored at 5% cut-off
    assert amp(g, "5")["tier"] == "clean"


# ---------------------------------------------------------------- repeats
def test_repeat_partly_covering_window_is_avoided():
    m0, g0 = run("PLUSGENE")
    a0 = amp(g0, "5")
    b0 = a0["pairs"][0]
    rep = {"start": b0["start"] - 5, "end": b0["start"] + 40, "name": "AluSx", "class": "SINE"}
    _, g = check_gene("PLUSGENE", ucsc=FakeUCSC(repeats=[rep]))
    a = amp(g, "5")
    assert a["tier"] == "clean" and a["pairs"][0]["f_repeat"] == 0


def test_repeat_only_5prime_allowed_as_fallback():
    m0, g0 = run("PLUSGENE")
    a0 = amp(g0, "5")
    # repeat over the left window except its last 12 bp: a primer can only
    # sit with its 5' part in the repeat and its 3' end in unique sequence
    w_end = a0["target_start"] - 50 - 1
    rep = {"start": a0["target_start"] - 400, "end": w_end - 12, "name": "L1MC", "class": "LINE"}
    m, g = check_gene("PLUSGENE", ucsc=FakeUCSC(repeats=[rep]))
    a = amp(g, "5")
    assert a["tier"] == "repeat-5p" and a["status"] == "CHECK"
    for b in a["candidates"]:
        f_end = b["start"] + b["f_len"] - 1
        assert b["f_repeat"] > 0 and f_end > rep["end"]       # 3' end outside the repeat
    assert any("L1MC" in w for w in a["warnings"])


# ---------------------------------------------------------------- homology
def test_segdup_flag():
    sd = {"start": 19000, "end": 21000, "other": "chr9:1000-3000", "identity": 98.7}
    _, g = run("PLUSGENE", ucsc=FakeUCSC(segdups=[sd]))
    a = amp(g, "5")
    assert a["status"] == "CHECK" and any("98.7%" in w for w in a["warnings"])
    assert amp(g, "1")["status"] == "OK"


def test_ispcr_off_target_demotes_pair():
    p = epb.Params()
    _, g = run("PLUSGENE", p)
    results = [g]
    a = amp(g, "5")
    first = a["candidates"][0]
    # fake isPcr: every candidate gives its on-target product; the best
    # candidate of exon 5 also amplifies chr9
    lines = {}
    for ai, x in enumerate(g["amplicons"]):
        for pi, b in enumerate(x["candidates"]):
            pid = f"q0_{ai}_{pi}"
            lines[pid] = [f"chr1\t{b['start'] - 1}\t{b['end']}\t{pid}\t1000\t+"]
            if b is first:
                lines[pid] += [f"chr9\t5000\t5600\t{pid}\t950\t+", f"chr1_KI270706v1_alt\t1\t500\t{pid}\t900\t+"]
    with tempfile.TemporaryDirectory() as td:
        mp = Path(td) / "map.json"
        mp.write_text(json.dumps(lines))
        fake = Path(td) / "isPcr"
        fake.write_text(f"""#!{sys.executable}
import json, sys
m = json.load(open({str(mp)!r}))
q, out = sys.argv[2], sys.argv[3]
with open(out, "w") as fh:
    for line in open(q):
        for l in m.get(line.split()[0], []):
            fh.write(l + "\\n")
""")
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
        assert epb.run_ispcr(results, str(fake), "genome.2bit", p, log=lambda *_: None)
    epb.finalize_gene(g, p)
    a = amp(g, "5")
    assert a["pairs"][0]["f_seq"] != first["f_seq"] or a["pairs"][0]["r_seq"] != first["r_seq"]
    assert a["pairs"][0]["ispcr"] == {"checked": True, "on_target": True, "off_target": []}
    demoted = next(b for b in a["candidates"] if b is first)
    assert demoted["ispcr"]["off_target"] == ["chr9:5001-5600 (600 bp)"]  # alt contig ignored
    assert a["status"] == "OK"
    assert epb.ispcr_text(a["pairs"][0]["ispcr"]) == "specific (1 product)"


def test_ispcr_failure_is_reported_not_fatal():
    p = epb.Params()
    _, g = run("PLUSGENE", p)
    assert epb.run_ispcr([g], "/nonexistent/isPcr", "x.2bit", p, log=lambda *_: None) is False


# ---------------------------------------------------------------- outputs
def test_outputs(out=Path("test_output")):
    out.mkdir(exist_ok=True)
    p = epb.Params()
    uc = FakeUCSC(snps=[snp(x, 0.05) for x in range(19800, 19950, 8)],
                  segdups=[{"start": 36000, "end": 37300, "other": "chr5:10-1310", "identity": 99.1}])
    res = [run(s, p, uc)[1] for s in GENES]
    errs = [{"gene": "NOTAGENE", "error": "LookupError: No valid lookup found for symbol NOTAGENE"}]
    epb.write_excel(res, errs, p, str(out / "test_primers.xlsx"), "GRCh38", False)
    epb.write_json(res, errs, p, str(out / "data"), "GRCh38", False)
    from openpyxl import load_workbook
    wb = load_workbook(out / "test_primers.xlsx")
    assert wb.sheetnames == ["Primers", "Alternatives", "Amplicons", "Genes", "Problems", "Parameters"]
    assert wb["Primers"].max_row - 1 == sum(len(g["amplicons"]) + len(g["failures"]) for g in res)
    j = json.loads((out / "data" / "genes" / "PLUSGENE.json").read_text())
    assert "candidates" not in j["amplicons"][0]


def test_cli_end_to_end():
    with tempfile.TemporaryDirectory() as td:
        gl = Path(td) / "genes.txt"
        gl.write_text("PLUSGENE\nMINUSGENE\nNOTAGENE\n# comment\n")
        rc_ = epb.main([str(gl), "-o", f"{td}/site/x.xlsx", "--json-dir", f"{td}/site/data", "--cache", ""],
                       ensembl=FakeEnsembl(), ucsc=FakeUCSC())
        assert rc_ == 0 and Path(f"{td}/site/x.xlsx").exists()
        meta = json.loads(Path(f"{td}/site/data/index.json").read_text())
        assert [x["gene"] for x in meta["genes"]] == ["PLUSGENE", "MINUSGENE"]
        assert meta["errors"][0]["gene"] == "NOTAGENE"


def test_mane_track_matches_the_ensembl_model():
    """The fast path must produce the same model as asking Ensembl per gene.

    This is the whole safety argument for the MANE track: 24 requests instead
    of thousands, but only if the structures are identical. Against the real
    data they agreed for 4,449 of 4,471 genes; here the equivalence is pinned
    for both strands.
    """
    p = epb.Params()
    pad = p.min_dist + p.primer_region + 100
    mane = epb.ManeTable(FakeUCSC(), log=lambda *_: None)
    assert set(mane.by_symbol) == {"PLUSGENE", "MINUSGENE"}
    for sym in GENES:
        a = mane.model(sym, pad)
        b = epb.build_model(FakeEnsembl(), sym, None, pad, with_seq=False)
        assert a is not None, sym
        assert a.strand == b.strand and a.chrom == b.chrom, sym
        assert a.transcript_version == b.transcript_version and a.refseq == b.refseq, sym
        assert a.region_start == b.region_start and a.region_end == b.region_end, sym
        assert [(e.number, e.start, e.end, e.cds_start, e.cds_end) for e in a.exons] == \
               [(e.number, e.start, e.end, e.cds_start, e.cds_end) for e in b.exons], sym
    assert mane.model("NOTAGENE", pad) is None


def test_unknown_gene_falls_back_to_ensembl():
    """A gene absent from the track is asked of Ensembl, not skipped."""
    mane = epb.ManeTable(FakeUCSC(), log=lambda *_: None)
    mane.by_symbol.pop("PLUSGENE")
    with tempfile.TemporaryDirectory() as td:
        gl = Path(td) / "genes.txt"
        gl.write_text("PLUSGENE\nMINUSGENE\n")
        ens = FakeEnsembl()
        rc_ = epb.main([str(gl), "-o", f"{td}/x.xlsx", "--json-dir", f"{td}/data", "--cache", ""],
                       ensembl=ens, ucsc=FakeUCSC())
        meta = json.loads(Path(f"{td}/data/index.json").read_text())
    assert rc_ == 0 and [x["gene"] for x in meta["genes"]] == ["PLUSGENE", "MINUSGENE"]


def test_sequences_are_fetched_in_one_batch():
    """Many genes, one POST - and matched by query string, not by position."""
    p = epb.Params()
    pad = p.min_dist + p.primer_region + 100
    ens = FakeEnsembl()
    mane = epb.ManeTable(FakeUCSC(), log=lambda *_: None)
    models = [mane.model(sym, pad) for sym in GENES]
    epb.fetch_sequences(ens, models, batch=50, log=lambda *_: None)
    assert getattr(ens, "posts", 0) == 1, "one POST should cover both genes"
    for m in models:
        assert len(m.seq) == m.region_end - m.region_start + 1
        direct = ens.sequence(m.chrom, m.region_start, m.region_end, m.strand)
        assert m.seq == direct, f"{m.gene}: batch and single-region sequence differ"


def test_batch_failure_falls_back_to_single_requests():
    class Broken(FakeEnsembl):
        def sequences(self, regions):
            raise RuntimeError("batch endpoint down")
    p = epb.Params()
    pad = p.min_dist + p.primer_region + 100
    mane = epb.ManeTable(FakeUCSC(), log=lambda *_: None)
    models = [mane.model(sym, pad) for sym in GENES]
    epb.fetch_sequences(Broken(), models, batch=50, log=lambda *_: None)
    assert all(m.seq for m in models), "a dead batch endpoint must not lose the sequences"


def test_both_tm_values_are_reported():
    _, g = run("PLUSGENE")
    for a in g["amplicons"]:
        for b in a["pairs"]:
            assert b["f_tm_na"] < b["f_tm"] and b["r_tm_na"] < b["r_tm"], \
                "Na+-only Tm must read lower than the PCR-condition Tm"
            assert 4 < b["f_tm"] - b["f_tm_na"] < 8, \
                f"the Mg2+ contribution should be about 6 C, got {b['f_tm'] - b['f_tm_na']:.1f}"
            assert b["d_tm"] == round(abs(b["f_tm"] - b["r_tm"]), 1)


def test_pairs_are_matched_within_the_tm_ceiling():
    p = epb.Params()
    assert p.max_tm_diff == 2.0 and p.wt_tm_diff > 0
    _, g = run("PLUSGENE", p)
    for a in g["amplicons"]:
        if a["relaxed"]:
            continue
        for b in a["pairs"]:
            assert b["d_tm"] <= p.max_tm_diff + 1e-9, f'{a["label"]} rank {b["rank"]}: dTm {b["d_tm"]}'


def test_salt_conditions_reach_primer3_and_change_the_tm():
    hot = epb.Params()
    cold = epb.Params(salt_div=0.0, dntp=0.0)
    assert hot.primer3_global()["PRIMER_SALT_DIVALENT"] == 1.5
    assert cold.primer3_global()["PRIMER_SALT_DIVALENT"] == 0.0
    assert hot.primer3_global()["PRIMER_PAIR_WT_DIFF_TM"] == hot.wt_tm_diff
    # Designing under Na+-only conditions picks different primers - Primer3
    # optimises towards 60 C under whatever conditions it is given - so the
    # two runs cannot be compared primer by primer. What must hold is that
    # when the design conditions ARE Na+-only, the reported Tm and the
    # Na+ column coincide, and that they diverge under PCR conditions.
    _, a = run("PLUSGENE", hot)
    _, b = run("PLUSGENE", cold)
    for x in amp(b, "1")["pairs"]:
        assert x["f_tm"] == x["f_tm_na"] and x["r_tm"] == x["r_tm_na"], \
            f'Na+-only design should report one Tm, got {x["f_tm"]} vs {x["f_tm_na"]}'
    for x in amp(a, "1")["pairs"]:
        assert x["f_tm"] > x["f_tm_na"], "PCR conditions must read higher than Na+ alone"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
