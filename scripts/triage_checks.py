#!/usr/bin/env python3
"""Rank the CHECK amplicons by how much they actually matter.

Every amplicon that is not OK is labelled CHECK, but the reasons carry very
different risk: an off-target PCR product can amplify a paralog, while a
repeat touching a primer's 5' end is usually harmless. Reviewing thousands of
CHECK rows in sheet order means spending most of the attention on the benign
majority.

This groups them by reason, and singles out the amplicons where *every*
alternative pair has the same problem - those are the ones the tool cannot fix
by choosing a different pair, so they are the real redesign queue.

  python scripts/triage_checks.py site/data            # a local build
  python scripts/triage_checks.py <base-url>/data      # the published site
"""
from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

# Worst first. An amplicon is reported under its most serious reason.
REASONS = [
    ("off-target",  "a second PCR product elsewhere in the genome"),
    ("segdup",      "target inside a segmental duplication"),
    ("no-product",  "in-silico PCR found no product at all"),
    ("snp-3prime",  "a common SNP under a primer"),
    ("relaxed",     "designed only after relaxing Tm/size/GC"),
    ("repeat-5p",   "a repeat touching a primer's 5' part"),
    ("other",       "flagged, reason not recognised"),
]


def load(src: str, rel: str):
    """Read one JSON file, from a directory or over HTTP.

    Thousands of sequential requests will draw 503s from a static host, so
    remote reads back off and retry. Running against a local directory - a
    build, or the unpacked data-cache release asset - avoids this entirely.
    """
    if not src.startswith(("http://", "https://")):
        return json.loads((Path(src) / rel).read_text())
    url = f"{src.rstrip('/')}/{rel}"
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504):
                raise
            wait = float(e.headers.get("Retry-After", 2 ** attempt))
        except urllib.error.URLError:
            wait = 2 ** attempt
        if attempt == 0:
            print(f"  {url}: retrying", file=sys.stderr)
        time.sleep(wait)
    raise RuntimeError(f"gave up on {url}")


def classify(pair: dict, amp: dict) -> str:
    """The most serious thing wrong with one candidate pair."""
    ip = pair.get("ispcr") or {}
    if ip.get("off_target"):
        return "off-target"
    if amp.get("segdups"):
        return "segdup"
    if ip.get("checked") and not ip.get("on_target"):
        return "no-product"
    if pair.get("f_snps") or pair.get("r_snps"):
        return "snp-3prime"
    if amp.get("relaxed"):
        return "relaxed"
    if pair.get("f_repeat") or pair.get("r_repeat"):
        return "repeat-5p"
    return "other"


def rank(reason: str) -> int:
    return next(i for i, (r, _) in enumerate(REASONS) if r == reason)


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "site/data"
    remote = src.startswith(("http://", "https://"))
    index = load(src, "index.json")
    genes = index["genes"]
    print(f"{index['assembly']} · built {index['generated']} · {len(genes)} genes\n")

    counts, stuck, totals = Counter(), [], Counter()
    for i, g in enumerate(genes, 1):
        if g.get("check", 0) == 0 and g.get("failures", 0) == 0:
            totals["amplicons"] += g.get("amplicons", 0)
            totals["ok"] += g.get("ok", 0)
            continue
        d = load(src, f"genes/{g['gene']}.json")
        if remote:
            time.sleep(0.05)
        for a in d["amplicons"]:
            totals["amplicons"] += 1
            if a.get("status") == "OK":
                totals["ok"] += 1
                continue
            pairs = a.get("pairs") or []
            if not pairs:
                continue
            best = classify(pairs[0], a)
            counts[best] += 1
            # Can a listed alternative avoid the problem entirely?
            alts = [classify(p, a) for p in pairs]
            if all(rank(x) <= rank(best) for x in alts) and rank(best) <= rank("no-product"):
                stuck.append((best, g["gene"], a["label"], pairs[0], a))
        if i % 500 == 0:
            print(f"  ...{i}/{len(genes)} genes", file=sys.stderr)

    print(f"{totals['amplicons']:,} amplicons · {totals['ok']:,} OK "
          f"({100*totals['ok']/max(1,totals['amplicons']):.1f}%) · "
          f"{sum(counts.values()):,} to check\n")
    print("By most serious reason, worst first:")
    for reason, what in REASONS:
        if counts[reason]:
            print(f"  {counts[reason]:6,}  {reason:<12} {what}")

    print(f"\nNo alternative pair avoids the problem - the redesign queue ({len(stuck)}):")
    if not stuck:
        print("  none")
    for reason, gene, label, b, a in sorted(stuck, key=lambda x: (rank(x[0]), x[1])):
        extra = ""
        if reason == "off-target":
            off = (b.get("ispcr") or {}).get("off_target") or []
            same = [o for o in off if f"({b['product']} bp)" in o]
            extra = f" -> {len(off)} off-target" + (
                f", {len(same)} the same size as the real product" if same else "")
        elif reason == "segdup":
            sd = (a.get("segdups") or [{}])[0]
            extra = f" -> {sd.get('identity','?')}% identical to {sd.get('other','?')}"
        print(f"  {reason:<12} {gene}_E{label:<6} chr{b['chrom']}:{b['start']}-{b['end']}{extra}")


if __name__ == "__main__":
    main()
