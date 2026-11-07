#!/usr/bin/env python3
"""Add the f_tm_ucsc / r_tm_ucsc fields to already-built site data.

The UCSC Tm is a pure function of the primer sequence, so the whole site does
not have to be redesigned to gain the column - every pair already in
site/data/genes/*.json can be filled in offline.

  python scripts/backfill_tm_ucsc.py site/data

Idempotent: pairs that already carry the field are left alone. New runs of
exonprimer_batch.py write it themselves, so this is only for data built
before the column existed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from exonprimer_batch import tm_ucsc  # noqa: E402


def main(argv: list[str]) -> int:
    root = Path(argv[1] if len(argv) > 1 else "site/data")
    files = sorted((root / "genes").glob("*.json"))
    if not files:
        print(f"no gene files under {root / 'genes'}", file=sys.stderr)
        return 1
    cache: dict[str, float] = {}
    changed = pairs = 0
    for f in files:
        g = json.loads(f.read_text())
        touched = False
        for a in g.get("amplicons", []):
            for b in a.get("pairs", []):
                pairs += 1
                if "f_tm_ucsc" in b and "r_tm_ucsc" in b:
                    continue
                for side in ("f", "r"):
                    seq = b[f"{side}_seq"]
                    if seq not in cache:
                        cache[seq] = tm_ucsc(seq)
                    b[f"{side}_tm_ucsc"] = cache[seq]
                touched = True
        if touched:
            f.write_text(json.dumps(g, separators=(",", ":")))
            changed += 1
    print(f"{changed}/{len(files)} gene files updated, {pairs} pairs seen, "
          f"{len(cache)} distinct primers computed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
