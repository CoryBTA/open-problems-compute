"""Merge the per-shard summaries of a run6b sweep (Actions matrix) into one labelled result.

Usage: python merge_shards.py DIR EXPECTED_SHARDS N OUT.json
       python merge_shards.py --finish OUT.json N PREFIX RUN_ID   (OUT.json is created if absent; prefix and run are recorded)
The label is `negative (this class only)` only when every expected shard reported, none timed out and none found a
witness; any timeout or missing shard makes it `inconclusive`. A witness is never accepted here: run6b re-checks it and
writes results/s6b-N.json itself, and that file must be re-checked by research/open-problems/schur.py before any claim.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def merge(directory: Path, expected: int, n: int) -> dict:
    """Combine the shard summaries in `directory` for value n."""
    rows = [json.loads(p.read_text()) for p in sorted(directory.glob("**/shard-*.json"))]
    rows = [r for r in rows if r.get("n") == n]
    total = {"sat": 0, "unsat": 0, "timeout": 0}
    for r in rows:
        for key in total:
            total[key] += r["stats"][key]
    complete = len(rows) == expected
    if total["sat"]:
        label = "witness-found (verify independently before any claim)"
    elif complete and not total["timeout"]:
        label = "negative (this class only)"
    else:
        label = "inconclusive"
    return {"n": n, "shards": len(rows), "expected_shards": expected, "starts": sum(r["starts"] for r in rows),
            "stats": total, "label": label}


def finish(out: Path, n: int, prefix: int, run: str) -> None:
    """Complete the result file: an absent one (no shard reported) becomes an inconclusive row; prefix and run are recorded."""
    row = json.loads(out.read_text()) if out.exists() and out.stat().st_size else {"n": n, "label": "inconclusive", "note": "no shard reported"}
    out.write_text(json.dumps({**row, "prefix": prefix, "run": run}, indent=1))


def main() -> None:
    if sys.argv[1] == "--finish":
        finish(Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
        return
    directory, expected, n, out = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4])
    result = merge(directory, expected, n)
    out.write_text(json.dumps(result, indent=1))
    print(json.dumps(result))
    if result["label"] == "inconclusive":
        sys.exit(1)


if __name__ == "__main__":
    main()
