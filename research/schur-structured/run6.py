"""S(6) gate and frontier: symmetric 6-set partitions with depth >= 161 (Fredricksen-Sweet recipe).

Usage: python run6.py N [N ...] [--timeout S]. Known: 536 exists (FS 2000), also 499..531 at the listed depths.
Anything above 536 is an attempt at a new record and is not claimed unless verify() passes AND the witness is re-checked.
"""
from __future__ import annotations

import json
import sys
import time

from search import RESULTS, is_symmetric, solve, verify

args = sys.argv[1:]
limit = float(args[args.index("--timeout") + 1]) if "--timeout" in args else 600.0
ns = [int(a) for a in args if a.isdigit()]
RESULTS.mkdir(exist_ok=True)
out = RESULTS / "s6.json"
rows = json.loads(out.read_text()) if out.exists() else []
for n in ns:
    t = time.time()
    case = {"n": n, "k": 6, "symmetric": True, "min_depth": 161}
    w = solve(case, timeout=limit)
    status = "timeout" if w == "timeout" else ("sat" if w is not None else "unsat")
    ok = isinstance(w, list) and verify(w, 6) and is_symmetric(w)
    row = {"n": n, "status": status, "verified": ok, "seconds": round(time.time() - t, 1)}
    if ok:
        row["coloring"] = w
    rows.append(row)
    print({k: v for k, v in row.items() if k != "coloring"}, flush=True)
    out.write_text(json.dumps(rows))
