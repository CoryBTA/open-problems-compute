"""Non-symmetric S(5)=160 calibration: how far does structure short of the full symmetry get?

Variants: depth44 (the last set has no element below 44, no symmetry), plain (colour-ordering only).
Usage: python run5_nonsym.py VARIANT [--timeout S]
"""
from __future__ import annotations

import json
import sys
import time

from search import RESULTS, solve, verify


def main() -> None:
    variant = sys.argv[1]
    limit = float(sys.argv[sys.argv.index("--timeout") + 1]) if "--timeout" in sys.argv else 1500.0
    case = {"n": 160, "k": 5, "symmetric": False, "min_depth": 44 if variant == "depth44" else 0}
    t = time.time()
    w = solve(case, timeout=limit)
    status = "timeout" if w == "timeout" else ("sat" if w is not None else "unsat")
    row = {"variant": variant, "status": status, "verified": isinstance(w, list) and verify(w, 5), "seconds": round(time.time() - t, 1)}
    print(row, flush=True)
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / f"s5-nonsym-{variant}.json").write_text(json.dumps({**row, "coloring": w if isinstance(w, list) else None}))


if __name__ == "__main__":
    main()
