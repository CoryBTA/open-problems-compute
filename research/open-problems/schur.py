"""Schur numbers via SAT: partition {1..n} into k sum-free classes (Erdos-Ramsey type problem).

S(k) is the largest n for which such a partition exists. Known: S(1..5) = 1, 4, 13, 44, 160;
S(6) >= 536 is open. verify_schur() is independent of the SAT encoding.
"""
from __future__ import annotations

import json
import sys
import threading
import time
from pathlib import Path

from pysat.solvers import Glucose4

RESULTS = Path(__file__).parent / "results"


def verify_schur(coloring: list[int], k: int) -> bool:
    """True when coloring[i-1] is the class of i, all in range(k), and no class has a+b=c (a may equal b)."""
    n = len(coloring)
    if any(not 0 <= c < k for c in coloring):
        return False
    for a in range(1, n + 1):
        for b in range(a, n + 1 - a):
            c = a + b
            if coloring[a - 1] == coloring[b - 1] == coloring[c - 1]:
                return False
    return True


def solve_schur(n: int, k: int, timeout: float | None = None):
    """Return a k-coloring of 1..n with no monochromatic a+b=c, None when unsatisfiable, 'timeout' if cut off.

    Symmetry breaking: colour classes are named by first appearance (colour c may only appear after colour c-1 has,
    which implies colour(1)=0). A palindromic constraint c(i)=c(n+1-i) was tried and is wrong: S(4)=44 becomes UNSAT.
    """
    var = lambda i, c: (i - 1) * k + c + 1  # noqa: E731
    s = Glucose4()
    for i in range(1, n + 1):
        s.add_clause([var(i, c) for c in range(k)])
    for c in range(k):
        for a in range(1, n + 1):
            for b in range(a, n + 1 - a):
                if a == b:
                    s.add_clause([-var(a, c), -var(2 * a, c)])
                else:
                    s.add_clause([-var(a, c), -var(b, c), -var(a + b, c)])
    s.add_clause([var(1, 0)])
    for i in range(1, n + 1):
        for c in range(1, k):
            # colour c at i needs colour c-1 somewhere before i
            s.add_clause([-var(i, c)] + [var(j, c - 1) for j in range(1, i)])
    timer = None
    if timeout:
        timer = threading.Timer(timeout, s.interrupt)
        timer.start()
    res = s.solve_limited(expect_interrupt=True)
    if timer:
        timer.cancel()
    if res is None:
        s.delete()
        return "timeout"
    if not res:
        s.delete()
        return None
    model = {x for x in s.get_model() if x > 0}
    s.delete()
    return [next(c for c in range(k) if var(i, c) in model) for i in range(1, n + 1)]


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    out = []
    for k, n, limit in [(2, 4, 60), (2, 5, 60), (3, 13, 60), (3, 14, 60), (4, 44, 60), (4, 45, 120), (5, 160, 1500)]:
        t = time.time()
        col = solve_schur(n, k, timeout=limit)
        status = "timeout" if col == "timeout" else ("sat" if col else "unsat")
        row = {"k": k, "n": n, "status": status, "verified": isinstance(col, list) and verify_schur(col, k),
               "seconds": round(time.time() - t, 1)}
        if isinstance(col, list) and k == 5:
            row["coloring"] = col
        out.append(row)
        print({**row, "coloring": "..."} if "coloring" in row else row, flush=True)
        (RESULTS / "schur.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
