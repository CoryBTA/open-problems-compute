"""schur-structured: Schur partitions in the symmetric class of Fredricksen and Sweet (EJC 7, 2000, #R32).

A partition of [1, n] into k sum-free sets is symmetric when i and (n+1)-i share a set, except that when 3 divides n+1
the points (n+1)/3 and 2(n+1)/3 must be in different sets. Depth d(P) is the largest set minimum: depth >= m means one set
has no element below m. Calibration: the class must be satisfiable where the paper says it is, and the verifier must agree
with `research/open-problems/schur.py` (independent code). Run long cases under a hard time limit
"""
from __future__ import annotations

import json
import time
from pathlib import Path

RESULTS = Path(__file__).parent / "results"


def verify(coloring: list[int], k: int) -> bool:
    """coloring[i-1] is the set of i; no set holds a, b (a == b allowed) and a+b."""
    n = len(coloring)
    if any(not 0 <= c < k for c in coloring):
        return False
    return not any(coloring[a - 1] == coloring[b - 1] == coloring[a + b - 1] for a in range(1, n + 1) for b in range(a, n + 1 - a))


def is_symmetric(coloring: list[int]) -> bool:
    n = len(coloring)
    return all(coloring[i - 1] == coloring[n - i] for i in range(1, n + 1) if 3 * i != n + 1 and 3 * i != 2 * (n + 1))


def depth(coloring: list[int], k: int) -> int:
    first = {}
    for i, c in enumerate(coloring, start=1):
        first.setdefault(c, i)
    return max(first.values()) if len(first) == k else 0


def _solve_inner(case, queue=None):
    """Colouring of 1..n into k sum-free sets in the requested class, None if unsatisfiable.

    case keys: n, k, symmetric (bool), min_depth (int or 0): the last set has no element below min_depth,
    fixed (colours of 1..m), solver (pysat backend, default Cadical153), merge (symmetric only: partners share variables
    instead of equality clauses), size_bounds ((lo, hi): every set has between lo and hi elements; not with merge).
    """
    import pysat.solvers as ps

    n, k = case["n"], case["k"]
    sym = bool(case.get("symmetric"))
    merge = sym and bool(case.get("merge"))
    if merge and case.get("size_bounds"):
        raise ValueError("size_bounds counts elements, so it needs one variable per element: do not combine with merge")

    def rep(i: int) -> int:
        return min(i, n + 1 - i) if merge and 3 * i != n + 1 and 3 * i != 2 * (n + 1) else i

    var = lambda i, c: (rep(i) - 1) * k + c + 1  # noqa: E731
    s = getattr(ps, case.get("solver", "Cadical153"))()

    def add(lits) -> None:
        s.add_clause(list(dict.fromkeys(lits)))

    for i in range(1, n + 1):
        add([var(i, c) for c in range(k)])
        for c in range(k):
            for d in range(c + 1, k):
                add([-var(i, c), -var(i, d)])
    for c in range(k):
        for a in range(1, n + 1):
            for b in range(a, n + 1 - a):
                add([-var(a, c), -var(2 * a, c)] if a == b else [-var(a, c), -var(b, c), -var(a + b, c)])
    if sym and not merge:
        for i in range(1, n + 1):
            j = n + 1 - i
            if i < j and 3 * i != n + 1 and 3 * i != 2 * (n + 1):
                for c in range(k):
                    add([-var(i, c), var(j, c)])
                    add([var(i, c), -var(j, c)])
    for i in range(1, min(case.get("min_depth", 0), n + 1)):
        add([-var(i, k - 1)])
    if case.get("size_bounds"):
        from pysat.card import CardEnc

        lo, hi = case["size_bounds"]
        top = n * k
        for c in range(k):
            lits = [var(i, c) for i in range(1, n + 1)]
            for enc in (CardEnc.atleast(lits=lits, bound=lo, top_id=top), CardEnc.atmost(lits=lits, bound=hi, top_id=top)):
                top = max(top, enc.nv)
                for cl in enc.clauses:
                    add(cl)
    fixed = case.get("fixed")
    if fixed:  # colours of 1..len(fixed) are given (set names as in the list); implies the symmetric mirror when symmetric
        for i, c in enumerate(fixed, start=1):
            add([var(i, c)])
    else:
        add([var(1, 0)])  # name set 0 as the set holding 1 (sets are interchangeable)
    res = s.solve()
    if not res:
        s.delete()
        return None
    model = {x for x in s.get_model() if x > 0}
    s.delete()
    out = [next(c for c in range(k) if var(i, c) in model) for i in range(1, n + 1)]
    if queue is not None:
        queue.put(out)
    return out


def _child(case, queue) -> None:
    queue.put(_solve_inner(case))


def solve(case, timeout: float | None = None):
    """`_solve_inner` with a hard time limit: CaDiCaL (via pysat) ignores interrupt(), so a timeout kills a child process.

    Returns the colouring, None when unsatisfiable, or "timeout".
    """
    if not timeout:
        return _solve_inner(case)
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(target=_child, args=(case, queue))
    proc.start()
    try:
        return queue.get(timeout=timeout)
    except Exception:  # queue.Empty: the child is still searching
        return "timeout"
    finally:
        proc.terminate()
        proc.join()


def max_n(k: int, symmetric: bool, lo: int, hi: int, limit: float = 120) -> tuple[int, str]:
    """Largest n in [lo, hi] with a partition in the class (scan from hi down; reports the first SAT)."""
    for n in range(hi, lo - 1, -1):
        w = solve({"n": n, "k": k, "symmetric": symmetric}, timeout=limit)
        if w == "timeout":
            return n, "timeout"
        if w is not None:
            return n, "sat"
    return 0, "none"


# Known: S(1..5) = 1, 4, 13, 44, 160. (case, expected) with expected "witness" or "none".
CALIBRATION = [
    ({"n": 4, "k": 2}, "witness"), ({"n": 5, "k": 2}, "none"),
    ({"n": 13, "k": 3}, "witness"), ({"n": 14, "k": 3}, "none"),
    ({"n": 44, "k": 4}, "witness"), ({"n": 45, "k": 4}, "none"),
    ({"n": 44, "k": 4, "symmetric": True}, "witness"),  # 45 = 3*15: the paper's exception leaves 15 and 30 unpaired
    ({"n": 45, "k": 4, "symmetric": True}, "none"),
    ({"n": 160, "k": 5, "symmetric": True, "min_depth": 44, "timeout": 1400}, "witness"),  # S(5)=160 (Exoo; Heule 2017 proved it is exact)
]
OPEN: list[dict] = []


def run_cases(cases) -> list[dict]:
    rows = []
    for case, expected in cases:
        t = time.time()
        w = solve(case, timeout=case.get("timeout"))
        status = "timeout" if w == "timeout" else ("witness" if w is not None else "none")
        ok = status == "witness" and verify(w, case["k"]) and (not case.get("symmetric") or is_symmetric(w))
        rows.append({"case": case, "expected": expected, "status": status, "verified": ok,
                     "matches_expected": status == expected and (status != "witness" or ok), "seconds": round(time.time() - t, 2)})
        print(rows[-1], flush=True)
    return rows


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    rows = run_cases([*CALIBRATION, *[(c, "open") for c in OPEN]])
    (RESULTS / "schur-structured.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
