"""Benchmark SAT backends on Schur instances with known answers. Usage: python bench.py [timeout_s]"""
from __future__ import annotations

import multiprocessing as mp
import sys
import time

CASES = [
    ("S(4)=44 sat", {"n": 44, "k": 4, "symmetric": False, "min_depth": 0}, True),
    ("S(4)=45 unsat", {"n": 45, "k": 4, "symmetric": False, "min_depth": 0}, False),
    ("S(5)=160 sym depth44", {"n": 160, "k": 5, "symmetric": True, "min_depth": 44}, True),
    ("S(5)=160 sym", {"n": 160, "k": 5, "symmetric": True, "min_depth": 0}, True),
    ("S(5)=161 sym depth44", {"n": 161, "k": 5, "symmetric": True, "min_depth": 44}, False),
]
SOLVERS = ["Cadical153", "Cadical195", "Kissat404", "Glucose4", "MapleChrono"]


def build(solver_name, case):
    import pysat.solvers as ps
    n, k = case["n"], case["k"]
    var = lambda i, c: (i - 1) * k + c + 1  # noqa: E731
    s = getattr(ps, solver_name)()
    for i in range(1, n + 1):
        s.add_clause([var(i, c) for c in range(k)])
        for c in range(k):
            for d in range(c + 1, k):
                s.add_clause([-var(i, c), -var(i, d)])
    for c in range(k):
        for a in range(1, n + 1):
            for b in range(a, n + 1 - a):
                s.add_clause([-var(a, c), -var(2 * a, c)] if a == b else [-var(a, c), -var(b, c), -var(a + b, c)])
    if case["symmetric"]:
        for i in range(1, n + 1):
            j = n + 1 - i
            if i < j and 3 * i != n + 1 and 3 * i != 2 * (n + 1):
                for c in range(k):
                    s.add_clause([-var(i, c), var(j, c)])
                    s.add_clause([var(i, c), -var(j, c)])
    for i in range(1, case["min_depth"]):
        s.add_clause([-var(i, k - 1)])
    s.add_clause([var(1, 0)])
    return s


def child(solver_name, case, q):
    t = time.time()
    s = build(solver_name, case)
    q.put((bool(s.solve()), round(time.time() - t, 2)))


if __name__ == "__main__":
    limit = float(sys.argv[1]) if len(sys.argv) > 1 else 120
    ctx = mp.get_context("spawn")
    print(f"{'case':26s}" + "".join(f"{s:>13s}" for s in SOLVERS), flush=True)
    for label, case, expect in CASES:
        row = f"{label:26s}"
        for sv in SOLVERS:
            q = ctx.Queue()
            p = ctx.Process(target=child, args=(sv, case, q))
            p.start()
            try:
                res, secs = q.get(timeout=limit)
                row += f"{('OK ' if res == expect else 'WRONG ') + str(secs) + 's':>13s}"
            except Exception:
                row += f"{'>' + str(int(limit)) + 's':>13s}"
            p.terminate()
            p.join()
        print(row, flush=True)
