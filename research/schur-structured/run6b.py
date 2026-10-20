"""S(6) by extension: fix a 5-set partition of [1,160] (from the symmetric depth-44 class) and solve only the middle.

Usage: python run6b.py N [--parts M] [--multiply] [--budget CONFLICTS] [--prefix M] [--sample K] [--workers W] [--shard I/S] [--summary FILE]. Gate: N=536 (known to exist, Fredricksen and Sweet 2000).
A found witness is re-checked by verify() and is_symmetric() before it is written to results/s6b-N.json.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from pysat.solvers import Cadical153
from search import RESULTS, is_symmetric, verify


def five_set_parts(count: int) -> list[list[int]]:
    """Distinct symmetric depth>=44 5-set partitions of [1,160] (blocking each found one)."""
    n, k = 160, 5
    var = lambda i, c: (i - 1) * k + c + 1  # noqa: E731
    s = Cadical153()
    for i in range(1, n + 1):
        s.add_clause([var(i, c) for c in range(k)])
        for c in range(k):
            for d in range(c + 1, k):
                s.add_clause([-var(i, c), -var(i, d)])
    for c in range(k):
        for a in range(1, n + 1):
            for b in range(a, n + 1 - a):
                s.add_clause([-var(a, c), -var(2 * a, c)] if a == b else [-var(a, c), -var(b, c), -var(a + b, c)])
    for i in range(1, n + 1):
        j = n + 1 - i
        if i < j:
            for c in range(k):
                s.add_clause([-var(i, c), var(j, c)])
                s.add_clause([var(i, c), -var(j, c)])
    for i in range(1, 44):
        s.add_clause([-var(i, k - 1)])
    s.add_clause([var(1, 0)])
    for i in range(1, n + 1):  # name sets 0..3 by first appearance, so each partition is listed once (set 4 has the depth)
        for c in range(1, k - 1):
            s.add_clause([-var(i, c)] + [var(j, c - 1) for j in range(1, i)])
    out = []
    while len(out) < count and s.solve():
        model = {x for x in s.get_model() if x > 0}
        col = [next(c for c in range(k) if var(i, c) in model) for i in range(1, n + 1)]
        out.append(col)
        s.add_clause([-var(i, col[i - 1]) for i in range(1, n + 1)])  # block this exact colouring
    s.delete()
    return out


def multiples(parts: list[list[int]]) -> list[list[int]]:
    """All distinct images of the partitions under multiplication by a unit mod 161 (a and -a give the same image)."""
    from math import gcd
    seen, out = set(), []
    for col in parts:
        for m in range(1, 81):
            if gcd(m, 161) != 1:
                continue
            new = [0] * 160
            for i in range(1, 161):
                new[(m * i) % 161 - 1] = col[i - 1]
            order, canon = {}, []
            for c in new:  # the five base sets are interchangeable in the 6-set problem: name them by first appearance
                if c not in order:
                    order[c] = len(order)
                canon.append(order[c])
            key = tuple(canon)
            if key not in seen:
                seen.add(key)
                out.append(canon)
    return out


_SOLVER = None


def _init(n: int, budget: int) -> None:
    """Per-worker: build the symmetric 6-set instance for n once; each task only adds assumptions."""
    global _SOLVER
    from pysat.solvers import Glucose4
    k = 6
    var = lambda i, c: (i - 1) * k + c + 1  # noqa: E731
    s = Glucose4()
    for i in range(1, n + 1):
        s.add_clause([var(i, c) for c in range(k)])
        for c in range(k):
            for d in range(c + 1, k):
                s.add_clause([-var(i, c), -var(i, d)])
    for c in range(k):
        for a in range(1, n + 1):
            for b in range(a, n + 1 - a):
                s.add_clause([-var(a, c), -var(2 * a, c)] if a == b else [-var(a, c), -var(b, c), -var(a + b, c)])
    for i in range(1, n + 1):
        j = n + 1 - i
        if i < j and 3 * i != n + 1 and 3 * i != 2 * (n + 1):
            for c in range(k):
                s.add_clause([-var(i, c), var(j, c)])
                s.add_clause([var(i, c), -var(j, c)])
    _SOLVER = (s, n, budget)


def _task(col: list[int]):
    s, n, budget = _SOLVER
    k = 6
    s.conf_budget(budget)
    res = s.solve_limited(assumptions=[(i - 1) * k + col[i - 1] + 1 for i in range(1, len(col) + 1)])
    if res is None:
        return "timeout"
    if not res:
        return None
    model = {x for x in s.get_model() if x > 0}
    return [next(c for c in range(k) if (i - 1) * k + c + 1 in model) for i in range(1, n + 1)]


def main() -> None:
    from multiprocessing import Pool
    a = sys.argv[1:]
    n = int(next(x for x in a if x.isdigit()))
    get = lambda f, d: type(d)(a[a.index(f) + 1]) if f in a else d  # noqa: E731
    parts, budget, workers = get("--parts", 100), get("--budget", 200000), get("--workers", 24)
    t0 = time.time()
    P = five_set_parts(parts)
    if "--multiply" in a:
        P = multiples(P)
    if "--prefix" in a:  # keep only the first m entries (a 5-set partition of [1, m]); the rest is left to the solver
        m = get("--prefix", 160)
        P = [list(t) for t in dict.fromkeys(tuple(c[:m]) for c in P)]
    if "--sample" in a:
        import random
        random.Random(1).shuffle(P)
        P = P[: get("--sample", 1000)]
    if "--shard" in a:  # shard I of S by a hash of the colouring, so the split does not depend on generation order
        import hashlib
        idx, total = (int(x) for x in a[a.index("--shard") + 1].split("/"))
        P = [c for c in P if int(hashlib.md5(bytes(c)).hexdigest(), 16) % total == idx]
    print(f"{len(P)} partitions of 160 in {time.time() - t0:.1f}s", flush=True)
    stats = {"sat": 0, "unsat": 0, "timeout": 0}
    RESULTS.mkdir(exist_ok=True)
    with Pool(workers, initializer=_init, initargs=(n, budget)) as pool:
        for i, w in enumerate(pool.imap_unordered(_task, P, chunksize=8)):
            st = "timeout" if w == "timeout" else ("unsat" if w is None else "sat")
            stats[st] += 1
            if i % 2000 == 0 or st != "unsat":
                print(i, stats, round(time.time() - t0), "s", flush=True)
            if isinstance(w, list):
                ok = verify(w, 6) and is_symmetric(w)
                (RESULTS / f"s6b-{n}.json").write_text(json.dumps({"n": n, "verified": ok, "coloring": w}))
                print("WITNESS", n, "verified" if ok else "FAILED VERIFY", flush=True)
                if ok:
                    pool.terminate()
                    break
    print("done", stats, round(time.time() - t0), "s", flush=True)
    if "--summary" in a:
        Path(a[a.index("--summary") + 1]).write_text(json.dumps({"n": n, "starts": len(P), "stats": stats, "seconds": round(time.time() - t0)}))


if __name__ == "__main__":
    main()
