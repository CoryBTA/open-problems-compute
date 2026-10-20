#!/usr/bin/env python3
"""Sanity test for BSD numerics: recover the known analytic rank of elliptic curves of prime conductor.

Curves (Cremona labels, minimal models y^2 + y = x^3 + a2 x^2 + a4 x + a6):
  11a rank 0, 37a rank 1, 389a rank 2, 5077a rank 3.
Method: coefficients a_n from point counts (a_p = p + 1 - #E(F_p), singular curve counted for p = N);
completed L-function Lambda(s) by the approximate functional equation with incomplete gamma functions;
order of vanishing at s = 1 read from the growth of |Lambda(1 + h)| ~ h^r.
Also checks the Hasse-Weil bound a_p^2 <= 4p (the Riemann Hypothesis for these curves over F_p)
and, for 11a, the BSD formula L(E,1)/Omega = Sha * c_p / |tors|^2 = 1/5.
CPU only, runs in seconds.
"""
import json
import math
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 30

CURVES = {
    "11a": ((0, -1, 1, -10, -20), 11, 0),
    "37a": ((0, 0, 1, -1, 0), 37, 1),
    "389a": ((0, 1, 1, -2, 0), 389, 2),
    "5077a": ((0, 0, 1, -7, 6), 5077, 3),
}


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i, v in enumerate(sieve) if v]


def ap(model, p):
    """a_p = p + 1 - #E(F_p) for y^2 + a3 y = x^3 + a2 x^2 + a4 x + a6 (a1 = 0)."""
    _a1, a2, a3, a4, a6 = model
    ycount = {}
    for y in range(p):
        v = (y * y + a3 * y) % p
        ycount[v] = ycount.get(v, 0) + 1
    count = 1 + sum(ycount.get((x**3 + a2 * x * x + a4 * x + a6) % p, 0) for x in range(p))
    return p + 1 - count


def coeffs(model, N, nmax):
    ps = primes_upto(nmax)
    a_p = {p: ap(model, p) for p in ps}
    a = [0] * (nmax + 1)
    a[1] = 1
    for n in range(2, nmax + 1):
        m, p = n, next(q for q in ps if n % q == 0)
        k = 0
        while m % p == 0:
            m //= p
            k += 1
        if m > 1:
            a[n] = a[p**k] * a[m]
            continue
        if k == 1:
            a[n] = a_p[p]
        elif p == N:
            a[n] = a_p[p] ** k
        elif k == 2:
            a[n] = a_p[p] * a_p[p] - p
        else:
            a[n] = a_p[p] * a[p ** (k - 1)] - p * a[p ** (k - 2)]
    return a, a_p


def lam(s, a, N, w):
    r = mp.sqrt(N) / (2 * mp.pi)
    total = mp.mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0:
            continue
        x = 2 * mp.pi * n / mp.sqrt(N)
        total += a[n] * (
            r**s * mp.gammainc(s, x) / mp.mpf(n) ** s
            + w * r ** (2 - s) * mp.gammainc(2 - s, x) / mp.mpf(n) ** (2 - s)
        )
    return total


def main():
    out = {}
    for label, (model, N, known_rank) in CURVES.items():
        nmax = int(12 * math.sqrt(N)) + 20
        a, a_p = coeffs(model, N, nmax)
        w = ap(model, N)  # root number of a prime-conductor curve equals a_N (+1 split, -1 nonsplit)
        hasse_ok = all(v * v <= 4 * p for p, v in a_p.items() if p != N)
        hs = [mp.mpf("0.02"), mp.mpf("0.04"), mp.mpf("0.08")]
        vals = [abs(lam(1 + h, a, N, w)) for h in hs]
        # |Lambda(1+h)| ~ c h^r  =>  r = log(v2/v1) / log(h2/h1)
        r_est = [float(mp.log(vals[i + 1] / vals[i]) / mp.log(hs[i + 1] / hs[i])) for i in range(2)]
        out[label] = {
            "N": N, "known_rank": known_rank, "root_number": w, "terms": nmax,
            "hasse_weil_ok_for_all_good_p<=nmax": hasse_ok,
            "order_estimates": r_est,
            "estimated_rank": round(sum(r_est) / 2),
            "parity_matches_root_number": (known_rank % 2 == 0) == (w == 1),
        }
        print(label, out[label])
    # BSD leading term for 11a: L(E,1) = 0.2538418608559107, Omega = 1.2692093042795534, ratio 1/5
    a11, _ = coeffs(CURVES["11a"][0], 11, 120)
    L1 = 2 * sum(a11[n] / n * mp.e ** (-2 * mp.pi * n / mp.sqrt(11)) for n in range(1, len(a11)))
    out["11a_L1"] = float(L1)
    out["11a_L1_over_Omega"] = float(L1 / mp.mpf("1.2692093042795534"))
    print("11a L(E,1) =", float(L1), " L/Omega =", out["11a_L1_over_Omega"], "(BSD predicts 0.2)")
    Path(__file__).with_name("bsd_rank_check.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    sys.exit(main())
