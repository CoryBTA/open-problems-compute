#!/usr/bin/env python3
"""BSD leading-term check on rank-0 elliptic curves of prime conductor, with every ingredient computed here.

For E: y^2 + y = x^3 + a2 x^2 + a4 x + a6 (a1 = 0, a3 = 1) found by a coefficient search with |Delta| = p^k, p not dividing c4
(then E is minimal and semistable at p with conductor N = p). For rank 0 BSD predicts
    L(E,1) / Omega = |Sha| * c_p / |E(Q)_tors|^2          (c_infinity is inside Omega = int over E(R) of |dx / (2y+1)|)
so |Sha|_pred = (L(E,1)/Omega) * |tors|^2 / c_p must be a perfect square integer. Ingredients:
  * L(E,1) = 2 sum a_n / n exp(-2 pi n / sqrt N) (root number +1), a_n from point counts (bsd_rank_check.py)
  * Omega by numerical integration over the real locus of Y^2 = 4x^3 + b2 x^2 + 2 b4 x + b6, Y = 2y + 1
  * c_p = k if split (a_p = +1); nonsplit (a_p = -1): 1 if k odd, 2 if k even
  * |tors| from the gcd over good primes l >= 3 of #E(F_l) (an upper bound that is usually exact); reported with the squares test
CPU only, a few seconds. Writes bsd_leading_term.json.
"""
import importlib.util
import json
import math
import sys
from math import gcd
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("bsd_rank_check", HERE / "bsd_rank_check.py")
bsd = importlib.util.module_from_spec(spec)
sys.modules["bsd_rank_check"] = bsd
spec.loader.exec_module(bsd)
mp.mp.dps = 25


def invariants(a2, a4, a6):
    """(b2, b4, b6, Delta, c4) for a1 = 0, a3 = 1."""
    b2, b4, b6 = 4 * a2, 2 * a4, 1 + 4 * a6
    b8 = 4 * a2 * a6 + a2 - a4 * a4
    delta = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    c4 = b2 * b2 - 24 * b4
    return b2, b4, b6, delta, c4


def prime_power(n):
    """(p, k) when n = p^k for a prime p, else None."""
    for p in range(2, int(math.isqrt(n)) + 1):
        if n % p == 0:
            k = 0
            while n % p == 0:
                n //= p
                k += 1
            return (p, k) if n == 1 else None
    return (n, 1) if n > 1 else None


def real_period(b2, b4, b6):
    """Integral of |dx / Y| over the real locus of Y^2 = 4x^3 + b2 x^2 + 2 b4 x + b6 (both branches of Y)."""
    def g(x):
        return 4 * x ** 3 + b2 * x ** 2 + 2 * b4 * x + b6

    roots = sorted(r.real for r in mp.polyroots([4, b2, 2 * b4, b6], maxsteps=200, extraprec=50) if abs(r.imag) < mp.mpf(10) ** -15)
    if len(roots) == 3:
        e1, e2, e3 = roots
        parts = [mp.quad(lambda x: 1 / mp.sqrt(abs(g(x))), [e1, (e1 + e2) / 2, e2]), mp.quad(lambda x: 1 / mp.sqrt(abs(g(x))), [e3, e3 + 1, mp.inf])]
    else:
        # one real root: a complex pair close to the real axis gives a narrow peak of 1/sqrt(g) at its real part, so
        # put quadrature breakpoints around it at multiples of its imaginary part
        allr = mp.polyroots([4, b2, 2 * b4, b6], maxsteps=500, extraprec=100)
        cplx = [r for r in allr if abs(mp.im(r)) >= mp.mpf(10) ** -15]
        pts = [roots[0] + 1]
        if cplx:
            cre, w = mp.re(cplx[0]), max(abs(mp.im(cplx[0])), mp.mpf(10) ** -8)
            pts += [cre + sgn * k * w for k in (1, 3, 10, 30, 100) for sgn in (-1, 1)]
        pts = sorted(set(p for p in pts if p > roots[0]))
        parts = [mp.quad(lambda x: 1 / mp.sqrt(abs(g(x))), [roots[0]] + pts + [mp.inf])]
    return mp.re(2 * sum(parts))


def analyse(a2, a4, a6, p, k, max_n=400):
    model = (0, a2, 1, a4, a6)
    N = p
    nmax = int(12 * math.sqrt(N)) + 20
    a, a_p = bsd.coeffs(model, N, nmax)
    w = bsd.ap(model, N)
    if w != 1:
        return None
    L1 = 2 * sum(a[n] / n * mp.e ** (-2 * mp.pi * n / mp.sqrt(N)) for n in range(1, nmax + 1))
    if abs(L1) < 1e-8:
        return None
    b2, b4, b6, delta, c4 = invariants(a2, a4, a6)
    omega = real_period(b2, b4, b6)
    c_p = k if w == 1 else (1 if k % 2 else 2)
    tors = 0
    for l in bsd.primes_upto(80):
        if l >= 3 and N % l:
            tors = gcd(tors, l + 1 - bsd.ap(model, l))
    ratio = L1 / omega
    sha = ratio * tors * tors / c_p
    return {"model": [0, a2, 1, a4, a6], "N": N, "delta": delta, "k": k, "split": w == 1, "c_p": c_p, "torsion_bound": tors,
            "L1": float(L1), "Omega": float(omega), "L_over_Omega": float(ratio), "sha_pred": float(sha),
            "sha_is_square": abs(sha - round(sha)) < 1e-6 and round(math.isqrt(max(round(sha), 0)) ** 2) == round(sha)}


def main():
    found, seen = [], set()
    for a2 in (-1, 0, 1):
        for a4 in range(-40, 41):
            for a6 in range(-80, 81):
                _, _, _, delta, c4 = invariants(a2, a4, a6)
                if delta == 0:
                    continue
                pk = prime_power(abs(delta))
                if not pk or c4 % pk[0] == 0 or pk[0] > 400 or pk[0] < 11:
                    continue
                p, k = pk
                key = (p, tuple(bsd.ap((0, a2, 1, a4, a6), q) for q in (2, 3, 5, 7, 13)))
                if key in seen:
                    continue
                seen.add(key)
                r = analyse(a2, a4, a6, p, k)
                if r:
                    found.append(r)
    found.sort(key=lambda r: (r["N"], r["model"]))
    print(f"{len(found)} rank-0 prime-conductor curves")
    print(f"{'N':>4} {'k':>2} {'c_p':>3} {'tors<=':>6} {'Omega':>10} {'L(E,1)':>10} {'L/Omega':>9} {'Sha_pred':>9} square")
    for r in found:
        print(f"{r['N']:>4} {r['k']:>2} {r['c_p']:>3} {r['torsion_bound']:>6} {r['Omega']:>10.6f} {r['L1']:>10.6f} {r['L_over_Omega']:>9.5f} {r['sha_pred']:>9.4f} {r['sha_is_square']}")
    Path(__file__).with_name("bsd_leading_term.json").write_text(json.dumps({"curves": found}, indent=1))


if __name__ == "__main__":
    main()
