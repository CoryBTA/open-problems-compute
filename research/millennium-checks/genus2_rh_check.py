#!/usr/bin/env python3
"""RH for genus-2 curves over finite fields: the check named in the RH finite-field template.

For y^2 = f(x) over F_p with f monic, squarefree, degree 5 (one point at infinity over every F_q):
  S_n = N_n - 1 - p^n,  P(T) = 1 + c1 T + c2 T^2 + p c1 T^3 + p^2 T^4,  c1 = S_1,  c2 = (S_1^2 + S_2) / 2
RH for the curve (Weil) says every root alpha of x^4 + c1 x^3 + c2 x^2 + p c1 x + p^2 has |alpha| = sqrt(p).
This reproduces a theorem; its purpose is a validated harness that a proposed method can be compared against.
N_1 counts points over F_p, N_2 over F_{p^2} (represented as a + b t with t^2 = r, r a non-residue mod p).
CPU only; seconds. Writes genus2_rh_check.json.
"""
import json
import random
from pathlib import Path

import numpy as np


def nonresidue(p):
    return next(r for r in range(2, p) if pow(r, (p - 1) // 2, p) == p - 1)


def poly_gcd_is_one(f, p):
    """True when f and f' are coprime over F_p (f squarefree). f: coefficient list, low degree first."""
    def trim(a):
        while a and a[-1] % p == 0:
            a = a[:-1]
        return a

    def mod(a, b):
        a = a[:]
        inv = pow(b[-1], -1, p)
        while len(a) >= len(b):
            c = a[-1] * inv % p
            for i in range(len(b)):
                a[len(a) - len(b) + i] = (a[len(a) - len(b) + i] - c * b[i]) % p
            a = trim(a)
            if not a:
                break
        return a

    a = trim([c % p for c in f])
    b = trim([(i * f[i]) % p for i in range(1, len(f))])
    while b:
        a, b = b, mod(a, b)
    return len(a) == 1


def count_over_Fp(f, p):
    sq = {}
    for y in range(p):
        sq[y * y % p] = sq.get(y * y % p, 0) + 1
    return 1 + sum(sq.get(sum(c * pow(x, i, p) for i, c in enumerate(f)) % p, 0) for x in range(p))


def count_over_Fp2(f, p):
    r = nonresidue(p)

    def mul(u, v):
        return ((u[0] * v[0] + r * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p)

    def power(u, e):
        out, base = (1, 0), u
        while e:
            if e & 1:
                out = mul(out, base)
            base = mul(base, base)
            e >>= 1
        return out

    count = 1  # point at infinity
    for a in range(p):
        for b in range(p):
            x = (a, b)
            val, xp = (0, 0), (1, 0)
            for c in f:
                val = ((val[0] + c * xp[0]) % p, (val[1] + c * xp[1]) % p)
                xp = mul(xp, x)
            if val == (0, 0):
                count += 1
            elif power(val, (p * p - 1) // 2) == (1, 0):
                count += 2
    return count


def check(f, p):
    n1, n2 = count_over_Fp(f, p), count_over_Fp2(f, p)
    s1, s2 = n1 - 1 - p, n2 - 1 - p * p
    assert (s1 * s1 + s2) % 2 == 0, "c2 must be an integer"
    c1, c2 = s1, (s1 * s1 + s2) // 2
    roots = np.roots([1, c1, c2, p * c1, p * p])
    return {"p": p, "f": f, "N1": n1, "N2": n2, "c1": c1, "c2": c2,
            "root_moduli_over_sqrt_p": sorted(round(float(abs(z)) / p ** 0.5, 9) for z in roots)}


def main():
    rng = random.Random(0)
    out = []
    for p in (5, 7, 11, 13):
        made = 0
        while made < 6:
            f = [rng.randrange(p) for _ in range(5)] + [1]
            if poly_gcd_is_one(f, p):
                out.append(check(f, p))
                made += 1
    worst = max(abs(m - 1) for r in out for m in r["root_moduli_over_sqrt_p"])
    print(len(out), "curves, worst | |alpha|/sqrt(p) - 1 | =", worst)
    Path(__file__).with_name("genus2_rh_check.json").write_text(json.dumps({"curves": out, "worst_deviation": worst}, indent=1))
    return 0 if worst < 1e-6 else 1


if __name__ == "__main__":
    raise SystemExit(main())
