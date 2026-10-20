#!/usr/bin/env python3
"""Reproduce the finite-prime construction in Connes, "The Riemann Hypothesis: Past, Present and a Letter Through Time"
(arXiv 2602.04022, Feb 2026), the 'Letter' (pp. 24-25), from the description there:

  * phi lives on [1, X] (X = 13: only the primes and prime powers up to 13 enter); f(y) = phi(e^y) on [0, L], L = log X.
  * Q(phi) = Weil's explicit formula applied to psi(v) = int phi(u) phi(uv) du/u, i.e. to g(x) = int f(y) f(y+x) dy.
  * minimise Q subject to int phi(u)^2 du/u = int f^2 = 1; take the Mellin transform of the minimiser eta; its zeros
    (it is stated they lie on the critical line) approximate the first zeros of zeta.

Explicit formula used (Rudnick-Sarnak form, g even, h(r) = int g(x) e^{irx} dx, zeros rho = 1/2 + i gamma):
  sum_gamma h(gamma) = h(i/2) + h(-i/2) - g(0) log(pi) + (1/2pi) int h(r) Re psi(1/4 + ir/2) dr - 2 sum_n Lambda(n) n^(-1/2) g(log n)
with the archimedean term rewritten in x-space (digamma integral representation) so only g appears. The paper's own
normalisation was not available to me; the implementation is therefore checked against the true zeros (test 1) and the
construction against the paper's reported accuracy (test 2). Basis: orthonormal Legendre polynomials on [0, L]. Double
precision: expect far fewer digits than the paper's 55.
Usage: connes_weil_form.py [X] [n_basis]   -> prints a table, writes connes_weil_form.json
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.optimize import brentq

EULER_GAMMA = float(mp.euler)


def prime_powers(X):
    """(n, Lambda(n)) for prime powers 2 <= n <= X."""
    out = []
    for n in range(2, int(X) + 1):
        for p in range(2, n + 1):
            if n % p == 0:
                m = n
                while m % p == 0:
                    m //= p
                if m == 1:
                    out.append((n, float(np.log(p))))
                break
    return out


class Form:
    """Matrices of Q on the Legendre basis of [0, L]; A is the quadratic form, M the Gram matrix (identity)."""

    def __init__(self, X, n_basis, n_y=120, n_x=160):
        self.X, self.L, self.n = X, float(np.log(X)), n_basis
        L = self.L
        gx, gw = leg.leggauss(n_x)
        self.xq, self.wx = 0.5 * L * (gx + 1), 0.5 * L * gw
        self.yq_nodes, self.yq_w = leg.leggauss(n_y)
        self.pp = prime_powers(X)
        self.A = self._build()

    def basis(self, y):
        """Orthonormal Legendre polynomials on [0, L] evaluated at y (array); shape (n, len(y))."""
        t = 2 * y / self.L - 1
        V = leg.legvander(t, self.n - 1).T
        return V * np.sqrt((2 * np.arange(self.n) + 1) / self.L)[:, None]

    def corr(self, x):
        """G[k, l] = int_0^{L-x} b_k(y) b_l(y+x) dy for a single x in [0, L]."""
        span = self.L - x
        y = 0.5 * span * (self.yq_nodes + 1)
        w = 0.5 * span * self.yq_w
        return (self.basis(y) * w) @ self.basis(y + x).T

    def _build(self):
        n, L = self.n, self.L
        yq = 0.5 * L * (self.yq_nodes + 1)
        wq = 0.5 * L * self.yq_w
        B = self.basis(yq)
        u = B @ (wq * np.exp(yq / 2))     # int b_k e^{y/2}
        v = B @ (wq * np.exp(-yq / 2))    # int b_k e^{-y/2}
        g0 = np.eye(n)                    # g(0) matrix = Gram = identity
        T1 = np.outer(u, v) + np.outer(v, u)              # h(i/2) + h(-i/2) = 2 ab, symmetrised
        T2 = -np.log(np.pi) * g0
        integ = np.zeros((n, n))
        for x, w in zip(self.xq, self.wx, strict=True):
            G = self.corr(x)
            G = 0.5 * (G + G.T)
            integ += w * 2 * (g0 * np.exp(-2 * x) - np.exp(-x / 2) * G) / (1 - np.exp(-2 * x))
        T3 = -EULER_GAMMA * g0 + integ - g0 * np.log(1 - np.exp(-2 * L))
        T4 = np.zeros((n, n))
        for m, lam in self.pp:
            G = self.corr(float(np.log(m)))
            T4 += -2 * lam / np.sqrt(m) * 0.5 * (G + G.T)
        return T1 + T2 + T3 + T4

    def value(self, c):
        """Q(f) for f = sum c_k b_k."""
        return float(c @ self.A @ c)


def zeros_sum(form, c, n_zeros=60, n_nodes=3000):
    """2 * sum_{gamma>0} |fhat(gamma)|^2 over the first n_zeros true zeros (h(r) = |fhat(r)|^2).

    The Fourier integral is done with n_nodes Gauss-Legendre nodes (enough for gamma up to a few hundred); too few
    nodes alias at large gamma and make this side wrong (found the hard way).
    """
    gn, gw = leg.leggauss(n_nodes)
    yq, wq = 0.5 * form.L * (gn + 1), 0.5 * form.L * gw
    f = c @ form.basis(yq)
    total = 0.0
    for k in range(1, n_zeros + 1):
        gam = float(mp.zetazero(k).imag)
        total += abs(np.sum(wq * f * np.exp(1j * gam * yq))) ** 2
    return 2 * total


def mellin_zeros(form, c, rmax=130.0):
    """Real zeros of the Mellin transform of eta on the critical line (eta assumed even about L/2)."""
    yq = 0.5 * form.L * (form.yq_nodes + 1)
    wq = 0.5 * form.L * form.yq_w
    f = c @ form.basis(yq)

    def F(r):
        return float(np.sum(wq * f * np.cos(r * (yq - form.L / 2))))

    grid = np.linspace(0.5, rmax, 6000)
    vals = np.array([F(r) for r in grid])
    return [brentq(F, grid[i], grid[i + 1], xtol=1e-14) for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0]


def main():
    X = float(sys.argv[1]) if len(sys.argv) > 1 else 13.0
    n_basis = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    form = Form(X, n_basis)
    # test 1: the explicit-formula side equals the zero side for a test function
    rng = np.random.default_rng(0)
    c = rng.standard_normal(n_basis) / (1 + np.arange(n_basis)) ** 2
    c = c / np.linalg.norm(c)
    lhs, rhs = form.value(c), zeros_sum(form, c)
    print(f"test 1 (explicit formula vs 60 true zeros, random smooth f): Q = {lhs:.10f}, zero side = {rhs:.10f}, diff = {abs(lhs - rhs):.2e}")
    # minimise Q on the unit sphere
    w, V = np.linalg.eigh(0.5 * (form.A + form.A.T))
    eta = V[:, 0]
    print(f"lowest eigenvalue {w[0]:.8f}, next {w[1]:.8f}, {w[2]:.8f}")
    zs = np.array(mellin_zeros(form, eta))
    true = [float(mp.zetazero(k).imag) for k in range(1, 51)]
    rows = []
    print(f"{len(zs)} real zeros of the Mellin transform below r = 130; each true zero is matched to the nearest one")
    print(f"{'k':>3} {'nearest Mellin zero':>20} {'true zero':>16} {'|diff|':>10}")
    for k, t in enumerate(true, 1):
        z = float(zs[np.argmin(abs(zs - t))])
        if k <= 10 or k % 10 == 0:
            print(f"{k:>3} {z:20.12f} {t:16.10f} {abs(z - t):10.2e}")
        rows.append({"k": k, "mellin": z, "true": t, "diff": abs(z - t)})
    extra = [float(z) for z in zs if min(abs(np.array(true) - z)) > 0.5]
    print(f"Mellin zeros not within 0.5 of any of the first 50 true zeros: {len(extra)}, e.g. {[round(z, 3) for z in extra[:6]]}")
    out = {"X": X, "n_basis": n_basis, "test1_Q": lhs, "test1_zero_side": rhs, "eigs": w[:4].tolist(), "zeros": rows,
           "n_found": int(len(zs)), "extra_mellin_zeros": extra}
    Path(__file__).with_name("connes_weil_form.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
