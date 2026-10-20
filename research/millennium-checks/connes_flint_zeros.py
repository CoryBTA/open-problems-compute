#!/usr/bin/env python3
"""Ball-arithmetic (python-flint / Arb) check of the Mellin zeros in connes_weil_form.Form.

What is certified: for f = sum c_k b_k with c_k the EXACT binary values of the double-precision lowest eigenvector of
Form.A (n_basis Legendre functions on [0, log X]), F(r) = int_0^L f(y) cos(r (y - L/2)) dy is evaluated with
acb.integral (rigorous), and a sign change of F across [z - d, z + d] proves a zero of F inside it (F continuous).
Not certified: that this f is the minimiser of Q (A and the eigenvector are double precision), nor that the zeros of F
approximate zeta zeros; the true zeros are Arb balls (acb.zeta_zero), compared with the double-precision Mellin zero.
Usage: connes_flint_zeros.py [n_basis] -> connes_flint_zeros.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx

sys.path.insert(0, str(Path(__file__).parent))
from connes_weil_form import Form, mellin_zeros  # noqa: E402

ctx.prec = 200


def F_ball(c_exact, L, r):
    """Rigorous enclosure of F(r) via acb.integral; c_exact = list of arb (exact doubles)."""
    n = len(c_exact)
    Lb = arb(L)  # L = log X evaluated as a double; the same L is used by Form, so the function class is the same

    def integrand(y, analytic):
        t = 2 * y / Lb - 1
        p0, p1 = acb(1), t
        s = c_exact[0] * p0 * acb(1 / L).sqrt() + (c_exact[1] * p1 * acb(3 / L).sqrt() if n > 1 else 0)
        for k in range(1, n - 1):
            p0, p1 = p1, ((2 * k + 1) * t * p1 - k * p0) / (k + 1)
            s += c_exact[k + 1] * p1 * acb((2 * k + 3) / L).sqrt()
        return s * (r * (y - Lb / 2)).cos()

    return acb.integral(integrand, 0, Lb).real


def run(Xs, nb, zeros_k=5, delta=1e-7):
    true_balls = [acb.zeta_zero(k).imag for k in range(1, zeros_k + 1)]
    out = []
    for X in Xs:
        t0 = time.time()
        form = Form(float(X), nb)
        w, V = np.linalg.eigh(0.5 * (form.A + form.A.T))
        eta = V[:, 0]
        zs = np.array(mellin_zeros(form, eta))
        c_exact = [arb(float(v)) for v in eta]
        rows = []
        for k, tb in enumerate(true_balls, 1):
            tmid = float(tb.mid())
            z = float(zs[np.argmin(abs(zs - tmid))])
            lo, hi = F_ball(c_exact, form.L, arb(z - delta)), F_ball(c_exact, form.L, arb(z + delta))
            bracket = (lo.upper() < 0 and hi.lower() > 0) or (lo.lower() > 0 and hi.upper() < 0)
            rows.append({"k": k, "double_mellin_zero": z, "bracket_halfwidth": delta, "F_lo": lo.str(5), "F_hi": hi.str(5),
                         "zero_certified_in_bracket": bool(bracket), "true_zero_ball_radius": float(tb.rad()),
                         "true_zero": tb.str(20, radius=False), "abs_diff_double": abs(z - tmid)})
            print(X, nb, k, f"diff={abs(z - tmid):.2e}", "bracket", bracket, flush=True)
        out.append({"X": X, "n_basis": nb, "lowest_eig": float(w[0]), "seconds": round(time.time() - t0, 1), "zeros": rows})
    return out


if __name__ == "__main__":
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 72
    res = {"n_basis": nb, "prec_bits": ctx.prec, "runs": run((13, 17, 23, 29, 37, 47), nb)}
    Path(__file__).with_name("connes_flint_zeros.json").write_text(json.dumps(res, indent=1))
