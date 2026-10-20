#!/usr/bin/env python3
"""Real period of the BSD-4 curves in Arb ball arithmetic (python-flint), compared with the mpmath quadrature in
bsd_leading_term.real_period. Same quantity: Omega = integral of |dx/Y| over both branches of the real locus of
Y^2 = 4x^3 + b2 x^2 + 2 b4 x + b6, via AGM:  int_{e3}^inf dx/sqrt(g) = pi / (2 agm(sqrt(e3-e1), sqrt(e3-e2))),
int_{e1}^{e2} dx/sqrt(g) = the same value (K(k) form, k^2 = (e2-e1)/(e3-e1)) for three real roots e1<e2<e3 of g; with one real
root e and a conjugate pair e', e'': int_e^inf = pi / (2 agm(sqrt(e-e'), sqrt(e-e''))).
Writes bsd_flint_period.json: per curve the ball (midpoint, radius) of Omega, the stored mpmath Omega, relative difference,
and L1/Omega from the stored L1 with the Arb Omega.
"""
import json
import sys
from pathlib import Path

from flint import acb, acb_poly, arb, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bsd_leading_term as bl  # noqa: E402

ctx.prec = 200


def omega_ball(a2, a4, a6):
    b2, b4, b6, *_ = bl.invariants(a2, a4, a6)
    roots = acb_poly([b6, 2 * b4, b2, 4]).roots(tol=arb(10) ** -60, maxprec=4000)
    pi = arb.pi()
    real = [r for r in roots if r.imag.contains(0)]
    # exact discriminant of g (integers b2, b4, b6): > 0 iff three real roots
    disc = -4 * b2 ** 3 * b6 + b2 ** 2 * (2 * b4) ** 2 - 4 * 4 * (2 * b4) ** 3 - 27 * 16 * b6 ** 2 + 18 * 4 * b2 * (2 * b4) * b6
    if (disc > 0) != (len(real) == 3):
        raise RuntimeError(f'root classification disagrees with discriminant {disc} for {(a2, a4, a6)}')
    if len(real) == 3:
        e = sorted((r.real for r in real), key=lambda x: float(x.mid()))
        e1, e2, e3 = e
        a = (e3 - e1).sqrt()
        return 2 * pi / a.agm((e3 - e2).sqrt())  # both components give pi / (2 agm(sqrt(e3-e1), sqrt(e3-e2)))
    e = real[0].real
    others = [r for r in roots if not r.imag.contains(0)]
    s1, s2 = (e - others[0]).sqrt(), (e - others[1]).sqrt()
    return (2 * pi / (2 * acb(s1).agm(s2))).real


def main():
    src = json.loads((HERE / "bsd_sha_search.json").read_text())
    rows = []
    for c in src["curves"]:
        _, a2, _, a4, a6 = c["model"]
        om = omega_ball(a2, a4, a6)
        rel = abs(float(om.mid()) - c["Omega"]) / c["Omega"]
        rows.append({"model": c["model"], "N": c["N"], "Omega_flint": om.str(25, radius=False), "Omega_radius": float(om.rad()),
                     "Omega_mpmath": c["Omega"], "rel_diff": rel, "L1": c["L1"], "L_over_Omega_flint": c["L1"] / float(om.mid())})
        print(c["N"], c["model"], f"rel_diff={rel:.1e} rad={float(om.rad()):.1e}", flush=True)
    out = {"prec_bits": ctx.prec, "n_curves": len(rows), "max_rel_diff": max(r["rel_diff"] for r in rows), "curves": rows}
    (HERE / "bsd_flint_period.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
