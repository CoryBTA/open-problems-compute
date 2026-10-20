#!/usr/bin/env python3
"""Resolution sweep of R = ||grad omega||^2 / ||grad|omega|||^2 with a stable integrator.

Differences from resolution_sweep.py (which used forward Euler, dt <= 0.005, and diverged at N=96):
  * integrating-factor RK4 (viscous term exact, nonlinear term 4th order in time)
  * time step from the CFL condition dt = CFL * dx / max|u|, recomputed every step
  * R is reported two ways: the original spectral |omega| gradient (R_spec) and a regularized one,
    grad sqrt(|omega|^2 + eps^2) with eps = 1e-2 * max|omega| (R_reg), because |omega| has cone
    singularities at zeros of omega that a spectral derivative resolves poorly
  * every run records whether it stayed finite, the peak max|omega|, steps and wall time
Config: Taylor-Green, amp=5, nu=0.01, T=2 (the 'TG5_Re500' configuration). CPU only (numpy).
Usage: ns_rk4_sweep.py [N ...]   -> writes ns_rk4_sweep.json next to this file.
"""
import importlib.util
import json
import sys
from pathlib import Path
from time import time as wall_time

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ns_ratio_R", HERE / "ns_ratio_R.py")
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

AMP, NU, T_END, CFL, EPS_REL = 5.0, 0.01, 2.0, 0.5, 1e-2


def nonlinear(s, w):
    """Dealiased nonlinear term of the vorticity equation, spectral space; w is a tuple of 3 spectral arrays."""
    wx_h, wy_h, wz_h = w
    ux_h, uy_h, uz_h = s.vel(wx_h, wy_h, wz_h)
    wx, wy, wz = s.p(wx_h), s.p(wy_h), s.p(wz_h)
    ux, uy, uz = s.p(ux_h), s.p(uy_h), s.p(uz_h)
    K = (s.KX, s.KY, s.KZ)
    dux = [s.g(ux_h, k) for k in K]
    duy = [s.g(uy_h, k) for k in K]
    duz = [s.g(uz_h, k) for k in K]
    dwx = [s.g(wx_h, k) for k in K]
    dwy = [s.g(wy_h, k) for k in K]
    dwz = [s.g(wz_h, k) for k in K]
    rx = wx * dux[0] + wy * dux[1] + wz * dux[2] - (ux * dwx[0] + uy * dwx[1] + uz * dwx[2])
    ry = wx * duy[0] + wy * duy[1] + wz * duy[2] - (ux * dwy[0] + uy * dwy[1] + uz * dwy[2])
    rz = wx * duz[0] + wy * duz[1] + wz * duz[2] - (ux * dwz[0] + uy * dwz[1] + uz * dwz[2])
    d = s.dealias
    return (np.fft.fftn(rx) * d, np.fft.fftn(ry) * d, np.fft.fftn(rz) * d), float(np.max(np.abs((ux, uy, uz))))


def ifrk4_step(s, w, dt):
    """One integrating-factor RK4 step; returns (new w, max|u| seen at the first stage)."""
    E = np.exp(-s.nu * s.K2 * dt / 2)
    E2 = E * E
    k1, umax = nonlinear(s, w)
    a = tuple(E * (wi + dt / 2 * ki) for wi, ki in zip(w, k1))
    k2, _ = nonlinear(s, a)
    b = tuple(E * wi + dt / 2 * ki for wi, ki in zip(w, k2))
    k3, _ = nonlinear(s, b)
    c = tuple(E2 * wi + dt * E * ki for wi, ki in zip(w, k3))
    k4, _ = nonlinear(s, c)
    new = tuple(E2 * wi + dt / 6 * (E2 * a1 + 2 * E * (a2 + a3) + a4)
                for wi, a1, a2, a3, a4 in zip(w, k1, k2, k3, k4))
    return new, umax


def measure(s, w):
    """R_spec (original estimator), R_reg (regularized |omega|), Omega, max|omega|."""
    wx, wy, wz = (s.p(c) for c in w)
    mag2 = wx ** 2 + wy ** 2 + wz ** 2
    mx = float(np.sqrt(np.max(mag2)))
    if not np.isfinite(mx) or mx > 1e6:
        return None
    dx3 = s.dx ** 3
    omega = float(np.sum(mag2) * dx3)
    grad_w = sum(float(np.sum(s.p(1j * K * c) ** 2) * dx3) for c in w for K in (s.KX, s.KY, s.KZ))

    def grad_mag(m):
        mh = np.fft.fftn(m)
        return sum(float(np.sum(s.p(1j * K * mh) ** 2) * dx3) for K in (s.KX, s.KY, s.KZ))

    g_spec = grad_mag(np.sqrt(mag2 + 1e-30))
    g_reg = grad_mag(np.sqrt(mag2 + (EPS_REL * mx) ** 2))
    return {"R_spec": grad_w / (g_spec + 1e-30), "R_reg": grad_w / (g_reg + 1e-30), "Omega": omega, "max_w": mx}


def run(N):
    s = base.NS3D(N, NU)
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    ux = AMP * np.sin(X) * np.cos(Y) * np.cos(Z)
    uy = -AMP * np.cos(X) * np.sin(Y) * np.cos(Z)
    w = tuple(np.asarray(c) for c in s.curl(np.fft.fftn(ux), np.fft.fftn(uy), np.fft.fftn(np.zeros_like(X))))
    t, steps, died = 0.0, 0, None
    rs, rr, peak = [], [], 0.0
    next_sample = 0.0
    u_prev = AMP
    t0 = wall_time()
    with np.errstate(all="ignore"):
        while t < T_END:
            if t >= next_sample:
                m = measure(s, w)
                if m is None:
                    died = round(t, 3)
                    break
                rs.append(m["R_spec"]); rr.append(m["R_reg"]); peak = max(peak, m["max_w"])
                next_sample += T_END / 100
            dt = min(CFL * s.dx / u_prev, T_END - t)
            w, u_prev = ifrk4_step(s, w, dt)
            u_prev = max(u_prev, 1e-6)
            t += dt
            steps += 1
    return {"N": N, "finite_to_T": died is None, "died_at_t": died, "steps": steps, "peak_max_omega": peak,
            "R_spec_max": max(rs) if rs else None, "R_reg_max": max(rr) if rr else None,
            "R_spec_mean": float(np.mean(rs)) if rs else None, "R_reg_mean": float(np.mean(rr)) if rr else None,
            "seconds": round(wall_time() - t0, 1)}


def main():
    sizes = [int(a) for a in sys.argv[1:]] or [24, 32, 48, 64, 96]
    out = []
    for N in sizes:
        r = run(N)
        out.append(r)
        print(r, flush=True)
        (HERE / "ns_rk4_sweep.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
