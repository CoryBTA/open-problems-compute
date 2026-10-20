#!/usr/bin/env python3
"""Does R = ||grad w||^2 / ||grad|w|||^2 stay bounded for strongly evolved flows, and is Step 6's constant uniform?

Uses the stable IFRK4 + CFL solver of ns_rk4_sweep.py on several initial conditions (not only Taylor-Green) at two or
more grid sizes, and records for each run:
  * R(t) (spectral |w| gradient), F(t) = 4 ||grad|w|||^2 / Omega, Omega(t), max|w|(t)
  * how much the flow evolved: max Omega / Omega(0), max|w| / max|w|(0)
  * Step 6 constant: C(t) = (dF/dt) / (F log F) from finite differences, reported overall and after t = 0.1
Usage: ns_ic_study.py N [N ...]   -> writes ns_ic_study.json. CPU only (numpy).
"""
import importlib.util
import json
import sys
from pathlib import Path
from time import time as wall_time

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


base = load("ns_ratio_R", "ns_ratio_R.py")
rk = load("ns_rk4_sweep", "ns_rk4_sweep.py")
NU, CFL = 0.01, 0.5
AMP = 1.0  # multiplies the initial velocity (amplitude scan); 1.0 reproduces all earlier runs


def velocity_ic(name, X, Y, Z, N):
    """Velocity field components for the named initial condition."""
    if name == "tg5":
        return 5 * np.sin(X) * np.cos(Y) * np.cos(Z), -5 * np.cos(X) * np.sin(Y) * np.cos(Z), np.zeros_like(X)
    if name == "roberts":  # u = (sin z, cos z, 0): Beltrami (curl u = u) with |u| = |w| = 1 constant
        return np.sin(Z), np.cos(Z), np.zeros_like(X)
    if name == "abc":
        A = B = C = 2.0
        return A * np.sin(Z) + C * np.cos(Y), B * np.sin(X) + A * np.cos(Z), C * np.sin(Y) + B * np.cos(X)
    if name == "kida":
        return 3 * (np.sin(X) * np.cos(Y) * np.cos(Z) + np.cos(X) * np.sin(Z)), 3 * (-np.cos(X) * np.sin(Y) * np.cos(Z) + np.sin(Y) * np.cos(Z)), 3 * np.cos(Z) * np.sin(X)
    if name.startswith("robp"):  # "robp_<k>_<eps>": Roberts at wavenumber k plus eps * Taylor-Green-like field (step7_family_scan.py)
        _, k, eps = name.split("_")
        k, eps = int(k), float(eps)
        return np.sin(k * Z) + eps * np.sin(X) * np.cos(Y) * np.cos(Z), np.cos(k * Z) - eps * np.cos(X) * np.sin(Y) * np.cos(Z), np.zeros_like(X)
    if name.startswith("rand"):
        seed = int(name[4:])
        rng = np.random.default_rng(seed)
        k1 = np.fft.fftfreq(N, d=1.0 / N)
        KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
        K = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
        env = np.where((K > 0) & (K <= 4), K ** -2.0, 0.0)  # only modes |k| <= 4, so every grid resolves the same field
        hats = [env * (rng.standard_normal((N, N, N)) + 1j * rng.standard_normal((N, N, N))) for _ in range(3)]
        kk = np.where(K == 0, 1, K ** 2)
        div = (KX * hats[0] + KY * hats[1] + KZ * hats[2]) / kk  # project to divergence free
        hats = [hats[0] - KX * div, hats[1] - KY * div, hats[2] - KZ * div]
        u = [np.real(np.fft.ifftn(h)) for h in hats]
        scale = 4.0 / max(float(np.max(np.abs(u))), 1e-12)
        return tuple(scale * c for c in u)
    raise ValueError(name)


def measure(s, w):
    """R, F, Omega, max|w| at one time; None when the field is not finite."""
    wx, wy, wz = (s.p(c) for c in w)
    mag2 = wx ** 2 + wy ** 2 + wz ** 2
    mx = float(np.sqrt(np.max(mag2)))
    if not np.isfinite(mx) or mx > 1e6:
        return None
    dx3 = s.dx ** 3
    omega = float(np.sum(mag2) * dx3)
    grad_w = sum(float(np.sum(s.p(1j * K * c) ** 2) * dx3) for c in w for K in (s.KX, s.KY, s.KZ))
    mh = np.fft.fftn(np.sqrt(mag2 + 1e-30))
    g = sum(float(np.sum(s.p(1j * K * mh) ** 2) * dx3) for K in (s.KX, s.KY, s.KZ))
    return {"R": grad_w / (g + 1e-30), "F": 4 * g / omega, "Omega": omega, "max_w": mx}


def _save(path, w, t, steps, nxt, u_prev, ts, ms):
    """Write a resumable snapshot of a run (atomic replace)."""
    tmp = str(path) + ".tmp.npz"
    np.savez(tmp, w0=w[0], w1=w[1], w2=w[2], scal=np.array([t, steps, nxt, u_prev]), ts=np.array(ts), ms=np.array([json.dumps(m) for m in ms]))
    Path(tmp).replace(path)


def _load(path):
    """Read a snapshot written by _save; returns the same values as the run loop state."""
    z = np.load(path)
    t, steps, nxt, u_prev = (float(v) for v in z["scal"])
    return (z["w0"], z["w1"], z["w2"]), t, int(steps), nxt, u_prev, list(z["ts"]), [json.loads(str(m)) for m in z["ms"]]


def run(name, N, T, ckpt=None, ckpt_every_s=120.0, series=False, budget_s=None):
    """One run; with `ckpt` (a file path) it resumes from that snapshot if present and refreshes it every ckpt_every_s.

    With `budget_s` (needs `ckpt`) the loop stops once that many wall seconds have passed in this call, saves the snapshot and
    returns {"partial": True, "t": ..., "steps": ...}; calling again resumes it. Without a budget the result is unchanged.
    """
    s = base.NS3D(N, NU)
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    ux, uy, uz = (AMP * c for c in velocity_ic(name, X, Y, Z, N))
    w = tuple(np.asarray(c) for c in s.curl(np.fft.fftn(ux), np.fft.fftn(uy), np.fft.fftn(uz)))
    t, steps, died, nxt, u_prev = 0.0, 0, None, 0.0, 6.0
    ts, ms = [], []
    if ckpt and Path(ckpt).exists():
        w, t, steps, nxt, u_prev, ts, ms = _load(ckpt)
    t0 = wall_time()
    last_save = t0
    with np.errstate(all="ignore"):
        while t < T:
            if t >= nxt:
                m = measure(s, w)
                if m is None:
                    died = round(t, 3)
                    break
                ts.append(t)
                ms.append(m)
                nxt += T / 120
            dt = min(CFL * s.dx / u_prev, T - t)
            w, u_prev = rk.ifrk4_step(s, w, dt)
            u_prev = max(u_prev, 1e-6)
            t += dt
            steps += 1
            if ckpt and wall_time() - last_save > ckpt_every_s:
                _save(ckpt, w, t, steps, nxt, u_prev, ts, ms)
                last_save = wall_time()
            if budget_s is not None and ckpt and t < T and wall_time() - t0 > budget_s:
                _save(ckpt, w, t, steps, nxt, u_prev, ts, ms)
                return {"partial": True, "ic": name, "N": N, "T": T, "t": round(t, 4), "steps": steps}
    out = {"ic": name, "N": N, "T": T, "finite": died is None, "died_at": died, "steps": steps, "seconds": round(wall_time() - t0, 1)}
    if len(ms) > 5:
        R = np.array([m["R"] for m in ms]); F = np.array([m["F"] for m in ms]); Om = np.array([m["Omega"] for m in ms])
        mw = np.array([m["max_w"] for m in ms]); tt = np.array(ts)
        dF = np.gradient(F, tt)
        denom = F * np.log(np.maximum(F, 1.0001))
        C = np.where(denom > 1e-9, dF / denom, np.nan)
        late = tt >= 0.1
        out.update({"R0": float(R[0]), "R_max": float(R.max()), "R_min": float(R.min()), "R_end": float(R[-1]),
                    "Omega_growth": float(Om.max() / Om[0]), "maxw_growth": float(mw.max() / mw[0]), "maxw_peak": float(mw.max()),
                    "F_max": float(F.max()), "C_max": float(np.nanmax(C)), "C_max_after_0.1": float(np.nanmax(C[late])) if late.any() else None,
                    "t_Omega_peak": float(tt[int(np.argmax(Om))])})
    if series and len(ms) > 5:
        out["series"] = {"t": [float(x) for x in tt], "R": R.tolist(), "F": F.tolist(), "Omega": Om.tolist(), "max_w": mw.tolist(), "C": [None if not np.isfinite(c) else float(c) for c in C]}
    return out


def run_robust(name, N, T, ckpt, tries=30, wait_s=60, series=False, budget_s=None):
    """run() that survives MemoryError (other processes exhausting memory): sleeps, then resumes from the snapshot."""
    import time
    for attempt in range(tries):
        try:
            extra = {**({"series": True} if series else {}), **({"budget_s": budget_s} if budget_s is not None else {})}
            return run(name, N, T, ckpt=ckpt, **extra)
        except MemoryError:
            print(f"MemoryError (attempt {attempt + 1}); resuming from the snapshot in {wait_s}s", flush=True)
            time.sleep(wait_s)
    raise RuntimeError(f"{name} N={N}: still out of memory after {tries} attempts; snapshot kept at {ckpt}")


def main():
    sizes = [int(a) for a in sys.argv[1:]] or [32, 48]
    results = []
    for name, T in [("tg5", 2.0), ("abc", 3.0), ("kida", 3.0), ("rand1", 3.0), ("rand2", 3.0), ("rand3", 3.0)]:
        for N in sizes:
            r = run(name, N, T)
            results.append(r)
            print(json.dumps(r), flush=True)
            (HERE / "ns_ic_study.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
