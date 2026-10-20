"""For each N: did the TG5 nu=0.01 run stay finite to T=2, and what did max|omega| and enstrophy do?"""
import importlib.util, sys
import numpy as np
from pathlib import Path
spec = importlib.util.spec_from_file_location("m", Path(__file__).with_name("ns_ratio_R.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
for N in [int(a) for a in sys.argv[1:]]:
    amp, nu, T = 5, 0.01, 2.0
    s = m.NS3D(N, nu)
    x = np.linspace(0, 2*np.pi, N, endpoint=False); X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    ux = amp*np.sin(X)*np.cos(Y)*np.cos(Z); uy = -amp*np.cos(X)*np.sin(Y)*np.cos(Z); uz = np.zeros_like(X)
    wh = s.curl(np.fft.fftn(ux), np.fft.fftn(uy), np.fft.fftn(uz))
    dt = min(0.4*s.dx/amp, 0.08*s.dx**2/nu, 0.005); n = int(T/dt)
    died = None; peak = 0.0
    with np.errstate(all='ignore'):
        for i in range(n):
            wh = s.step(*wh, dt)
            if i % 20 == 0:
                d = m.compute_R(s, *wh)
                if d is None or not np.isfinite(d["max_w"]) or d["max_w"] > 1e6:
                    died = round(i*dt, 3); break
                peak = max(peak, d["max_w"])
    print(f"N={N} dt={dt:.5f} steps={n} finite_to_T={died is None} died_at_t={died} peak_max_omega={peak:.1f}", flush=True)
