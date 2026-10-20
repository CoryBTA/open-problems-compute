#!/usr/bin/env python3
"""NS: The FINAL gap — bound R = ||∇ω||² / ||∇|ω|||².

The proof chain:
  BG: ||ω||_∞ ≤ C · Ω^{1/2} · √(1 + log(||∇ω||²/Ω))
  Kato: ||∇|ω|||² ≤ ||∇ω||²
  Fisher: ||∇|ω|||² = FΩ/4

  Need: log(||∇ω||²/Ω) ≤ C · log(F)
  i.e.: ||∇ω||² ≤ C · Ω · F^C
  i.e.: R = ||∇ω||²/||∇|ω|||² ≤ C · F^{C-1}

  If R is BOUNDED (R ≤ R_max): then ||∇ω||² ≤ R_max · FΩ/4
  and log(||∇ω||²/Ω) ≤ log(R_max·F/4) = log(R_max) + log(F/4)
  So: ||ω||_∞ ≤ C · Ω^{1/2} · √(log R_max + log F) ≤ C' · Ω^{1/2} · √(log F)
  DONE!

  R = ||∇ω||²/||∇|ω|||² = 1 + ||ω||²·||∇ω̂||²_{L²} / ||∇|ω|||²_{L²}
  where ω̂ = ω/|ω| is the vorticity direction.

  The extra term ||∇ω̂|| measures how fast the vorticity DIRECTION changes.
  For a scalar field: R = 1 exactly.
  For a vector field: R > 1 due to directional gradients.

  R is bounded iff ∇ω̂ is controlled.
  On NS solutions: ∇ω̂ IS controlled because the vorticity equation
  constrains the rotation rate of ω̂.

  THIS SCRIPT: Measure R across all 3D configurations and determine
  if it's bounded.
"""

import json
from pathlib import Path
from time import time as wall_time

import numpy as np

OUT_DIR = Path(__file__).resolve().parent.parent / "results"


class NS3D:
    def __init__(self, N, nu):
        self.N, self.nu, self.dx = N, nu, 2*np.pi/N
        k1d = np.fft.fftfreq(N, d=1.0/N)
        self.KX,self.KY,self.KZ = np.meshgrid(k1d,k1d,k1d,indexing='ij')
        self.K2 = self.KX**2+self.KY**2+self.KZ**2
        self.K2[0,0,0] = 1
        kmax = N//3
        self.dealias = ((np.abs(self.KX)<kmax)&(np.abs(self.KY)<kmax)&
                        (np.abs(self.KZ)<kmax)).astype(float)
    def curl(self, ux_h,uy_h,uz_h):
        return (1j*(self.KY*uz_h-self.KZ*uy_h),1j*(self.KZ*ux_h-self.KX*uz_h),
                1j*(self.KX*uy_h-self.KY*ux_h))
    def vel(self, wx_h,wy_h,wz_h):
        ux=1j*(self.KY*wz_h-self.KZ*wy_h)/self.K2
        uy=1j*(self.KZ*wx_h-self.KX*wz_h)/self.K2
        uz=1j*(self.KX*wy_h-self.KY*wx_h)/self.K2
        ux[0,0,0]=0;uy[0,0,0]=0;uz[0,0,0]=0
        return ux,uy,uz
    def p(self, f_h): return np.real(np.fft.ifftn(f_h*self.dealias))
    def g(self, f_h, K): return self.p(1j*K*f_h)
    def step(self, wx_h,wy_h,wz_h, dt):
        ux_h,uy_h,uz_h = self.vel(wx_h,wy_h,wz_h)
        wx,wy,wz = self.p(wx_h),self.p(wy_h),self.p(wz_h)
        ux,uy,uz = self.p(ux_h),self.p(uy_h),self.p(uz_h)
        dux=[self.g(ux_h,K) for K in [self.KX,self.KY,self.KZ]]
        duy=[self.g(uy_h,K) for K in [self.KX,self.KY,self.KZ]]
        duz=[self.g(uz_h,K) for K in [self.KX,self.KY,self.KZ]]
        sx=wx*dux[0]+wy*dux[1]+wz*dux[2]
        sy=wx*duy[0]+wy*duy[1]+wz*duy[2]
        sz=wx*duz[0]+wy*duz[1]+wz*duz[2]
        dwx=[self.g(wx_h,K) for K in [self.KX,self.KY,self.KZ]]
        dwy=[self.g(wy_h,K) for K in [self.KX,self.KY,self.KZ]]
        dwz=[self.g(wz_h,K) for K in [self.KX,self.KY,self.KZ]]
        ax=-(ux*dwx[0]+uy*dwx[1]+uz*dwx[2])
        ay=-(ux*dwy[0]+uy*dwy[1]+uz*dwy[2])
        az=-(ux*dwz[0]+uy*dwz[1]+uz*dwz[2])
        rx=np.fft.fftn(sx+ax)*self.dealias
        ry=np.fft.fftn(sy+ay)*self.dealias
        rz=np.fft.fftn(sz+az)*self.dealias
        e=np.exp(-self.nu*self.K2*dt)
        return (wx_h+rx*dt)*e,(wy_h+ry*dt)*e,(wz_h+rz*dt)*e


def compute_R(solver, wx_h, wy_h, wz_h):
    """Compute R = ||∇ω||² / ||∇|ω|||²."""
    _N, dx = solver.N, solver.dx

    wx,wy,wz = solver.p(wx_h),solver.p(wy_h),solver.p(wz_h)
    w_mag = np.sqrt(wx**2+wy**2+wz**2+1e-30)
    Omega = float(np.sum(w_mag**2)*dx**3)
    if Omega < 1e-20:
        return None

    # ||∇ω||² = Σ_i Σ_j (∂ω_i/∂x_j)²
    grad_omega_sq = 0
    for w_h in [wx_h,wy_h,wz_h]:
        for K in [solver.KX,solver.KY,solver.KZ]:
            dw = solver.p(1j*K*w_h)
            grad_omega_sq += float(np.sum(dw**2)*dx**3)

    # ||∇|ω|||² = Σ_j (∂|ω|/∂x_j)²
    w_fft = np.fft.fftn(w_mag)
    grad_wmag_sq = 0
    for K in [solver.KX,solver.KY,solver.KZ]:
        dw = solver.p(1j*K*w_fft)
        grad_wmag_sq += float(np.sum(dw**2)*dx**3)

    R = grad_omega_sq / (grad_wmag_sq + 1e-30)

    # Fisher
    w_mag**2 / Omega
    F = 4 * grad_wmag_sq / Omega

    # Also compute ||∇ω̂||² (directional gradient)
    # ∇ω̂_i = (∂ω_i/∂x_j)/|ω| - ω_i(ω·∂ω/∂x_j)/|ω|³
    # ||∇ω̂||² = (||∇ω||² - ||∇|ω|||²) / ||ω||²_{appropriate}
    # By Pythagorean: ||∇ω||² = ||∇|ω|||² + ∫|ω|²|∇ω̂|² dx
    directional_sq = grad_omega_sq - grad_wmag_sq
    # Normalized by Ω
    directional_normalized = directional_sq / (Omega + 1e-30)

    return {
        "R": R,
        "grad_omega_sq": grad_omega_sq,
        "grad_wmag_sq": grad_wmag_sq,
        "F": F,
        "Omega": Omega,
        "directional_sq": directional_sq,
        "directional_normalized": directional_normalized,
        "max_w": float(np.max(w_mag)),
    }


def run_R_test(name, amp, nu, T, N=48):
    solver = NS3D(N, nu)
    x = np.linspace(0,2*np.pi,N,endpoint=False)
    X,Y,Z = np.meshgrid(x,x,x,indexing='ij')
    ux = amp*np.sin(X)*np.cos(Y)*np.cos(Z)
    uy = -amp*np.cos(X)*np.sin(Y)*np.cos(Z)
    uz = np.zeros_like(X)
    ux_h,uy_h,uz_h = np.fft.fftn(ux),np.fft.fftn(uy),np.fft.fftn(uz)
    wx_h,wy_h,wz_h = solver.curl(ux_h,uy_h,uz_h)

    dt = min(0.4*solver.dx/amp, 0.08*solver.dx**2/nu, 0.005)
    n_steps = int(T/dt)
    mi = max(1, n_steps//100)

    print(f"  {name}: N={N}, ν={nu}, amp={amp}, T={T}")
    t0 = wall_time()

    R_history = []
    F_history = []

    for step in range(n_steps):
        if step % mi == 0:
            data = compute_R(solver, wx_h, wy_h, wz_h)
            if data is None:
                break
            R_history.append(data["R"])
            F_history.append(data["F"])

        wx_h,wy_h,wz_h = solver.step(wx_h,wy_h,wz_h,dt)
        if np.any(np.isnan(wx_h)):
            break

    elapsed = wall_time()-t0

    if not R_history:
        return {"name":name,"R_max":0,"R_mean":0,"R_max_below_100":False}

    R_arr = np.array(R_history)
    F_arr = np.array(F_history)

    # Does R grow with F?
    if len(R_arr) > 10 and np.std(np.log(F_arr+1)) > 0.1:
        corr = float(np.corrcoef(np.log(F_arr+1), np.log(R_arr+1e-10))[0,1])
    else:
        corr = 0

    R_max = float(np.max(R_arr))
    R_mean = float(np.mean(R_arr))

    print(f"    {elapsed:.1f}s | R_max={R_max:.4f} R_mean={R_mean:.4f} "
          f"R-F corr={corr:.3f}")

    return {
        "name": name, "R_max": R_max, "R_mean": R_mean,
        "R_F_correlation": corr,
        "R_max_below_100": R_max < 100,
        "F_max": float(np.max(F_arr)),
    }


def main():
    print("="*70)
    print("NS: THE FINAL GAP — RATIO R = ||∇ω||²/||∇|ω|||²")
    print("="*70)

    results = []

    configs = [
        ("TG1_Re100", 1, 0.01, 10.0),
        ("TG5_Re100", 5, 0.05, 5.0),
        ("TG10_Re100", 10, 0.1, 3.0),
        ("TG1_Re200", 1, 0.005, 5.0),
        ("TG5_Re200", 5, 0.025, 3.0),
        ("TG1_Re500", 1, 0.002, 2.0),
        ("TG10_Re200", 10, 0.05, 2.0),
        ("TG5_Re500", 5, 0.01, 2.0),
    ]

    for name, amp, nu, T in configs:
        r = run_R_test(name, amp, nu, T)
        results.append(r)

    # Also test non-TG: Kida, shear
    N = 48
    x = np.linspace(0,2*np.pi,N,endpoint=False)
    X,Y,Z = np.meshgrid(x,x,x,indexing='ij')

    # Kida
    solver = NS3D(N, 0.01)
    ux = np.sin(X)*np.cos(Y)*np.cos(Z)+np.cos(X)*np.sin(Z)
    uy = -np.cos(X)*np.sin(Y)*np.cos(Z)+np.sin(Y)*np.cos(Z)
    uz = np.cos(Z)*np.sin(X)
    ux_h,uy_h,uz_h = np.fft.fftn(ux),np.fft.fftn(uy),np.fft.fftn(uz)
    wx_h,wy_h,wz_h = solver.curl(ux_h,uy_h,uz_h)
    dt = 0.003
    R_vals = []
    for step in range(1000):
        if step % 20 == 0:
            d = compute_R(solver, wx_h,wy_h,wz_h)
            if d: R_vals.append(d["R"])
        wx_h,wy_h,wz_h = solver.step(wx_h,wy_h,wz_h,dt)
        if np.any(np.isnan(wx_h)): break
    if R_vals:
        results.append({"name":"Kida_Re100","R_max":float(max(R_vals)),
                        "R_mean":float(np.mean(R_vals)),"R_max_below_100":max(R_vals)<100,
                        "R_F_correlation":0,"F_max":0})
        print(f"  Kida_Re100: R_max={max(R_vals):.4f} R_mean={np.mean(R_vals):.4f}")

    # Summary
    print("\n"+"="*70)
    print("R = ||∇ω||²/||∇|ω|||² — THE FINAL RATIO")
    print("="*70)

    print(f"\n  {'Config':<18} {'R_max':>8} {'R_mean':>8} {'R-F corr':>8} {'OK':>4}")
    print("  "+"-"*48)
    for r in results:
        ok = "✓" if r["R_max_below_100"] else "✗"
        print(f"  {r['name']:<18} {r['R_max']:>8.4f} {r['R_mean']:>8.4f} "
              f"{r.get('R_F_correlation',0):>8.3f} {ok:>4}")

    R_max_all = max(r["R_max"] for r in results)
    all_below_100 = all(r["R_max_below_100"] for r in results)
    R_F_corrs = [r.get("R_F_correlation",0) for r in results]

    print(f"\n  Max R across all configs: {R_max_all:.4f}")
    print(f"  All R_max below 100: {all_below_100}")
    print(f"  R-F correlation range: [{min(R_F_corrs):.3f}, {max(R_F_corrs):.3f}]")

    print("  Measured at one grid size (N=48) only. This is not a bound on R: see resolution_sweep.py,")
    print("  which shows R_max changing with N for the amp=5, nu=0.01 configuration.")

    output = {"results":results,"R_max":R_max_all,"all_below_100":all_below_100}
    out_path = OUT_DIR/"ns-ratio-R-results.json"
    with open(out_path,"w") as f:
        json.dump(output,f,indent=2,default=str)
    print(f"\n  Results: {out_path}")


if __name__ == "__main__":
    main()
