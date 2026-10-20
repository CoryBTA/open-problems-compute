import Mathlib

/-!
# Pointwise decomposition behind R = ||grad w||^2 / ||grad|w|||^2

Source: `docs/ns-step7-rigorous-proof.md` ("Pythagorean decomposition ... exact"). For a nonzero vector `v`
(the vorticity at a point) and a vector `w` (one directional derivative of the vorticity), with
`u = v / |v|`:
  * `⟨u, w⟩` is the derivative of `|v|` in that direction (the component of `w` along `v`),
  * `w - ⟨u, w⟩ u` is `|v|` times the derivative of the direction `u` in that direction.
So `|w|^2 = (d|v|)^2 + |v|^2 |d u|^2` pointwise, and in particular `|d|v|| <= |w|` (Kato). Summing over
directions and integrating gives `R = 1 + (weighted ||grad w-hat||^2) / ||grad|w|||^2 >= 1`.

This file proves only the pointwise algebra. It does not touch Navier-Stokes.
-/

open RealInnerProductSpace

namespace Pythagoras

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- Pythagoras for the splitting of `w` along the unit vector `u`. -/
theorem norm_sq_split (u w : E) (hu : ‖u‖ = 1) :
    ‖w‖ ^ 2 = ⟪u, w⟫ ^ 2 + ‖w - ⟪u, w⟫ • u‖ ^ 2 := by
  have horth : ⟪(⟪u, w⟫ • u), w - ⟪u, w⟫ • u⟫ = 0 := by
    rw [inner_smul_left, inner_sub_right, inner_smul_right, real_inner_self_eq_norm_sq, hu]
    simp
  have h := norm_add_sq_eq_norm_sq_add_norm_sq_real (x := ⟪u, w⟫ • u) (y := w - ⟪u, w⟫ • u) horth
  have hsum : ⟪u, w⟫ • u + (w - ⟪u, w⟫ • u) = w := by abel
  rw [hsum, norm_smul, hu] at h
  have h' : ‖w‖ * ‖w‖ = ⟪u, w⟫ * ⟪u, w⟫ + ‖w - ⟪u, w⟫ • u‖ * ‖w - ⟪u, w⟫ • u‖ := by
    simpa using h
  rw [sq, sq, sq]
  exact h'

/-- The same statement with `u = v / |v|` for a nonzero `v`: derivative of `|v|` plus the transverse part. -/
theorem norm_sq_split_of_ne_zero (v w : E) (hv : v ≠ 0) :
    ‖w‖ ^ 2 = ⟪‖v‖⁻¹ • v, w⟫ ^ 2 + ‖w - ⟪‖v‖⁻¹ • v, w⟫ • (‖v‖⁻¹ • v)‖ ^ 2 := by
  apply norm_sq_split
  rw [norm_smul, norm_inv, norm_norm]
  exact inv_mul_cancel₀ (norm_ne_zero_iff.mpr hv)

/-- Kato's inequality, pointwise: the derivative of the norm is at most the norm of the derivative. -/
theorem kato_pointwise (v w : E) (hv : v ≠ 0) : ⟪‖v‖⁻¹ • v, w⟫ ^ 2 ≤ ‖w‖ ^ 2 := by
  have := norm_sq_split_of_ne_zero v w hv
  nlinarith [sq_nonneg ‖w - ⟪‖v‖⁻¹ • v, w⟫ • (‖v‖⁻¹ • v)‖]

end Pythagoras

#print axioms Pythagoras.norm_sq_split
#print axioms Pythagoras.kato_pointwise
