import Mathlib

/-!
# NS regularity chain: the closing ODE step (Step 6 consequence)

Source: `docs/ns-step7-rigorous-proof.md`. The chain there is
  Step 6: dF/dt <= C * F * log F   (numerical, NOT proved analytically)
  Step 7: R <= R_*                 (numerical, NOT proved)
  => ||omega||_inf controlled => Beale-Kato-Majda => regularity.

What this file proves: IF the differential inequality of Step 6 holds, THEN F cannot blow up
in finite time, with the explicit double-exponential bound
  log F(t) <= log F(0) * exp (C * t).
It does not prove Step 6 itself. Hypothesis `hF` below is exactly the unproved input.
-/

open Real Set

namespace NSChain

/-- Log-Gronwall closure: `F' <= C F log F` and `F > 1` give `log F t <= log F 0 * exp (C t)`. -/
theorem log_gronwall (F : ℝ → ℝ) (C T : ℝ) (hT : 0 ≤ T)
    (hdiff : ∀ t ∈ Icc (0 : ℝ) T, HasDerivAt F (deriv F t) t)
    (hpos : ∀ t ∈ Icc (0 : ℝ) T, 1 < F t)
    (hF : ∀ t ∈ Icc (0 : ℝ) T, deriv F t ≤ C * F t * Real.log (F t)) :
    ∀ t ∈ Icc (0 : ℝ) T, Real.log (F t) ≤ Real.log (F 0) * Real.exp (C * t) := by
  -- G t = log (F t) * exp (-C t) is antitone on [0, T].
  set G : ℝ → ℝ := fun t => Real.log (F t) * Real.exp (-C * t) with hG
  have hGderiv : ∀ t ∈ Icc (0 : ℝ) T,
      HasDerivAt G ((deriv F t / F t) * Real.exp (-C * t)
        + Real.log (F t) * (Real.exp (-C * t) * (-C))) t := by
    intro t ht
    have h1 : HasDerivAt (fun s => Real.log (F s)) (deriv F t / F t) t :=
      (hdiff t ht).log (ne_of_gt (lt_trans zero_lt_one (hpos t ht)))
    have h2 : HasDerivAt (fun s => Real.exp (-C * s)) (Real.exp (-C * t) * (-C)) t := by
      have := ((hasDerivAt_id t).const_mul (-C)).exp
      simpa using this
    exact h1.mul h2
  have hnonpos : ∀ t ∈ Icc (0 : ℝ) T,
      (deriv F t / F t) * Real.exp (-C * t)
        + Real.log (F t) * (Real.exp (-C * t) * (-C)) ≤ 0 := by
    intro t ht
    have hFt : 0 < F t := lt_trans zero_lt_one (hpos t ht)
    have h := hF t ht
    have hdiv : deriv F t / F t ≤ C * Real.log (F t) := by
      rw [div_le_iff₀ hFt]; nlinarith
    have hexp := Real.exp_pos (-C * t)
    nlinarith
  have hanti : AntitoneOn G (Icc (0 : ℝ) T) := by
    apply antitoneOn_of_deriv_nonpos (convex_Icc 0 T)
    · exact fun t ht => (hGderiv t ht).continuousAt.continuousWithinAt
    · intro t ht
      exact (hGderiv t (interior_subset ht)).differentiableAt.differentiableWithinAt
    · intro t ht
      rw [(hGderiv t (interior_subset ht)).deriv]
      exact hnonpos t (interior_subset ht)
  intro t ht
  have h0 : (0 : ℝ) ∈ Icc (0 : ℝ) T := ⟨le_refl 0, hT⟩
  have := hanti h0 ht ht.1
  simp only [hG, mul_zero, Real.exp_zero, mul_one] at this
  have hexp := Real.exp_pos (C * t)
  calc Real.log (F t) = (Real.log (F t) * Real.exp (-C * t)) * Real.exp (C * t) := by
        rw [mul_assoc, ← Real.exp_add]; simp
    _ ≤ Real.log (F 0) * Real.exp (C * t) := by
        exact mul_le_mul_of_nonneg_right this hexp.le

/-- Corollary: `F` stays finite on every `[0, T]` (no finite-time blowup), under Step 6's inequality. -/
theorem no_finite_time_blowup (F : ℝ → ℝ) (C T : ℝ) (hT : 0 ≤ T)
    (hdiff : ∀ t ∈ Icc (0 : ℝ) T, HasDerivAt F (deriv F t) t)
    (hpos : ∀ t ∈ Icc (0 : ℝ) T, 1 < F t)
    (hF : ∀ t ∈ Icc (0 : ℝ) T, deriv F t ≤ C * F t * Real.log (F t)) :
    ∃ M : ℝ, ∀ t ∈ Icc (0 : ℝ) T, F t ≤ M := by
  refine ⟨Real.exp (Real.log (F 0) * Real.exp (|C| * T)), fun t ht => ?_⟩
  have h := log_gronwall F C T hT hdiff hpos hF t ht
  have hFt : 0 < F t := lt_trans zero_lt_one (hpos t ht)
  have hle : Real.exp (C * t) ≤ Real.exp (|C| * T) := by
    apply Real.exp_le_exp.mpr
    calc C * t ≤ |C| * t := mul_le_mul_of_nonneg_right (le_abs_self C) ht.1
      _ ≤ |C| * T := mul_le_mul_of_nonneg_left ht.2 (abs_nonneg C)
  have hlog0 : 0 ≤ Real.log (F 0) := Real.log_nonneg (hpos 0 ⟨le_refl 0, hT⟩).le
  rw [← Real.exp_log hFt]
  apply Real.exp_le_exp.mpr
  exact h.trans (mul_le_mul_of_nonneg_left hle hlog0)

end NSChain
