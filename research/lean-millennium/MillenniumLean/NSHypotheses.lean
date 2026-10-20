import Mathlib
import MillenniumLean.NSChain

/-!
# NS regularity chain as explicit hypotheses

Abstracts the chain in `docs/ns-step7-rigorous-proof.md` to functions of time:
`Ω` enstrophy, `F` Fisher information, `R` the R-ratio, `W` = ||ω||_∞.

Every input that is NOT proved anywhere is a named hypothesis of `chain_gives_bounded_vorticity`:
* `hEnergy`  (Step 1, known: energy dissipation bound on enstrophy),
* `hBG`      (Brezis-Gallouet-type bound on ||ω||_∞; known theorem, not in Mathlib),
* `hR`       (Step 7: R <= R*, NUMERICAL ONLY),
* `hStep6`   (Step 6: dF/dt <= C F log F, NUMERICAL ONLY).
The conclusion is an explicit bound on ||ω||_∞ over [0,T]; Beale-Kato-Majda (not in Mathlib) then
turns that bound into regularity, and is stated as `BKM` below instead of being assumed silently.

This file proves the implication only. It does NOT prove any of the hypotheses.
-/

open Real Set

namespace NSChain

/-- Beale-Kato-Majda, stated abstractly: a bounded vorticity sup-norm on every `[0,T]`
rules out blowup. A known theorem, not formalized here. -/
def BKM (W : ℝ → ℝ) (Regular : Prop) : Prop :=
  (∀ T : ℝ, 0 ≤ T → ∃ M : ℝ, ∀ t ∈ Icc (0 : ℝ) T, W t ≤ M) → Regular

theorem chain_gives_bounded_vorticity
    (Ω F R W : ℝ → ℝ) (E0 ν C C1 Rstar : ℝ) (hC1 : 0 ≤ C1)
    (hRstar : 0 < Rstar)
    -- Step 1 (known): enstrophy bound
    (hEnergy : ∀ t, 0 ≤ t → 0 ≤ Ω t ∧ Ω t ≤ E0 / ν)
    -- Brezis-Gallouet-type bound (known theorem, not in Mathlib)
    (hBG : ∀ t, 0 ≤ t → 1 < F t → 0 < R t →
      W t ≤ C1 * Real.sqrt (Ω t) * Real.sqrt (1 + Real.log (R t * F t / 4)))
    -- Step 7 (NUMERICAL ONLY)
    (hR : ∀ t, 0 ≤ t → 0 < R t ∧ R t ≤ Rstar)
    -- Step 6 (NUMERICAL ONLY), with F differentiable and above 1
    (T : ℝ) (hT : 0 ≤ T)
    (hdiff : ∀ t ∈ Icc (0 : ℝ) T, HasDerivAt F (deriv F t) t)
    (hpos : ∀ t ∈ Icc (0 : ℝ) T, 1 < F t)
    (hStep6 : ∀ t ∈ Icc (0 : ℝ) T, deriv F t ≤ C * F t * Real.log (F t)) :
    ∃ M : ℝ, ∀ t ∈ Icc (0 : ℝ) T, W t ≤ M := by
  obtain ⟨MF, hMF⟩ := no_finite_time_blowup F C T hT hdiff hpos hStep6
  refine ⟨C1 * Real.sqrt (E0 / ν) * Real.sqrt (1 + Real.log (Rstar * MF / 4)), fun t ht => ?_⟩
  have ht0 : 0 ≤ t := ht.1
  have hFt := hpos t ht
  obtain ⟨hRpos, hRle⟩ := hR t ht0
  have hbg := hBG t ht0 hFt hRpos
  have hΩ := hEnergy t ht0
  have hFM : F t ≤ MF := hMF t ht
  have hprod : R t * F t / 4 ≤ Rstar * MF / 4 := by
    have h1 : R t * F t ≤ Rstar * MF :=
      mul_le_mul hRle hFM (by linarith) hRstar.le
    linarith
  have hprodpos : 0 < R t * F t / 4 := by
    have : 0 < F t := lt_trans zero_lt_one hFt
    positivity
  have hlog : Real.log (R t * F t / 4) ≤ Real.log (Rstar * MF / 4) :=
    Real.log_le_log hprodpos hprod
  have hsqrt1 : Real.sqrt (Ω t) ≤ Real.sqrt (E0 / ν) := Real.sqrt_le_sqrt hΩ.2
  have hsqrt2 : Real.sqrt (1 + Real.log (R t * F t / 4))
      ≤ Real.sqrt (1 + Real.log (Rstar * MF / 4)) := Real.sqrt_le_sqrt (by linarith)
  calc W t ≤ C1 * Real.sqrt (Ω t) * Real.sqrt (1 + Real.log (R t * F t / 4)) := hbg
    _ ≤ C1 * Real.sqrt (E0 / ν) * Real.sqrt (1 + Real.log (Rstar * MF / 4)) := by
        apply mul_le_mul
        · exact mul_le_mul_of_nonneg_left hsqrt1 hC1
        · exact hsqrt2
        · exact Real.sqrt_nonneg _
        · exact mul_nonneg hC1 (Real.sqrt_nonneg _)

end NSChain

#print axioms NSChain.chain_gives_bounded_vorticity
