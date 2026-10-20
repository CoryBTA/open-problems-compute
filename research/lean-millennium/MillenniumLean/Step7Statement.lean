import Mathlib
import MillenniumLean.NSChain

/-!
# Step 7 restated: an affine bound on G = ||grad w||^2 / Omega instead of a bound on R

`NSHypotheses.lean` assumes `R <= R*` with `R = ||grad w||^2 / ||grad|w|||^2`. That hypothesis is false as stated: the
Roberts flow `u = (sin z, cos z, 0)` is an exact, smooth, steady Euler solution (a decaying one for Navier-Stokes) with
`|w| = 1` constant, so `grad|w| = 0` and `R` is infinite, while `||grad w||^2 / Omega = 1`.

What the Brezis-Gallouet step actually uses is `G = ||grad w||^2 / Omega` (note `G = R * F / 4` where defined, with
`F = 4 ||grad|w|||^2 / Omega`). The weaker, correct hypothesis is the affine bound
  `0 < G t <= A * F t`   (from `G <= (R*/4) F + K` and `F > 1`, with `A = R*/4 + K`),
which holds for the Roberts flow (`G = 1`, `F` tiny, `K = 1`). This file proves that, with the same Step 6 input, the
chain still gives a bound on `||w||_inf`. Neither Step 6 nor the affine bound is proved for Navier-Stokes.
-/

open Real Set

namespace NSChain

theorem chain_affine
    (Ω F G W : ℝ → ℝ) (Emax C C1 A : ℝ) (hC1 : 0 ≤ C1) (hA : 0 < A) (T : ℝ) (hT : 0 ≤ T)
    -- energy-type bound on the enstrophy (known; not in Mathlib)
    (hEnergy : ∀ t ∈ Icc (0 : ℝ) T, 0 ≤ Ω t ∧ Ω t ≤ Emax)
    -- Brezis-Gallouet-type bound in terms of G (known inequality; not in Mathlib)
    (hBG : ∀ t ∈ Icc (0 : ℝ) T, 0 < G t →
      W t ≤ C1 * Real.sqrt (Ω t) * Real.sqrt (1 + Real.log (G t)))
    -- Step 7, corrected (NOT proved): affine bound of G by F
    (hG : ∀ t ∈ Icc (0 : ℝ) T, 0 < G t ∧ G t ≤ A * F t)
    -- Step 6 (NUMERICAL ONLY), with F differentiable and above 1
    (hdiff : ∀ t ∈ Icc (0 : ℝ) T, HasDerivAt F (deriv F t) t)
    (hpos : ∀ t ∈ Icc (0 : ℝ) T, 1 < F t)
    (hStep6 : ∀ t ∈ Icc (0 : ℝ) T, deriv F t ≤ C * F t * Real.log (F t)) :
    ∃ M : ℝ, ∀ t ∈ Icc (0 : ℝ) T, W t ≤ M := by
  obtain ⟨MF, hMF⟩ := no_finite_time_blowup F C T hT hdiff hpos hStep6
  refine ⟨C1 * Real.sqrt Emax * Real.sqrt (1 + Real.log (A * MF)), fun t ht => ?_⟩
  obtain ⟨hGpos, hGle⟩ := hG t ht
  have hFM : F t ≤ MF := hMF t ht
  have hGM : G t ≤ A * MF := le_trans hGle (mul_le_mul_of_nonneg_left hFM hA.le)
  have hlog : Real.log (G t) ≤ Real.log (A * MF) := Real.log_le_log hGpos hGM
  have h1 : Real.sqrt (Ω t) ≤ Real.sqrt Emax := Real.sqrt_le_sqrt (hEnergy t ht).2
  have h2 : Real.sqrt (1 + Real.log (G t)) ≤ Real.sqrt (1 + Real.log (A * MF)) :=
    Real.sqrt_le_sqrt (by linarith)
  calc W t ≤ C1 * Real.sqrt (Ω t) * Real.sqrt (1 + Real.log (G t)) := hBG t ht hGpos
    _ ≤ C1 * Real.sqrt Emax * Real.sqrt (1 + Real.log (A * MF)) := by
        apply mul_le_mul
        · exact mul_le_mul_of_nonneg_left h1 hC1
        · exact h2
        · exact Real.sqrt_nonneg _
        · exact mul_nonneg hC1 (Real.sqrt_nonneg _)

/-- The old hypothesis implies the new one: `R <= R*` gives `G = R F / 4 <= (R*/4) F` (for `F >= 0`). -/
theorem affine_of_ratio_bound (R F Rstar : ℝ) (hF : 0 ≤ F) (hR : R ≤ Rstar) :
    R * F / 4 ≤ (Rstar / 4) * F := by
  nlinarith

end NSChain

#print axioms NSChain.chain_affine
#print axioms NSChain.affine_of_ratio_bound
