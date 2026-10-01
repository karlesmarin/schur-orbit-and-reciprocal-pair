/- ParidadNucleo.lean — two finite lemmas of `PRUEBA_VUELTAS.md`:
   (1) "Lema P" (parity): an even function `ℤ → K` whose Mahler coefficients vanish at every even
       index is zero (char K ≠ 2, in the form `2a = 0 → a = 0`);
   (2) step 7 of "Teorema F": on `ZMod (p^b) → K`, `K` a field of characteristic `p`, the kernel of
       `Δ^r` (`Δ f x = f (x+1) - f x`) is spanned by the binomial functions `C(x.val, i)`, `i < r`.
   Author: Carles Marín  <karlesmarin@gmail.com>   (with Claude, Anthropic, as assistant)

   Main theorems: `mahlerZ_coeff_eq_zero`, `mahlerZ_coeff_eq_zero_field` (Lema P),
   `ker_diffₗ_pow`, `ker_diffₗ_pow'` (kernel of Δ^r).
   Fast-loop build: lake env lean <this file>  (from E:\proyectos\godsil-gutman-lean). -/
import Mathlib.Data.ZMod.Basic
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Algebra.Group.ForwardDiff
import Mathlib.RingTheory.Binomial
import Mathlib.Data.Nat.Multiplicity
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Algebra.CharP.Defs
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination

namespace ConjeturaC

open Finset fwdDiff Module

/-! ## Lema P (parity)

Mapping to the proof: a row `J < p^a` is the full Mahler vector `(c_i)_{i ≤ D}` of an even
function `G(x) = ∑_{i ≤ D} c_i C(x, i)` on `ℤ` (values in `K = 𝔽_p`, `p` odd). Here `C(x, i)` for
`x : ℤ` is Mathlib's `Ring.choose x i` (the integer-valued binomial, `C(-x,i) = (-1)^i C(x+i-1,i)`),
cast to `K`. This is "Version A" of the task restricted to integer arguments (which is all the
lemma uses); `K` is any commutative ring in which `2a = 0 → a = 0`.

Proof (differs from the informal descending recursion `2 c_k = ∑_{i>k} …`, same content):
`Δ^D G` is the constant `c_D`; by Newton's formula and evenness, `(Δ^D G)(-D) = (-1)^D (Δ^D G)(0)`,
so for odd `D`, `c_D = -c_D`, hence `c_D = 0`; for even `D`, `c_D = 0` by hypothesis; induct on `D`. -/

section Parity

variable {K : Type*} [CommRing K]

/-- `mahlerZ D c x = ∑_{i ≤ D} c_i · C(x, i)` for `x : ℤ`. -/
def mahlerZ (D : ℕ) (c : ℕ → K) (x : ℤ) : K :=
  ∑ i ∈ range (D + 1), c i * ((Ring.choose x i : ℤ) : K)

lemma fwdDiff_ringChoose_succ (j : ℕ) :
    Δ_[(1 : ℤ)] (fun x : ℤ => ((Ring.choose x (j + 1) : ℤ) : K)) =
      fun x => ((Ring.choose x j : ℤ) : K) := by
  ext x
  simp only [fwdDiff, Ring.choose_succ_succ, Int.cast_add, add_sub_cancel_right]

lemma fwdDiff_iter_ringChoose (n j : ℕ) :
    Δ_[(1 : ℤ)]^[n] (fun x : ℤ => ((Ring.choose x (n + j) : ℤ) : K)) =
      fun x => ((Ring.choose x j : ℤ) : K) := by
  induction n generalizing j with
  | zero => simp
  | succ n ih =>
    rw [Function.iterate_succ_apply, show n + 1 + j = (n + j) + 1 by omega,
      fwdDiff_ringChoose_succ, ih]

lemma fwdDiff_iter_zero {M G : Type*} [AddCommMonoid M] [AddCommGroup G] (h : M) (k : ℕ) :
    Δ_[h]^[k] (0 : M → G) = 0 := by
  induction k with
  | zero => rfl
  | succ k ih => rw [Function.iterate_succ_apply']; rw [ih]; ext; simp [fwdDiff]

lemma fwdDiff_iter_ringChoose_lt {n m : ℕ} (h : m < n) :
    Δ_[(1 : ℤ)]^[n] (fun x : ℤ => ((Ring.choose x m : ℤ) : K)) = 0 := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_lt h
  have h1 := fwdDiff_iter_ringChoose (K := K) m 0
  simp only [add_zero, Ring.choose_zero_right, Int.cast_one] at h1
  rw [show m + k + 1 = (k + 1) + m by omega, Function.iterate_add_apply, h1,
    Function.iterate_succ_apply]
  have : Δ_[(1 : ℤ)] (fun _ : ℤ => (1 : K)) = 0 := by ext; simp [fwdDiff]
  rw [this, fwdDiff_iter_zero]

/-- `Δ^D` of `mahlerZ D c` is the constant `c D`. -/
lemma fwdDiff_iter_mahlerZ (D : ℕ) (c : ℕ → K) (y : ℤ) :
    Δ_[(1 : ℤ)]^[D] (mahlerZ D c) y = c D := by
  have hG : mahlerZ D c =
      ∑ i ∈ range (D + 1), c i • (fun x : ℤ => ((Ring.choose x i : ℤ) : K)) := by
    ext x; simp [mahlerZ, Finset.sum_apply, smul_eq_mul]
  rw [hG, fwdDiff_iter_finsetSum, Finset.sum_apply, Finset.sum_range_succ, Finset.sum_eq_zero]
  · have h1 := fwdDiff_iter_ringChoose (K := K) D 0
    simp only [add_zero, Ring.choose_zero_right, Int.cast_one] at h1
    simp [fwdDiff_iter_const_smul, h1]
  · intro i hi
    rw [fwdDiff_iter_const_smul, fwdDiff_iter_ringChoose_lt (Finset.mem_range.mp hi)]
    simp

/-- Newton + evenness: `(Δ^D G)(-D) = (-1)^D (Δ^D G)(0)` for even `G`. -/
lemma fwdDiff_iter_even (G : ℤ → K) (heven : ∀ x, G (-x) = G x) (D : ℕ) :
    Δ_[(1 : ℤ)]^[D] G (-(D : ℤ)) = (-1) ^ D * Δ_[(1 : ℤ)]^[D] G 0 := by
  rw [fwdDiff_iter_eq_sum_shift, fwdDiff_iter_eq_sum_shift, Finset.mul_sum]
  -- reflect `k ↦ D - k` on the left
  have hL : ∀ k ∈ range (D + 1),
      ((-1 : ℤ) ^ (D - k) * (D.choose k : ℤ)) • G (-(D : ℤ) + (k : ℕ) • (1 : ℤ)) =
      (fun j : ℕ => ((-1 : ℤ) ^ j * (D.choose j : ℤ)) • G (j : ℤ)) (D + 1 - 1 - k) := by
    intro k hk
    have hkD : k ≤ D := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
    simp only [Nat.add_sub_cancel]
    rw [Nat.choose_symm hkD, Nat.cast_sub hkD, ← heven, nsmul_eq_mul, mul_one]
    congr 2
    ring
  rw [Finset.sum_congr rfl hL, Finset.sum_range_reflect (fun j : ℕ =>
      ((-1 : ℤ) ^ j * (D.choose j : ℤ)) • G (j : ℤ)) (D + 1)]
  refine Finset.sum_congr rfl fun k hk => ?_
  have hkD : k ≤ D := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
  have hsign : ((-1 : K) ^ D * (-1) ^ (D - k)) = (-1) ^ k := by
    rw [← pow_add, show D + (D - k) = k + 2 * (D - k) by omega, pow_add, pow_mul]
    simp
  simp only [zero_add, nsmul_eq_mul, mul_one, zsmul_eq_mul, Int.cast_mul, Int.cast_pow,
    Int.cast_neg, Int.cast_one, Int.cast_natCast]
  rw [← mul_assoc, ← mul_assoc, hsign]

/-- The top coefficient of an even `mahlerZ` vanishes when `D` is odd. -/
lemma mahlerZ_top_odd (h2 : ∀ a : K, 2 * a = 0 → a = 0) (D : ℕ) (c : ℕ → K)
    (heven : ∀ x, mahlerZ D c (-x) = mahlerZ D c x) (hD : Odd D) : c D = 0 := by
  have h := fwdDiff_iter_even (mahlerZ D c) heven D
  rw [fwdDiff_iter_mahlerZ, fwdDiff_iter_mahlerZ, hD.neg_one_pow] at h
  apply h2
  rw [two_mul]
  nth_rewrite 1 [h]
  ring

/-- **Lema P.** If `G(x) = ∑_{i ≤ D} c_i C(x,i)` is even on `ℤ` and `c_i = 0` for every even `i`,
then `c_i = 0` for all `i ≤ D` (hence `G = 0`). `K`: `2a = 0 → a = 0`. -/
theorem mahlerZ_coeff_eq_zero (h2 : ∀ a : K, 2 * a = 0 → a = 0) (D : ℕ) (c : ℕ → K)
    (heven : ∀ x : ℤ, mahlerZ D c (-x) = mahlerZ D c x) (hc : ∀ i, Even i → c i = 0) :
    ∀ i ≤ D, c i = 0 := by
  induction D with
  | zero => intro i hi; exact hc i (by rw [Nat.le_zero.mp hi]; exact ⟨0, rfl⟩)
  | succ D ih =>
    have htop : c (D + 1) = 0 := by
      rcases Nat.even_or_odd (D + 1) with he | ho
      · exact hc _ he
      · exact mahlerZ_top_odd h2 (D + 1) c heven ho
    have hred : mahlerZ (D + 1) c = mahlerZ D c := by
      ext x; rw [mahlerZ, Finset.sum_range_succ, htop, zero_mul, add_zero, mahlerZ]
    rw [hred] at heven
    intro i hi
    rcases Nat.lt_or_eq_of_le hi with hlt | rfl
    · exact ih heven i (Nat.lt_succ_iff.mp hlt)
    · exact htop

/-- Field form of Lema P (`(2 : K) ≠ 0`, e.g. `K = ZMod p`, `p` odd). -/
theorem mahlerZ_coeff_eq_zero_field {F : Type*} [Field F] (h2 : (2 : F) ≠ 0) (D : ℕ)
    (c : ℕ → F) (heven : ∀ x : ℤ, mahlerZ D c (-x) = mahlerZ D c x)
    (hc : ∀ i, Even i → c i = 0) : ∀ i ≤ D, c i = 0 :=
  mahlerZ_coeff_eq_zero (fun _ ha => (mul_eq_zero.mp ha).resolve_left h2) D c heven hc

end Parity

/-! ## Kernel of `Δ^r` on `ZMod (p^b) → K`, `char K = p` (step 7 of Teorema F) -/

/-- Copied verbatim from `MahlerModP.lean` (same namespace; kept local because the fast env
compiles single files and cannot import sibling files). -/
def chooseFun (N : ℕ) (R : Type*) [CommRing R] (i : Fin N) : ZMod N → R :=
  fun x => ((x.val.choose i : ℕ) : R)

/-- Copied verbatim from `MahlerModP.lean`. -/
theorem linearIndependent_choose (N : ℕ) [NeZero N] (R : Type*) [CommRing R] [Nontrivial R] :
    LinearIndependent R (chooseFun N R) := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  have key : ∀ n : ℕ, ∀ i : Fin N, i.val = n → g i = 0 := by
    intro n
    induction n using Nat.strong_induction_on with
    | _ n ih =>
      intro i hi
      have h := congrFun hg (i.val : ZMod N)
      simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply, chooseFun,
        ZMod.val_natCast_of_lt i.isLt] at h
      rw [Finset.sum_eq_single i] at h
      · simpa using h
      · intro j _ hji
        rcases lt_or_gt_of_ne (Fin.val_ne_of_ne hji) with hlt | hgt
        · rw [ih j.val (by omega) j rfl, zero_mul]
        · rw [Nat.choose_eq_zero_of_lt hgt, Nat.cast_zero, mul_zero]
      · simp
  exact fun i => key i.val i rfl

section Kernel

variable (N : ℕ) (K : Type*) [Field K]

/-- `Δ f x = f (x + 1) - f x` as a `K`-linear endomorphism of `ZMod N → K`. -/
def diffₗ : Module.End K (ZMod N → K) where
  toFun := fwdDiff (1 : ZMod N)
  map_add' f g := fwdDiff_add _ f g
  map_smul' a f := fwdDiff_const_smul _ a f

lemma diffₗ_pow_apply (r : ℕ) (f : ZMod N → K) :
    (diffₗ N K ^ r) f = Δ_[(1 : ZMod N)]^[r] f := by
  rw [Module.End.pow_apply]; rfl

/-- `ℕ`-indexed binomial function `x ↦ C(x.val, i)`. -/
def binF (i : ℕ) : ZMod N → K := fun x => ((x.val.choose i : ℕ) : K)

lemma chooseFun_eq_binF (i : Fin N) : chooseFun N K i = binF N K i := rfl

variable {N K}
variable {p : ℕ} [hp : Fact p.Prime] [CharP K p]

/-- Pascal with wrap-around: `Δ C(·, i+1) = C(·, i)` on `ZMod (p^b)`, for `i + 1 ≠ p^b`
(the wrap `x = p^b - 1 ↦ 0` uses `p ∣ C(p^b, i+1)`). -/
lemma fwdDiff_binF_succ (b : ℕ) [NeZero (p ^ b)] (i : ℕ) (hi : i + 1 ≠ p ^ b) :
    Δ_[(1 : ZMod (p ^ b))] (binF (p ^ b) K (i + 1)) = binF (p ^ b) K i := by
  ext x
  have hN : 0 < p ^ b := Nat.pos_of_ne_zero (NeZero.ne _)
  have hx : x.val < p ^ b := ZMod.val_lt x
  have hx1 : (x + 1 : ZMod (p ^ b)) = ((x.val + 1 : ℕ) : ZMod (p ^ b)) := by
    rw [Nat.cast_add, ZMod.natCast_zmod_val, Nat.cast_one]
  simp only [fwdDiff, binF, hx1, ZMod.val_natCast]
  have hx' : x.val + 1 ≤ p ^ b := hx
  rcases Nat.lt_or_eq_of_le hx' with hlt | heq
  · rw [Nat.mod_eq_of_lt hlt, Nat.choose_succ_succ', Nat.cast_add, add_sub_cancel_right]
  · rw [heq, Nat.mod_self, Nat.choose_zero_succ, Nat.cast_zero, zero_sub]
    have hdvd : p ∣ (p ^ b).choose (i + 1) :=
      hp.out.dvd_choose_pow (Nat.succ_ne_zero i) hi
    have h0 : (((p ^ b).choose (i + 1) : ℕ) : K) = 0 := (CharP.cast_eq_zero_iff K p _).mpr hdvd
    rw [← heq, Nat.choose_succ_succ', Nat.cast_add] at h0
    linear_combination (-1 : K) * h0

lemma fwdDiff_binF_zero (N : ℕ) :
    Δ_[(1 : ZMod N)] (binF N K 0) = 0 := by
  ext x; simp [fwdDiff, binF]

lemma fwdDiff_iter_binF (b : ℕ) [NeZero (p ^ b)] (n j : ℕ) (h : n + j < p ^ b) :
    Δ_[(1 : ZMod (p ^ b))]^[n] (binF (p ^ b) K (n + j)) = binF (p ^ b) K j := by
  induction n generalizing j with
  | zero => simp
  | succ n ih =>
    rw [Function.iterate_succ_apply, show n + 1 + j = (n + j) + 1 by omega,
      fwdDiff_binF_succ b (n + j) (by omega), ih j (by omega)]

lemma fwdDiff_iter_binF_lt (b : ℕ) [NeZero (p ^ b)] {n m : ℕ} (h : m < n)
    (hm : m < p ^ b) : Δ_[(1 : ZMod (p ^ b))]^[n] (binF (p ^ b) K m) = 0 := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_lt h
  have h1 := fwdDiff_iter_binF (K := K) b m 0 (by omega)
  rw [add_zero] at h1
  rw [show m + k + 1 = (k + 1) + m by omega, Function.iterate_add_apply, h1,
    Function.iterate_succ_apply, fwdDiff_binF_zero, fwdDiff_iter_zero]

/-- **Kernel of `Δ^r`** (step 7 of Teorema F). `N = p^b`, `K` a field of characteristic `p`,
`r ≤ N`: `ker Δ^r = span {C(x.val, i) : i < r}`. -/
theorem ker_diffₗ_pow (b : ℕ) [NeZero (p ^ b)] (r : ℕ) (hr : r ≤ p ^ b) :
    LinearMap.ker (diffₗ (p ^ b) K ^ r) =
      Submodule.span K (Set.range fun i : Fin r => chooseFun (p ^ b) K (Fin.castLE hr i)) := by
  set N := p ^ b with hNdef
  have hli := linearIndependent_choose N K
  -- span ≤ ker
  have hle : Submodule.span K (Set.range fun i : Fin r => chooseFun N K (Fin.castLE hr i)) ≤
      LinearMap.ker (diffₗ N K ^ r) := by
    rw [Submodule.span_le]
    rintro _ ⟨i, rfl⟩
    simp only [SetLike.mem_coe, LinearMap.mem_ker, diffₗ_pow_apply, chooseFun_eq_binF]
    exact fwdDiff_iter_binF_lt b i.isLt (by simp only [Fin.val_castLE]; omega)
  refine (Submodule.eq_of_le_of_finrank_le hle ?_).symm
  -- finrank ker ≤ r
  have hspan : finrank K (Submodule.span K
      (Set.range fun i : Fin r => chooseFun N K (Fin.castLE hr i))) = r :=
    (finrank_span_eq_card (hli.comp _ (Fin.castLE_injective hr))).trans (Fintype.card_fin r)
  rw [hspan]
  have hrange : Submodule.span K
      (Set.range fun j : Fin (N - r) => chooseFun N K (Fin.castLE (Nat.sub_le N r) j)) ≤
      LinearMap.range (diffₗ N K ^ r) := by
    rw [Submodule.span_le]
    rintro _ ⟨j, rfl⟩
    refine ⟨binF N K (r + j), ?_⟩
    rw [diffₗ_pow_apply]
    exact fwdDiff_iter_binF b r j (by have := j.isLt; omega)
  have h2 : finrank K (Submodule.span K
      (Set.range fun j : Fin (N - r) => chooseFun N K (Fin.castLE (Nat.sub_le N r) j))) =
      N - r :=
    (finrank_span_eq_card (hli.comp _ (Fin.castLE_injective (Nat.sub_le N r)))).trans
      (Fintype.card_fin _)
  have hrk : N - r ≤ finrank K (LinearMap.range (diffₗ N K ^ r)) := by
    have := Submodule.finrank_mono hrange
    rw [h2] at this
    exact this
  have hrn := LinearMap.finrank_range_add_finrank_ker (diffₗ N K ^ r)
  rw [Module.finrank_fintype_fun_eq_card, ZMod.card] at hrn
  omega

/-- Same, with the spanning set written as `{chooseFun i : i < r}`. -/
theorem ker_diffₗ_pow' (b : ℕ) [NeZero (p ^ b)] (r : ℕ) (hr : r ≤ p ^ b) :
    LinearMap.ker (diffₗ (p ^ b) K ^ r) =
      Submodule.span K (chooseFun (p ^ b) K '' {i | (i : ℕ) < r}) := by
  rw [ker_diffₗ_pow b r hr]
  congr 1
  ext f
  constructor
  · rintro ⟨i, rfl⟩
    exact ⟨Fin.castLE hr i, by simp, rfl⟩
  · rintro ⟨i, hi, rfl⟩
    exact ⟨⟨i, hi⟩, rfl⟩

end Kernel

end ConjeturaC

#print axioms ConjeturaC.mahlerZ_coeff_eq_zero
#print axioms ConjeturaC.mahlerZ_coeff_eq_zero_field
#print axioms ConjeturaC.ker_diffₗ_pow
#print axioms ConjeturaC.ker_diffₗ_pow'
