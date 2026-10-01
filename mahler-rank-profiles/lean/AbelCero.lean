/- AbelCero.lean — Abel's lemma at j = 0, purely algebraic (finite Fourier, no L-functions).
   Author: Carles Marín  <karlesmarin@gmail.com>   (with Claude, Anthropic, as assistant)

   Setting: `K` a field of characteristic 0, `N = 2m+1` odd, `N ≥ 3`, `ζ` a primitive `N`-th root
   of unity, `U(Z) = 2Z/(Z^2-1)`, `B1 χ = (1/N) Σ_a χ(a) a.val`.

   (a) `sub_one_mul_inv_eq`   : `w^N = 1`, `w ≠ 1` ⇒ `(w-1)⁻¹ = N⁻¹ Σ_{k<N} k w^k`.
   (b) `U_pow_eq`             : `c ≠ 0` ⇒ `U(ζ^c) = (2/N) Σ_{k<N} k ζ^{c(2k+1)}`.
   (c) `sum_mul_char_two_mul_add_one` : `Σ_{k<N} k χ(2k+1) = (1 - χ 2) Σ_a χ(a) a.val` (χ ≠ 1).
   (d) `abel_zero_of_gauss`   : with an explicit Gauss-sum hypothesis
         `hgauss : ∀ a, Σ_c ψ c ζ^{(a c).val} = χ a τ`,
         `Σ_c ψ c U(ζ^c) = 2 τ (1 - χ 2) B1 χ`.
   (e) `abel_zero`            : `χ` a primitive Dirichlet character; `ψ = χ⁻¹`,
         `τ = gaussSum χ⁻¹ (zmodChar N ζ)`; `hgauss` discharged by Mathlib
         `gaussSum_mulShift_of_isPrimitive` + `conductor_inv`.

   Fast-loop build: lake env lean <this file>  (from E:\proyectos\godsil-gutman-lean). -/
import Mathlib.NumberTheory.DirichletCharacter.GaussSum
import Mathlib.NumberTheory.LegendreSymbol.AddCharacter
import Mathlib.RingTheory.RootsOfUnity.PrimitiveRoots
import Mathlib.Algebra.Field.GeomSum
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination

open Finset

namespace ConjeturaC

/-- `U(Z) = 2Z/(Z^2 - 1)`. -/
def U {K : Type*} [Field K] (z : K) : K := 2 * z / (z ^ 2 - 1)

/-- `B_{1,χ} = (1/N) Σ_a χ(a) a.val`. -/
noncomputable def B1 {K : Type*} [Field K] {N : ℕ} [NeZero N] (χ : MulChar (ZMod N) K) : K :=
  (N : K)⁻¹ * ∑ a : ZMod N, χ a * (a.val : K)

/-! ### (a) finite Fourier identity for `1/(w-1)` -/

theorem sub_one_mul_sum_mul_pow {K : Type*} [CommRing K] (w : K) (n : ℕ) :
    (w - 1) * ∑ k ∈ range n, (k : K) * w ^ k
      = (n : K) * w ^ n - w ^ n - ∑ k ∈ range n, w ^ k + 1 := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [sum_range_succ, mul_add, ih, sum_range_succ]
    push_cast; ring

theorem sub_one_mul_inv_eq {K : Type*} [Field K] [CharZero K] {N : ℕ} (hN : 0 < N) {w : K}
    (hw : w ^ N = 1) (hw1 : w ≠ 1) :
    (w - 1)⁻¹ = (N : K)⁻¹ * ∑ k ∈ range N, (k : K) * w ^ k := by
  have hgeom : ∑ k ∈ range N, w ^ k = 0 := by
    rw [geom_sum_eq hw1, hw, sub_self, zero_div]
  have key := sub_one_mul_sum_mul_pow w N
  rw [hw, hgeom] at key
  have hN' : (N : K) ≠ 0 := by exact_mod_cast hN.ne'
  have hw1' : w - 1 ≠ 0 := sub_ne_zero.mpr hw1
  rw [eq_comm, inv_mul_eq_iff_eq_mul₀ hN', eq_mul_inv_iff_mul_eq₀ hw1']
  linear_combination key

/-! ### sums over `ZMod N` as sums over `range N` -/

theorem sum_zmod_eq_sum_range {M : Type*} [AddCommMonoid M] (N : ℕ) [NeZero N]
    (f : ZMod N → M) : ∑ a : ZMod N, f a = ∑ n ∈ range N, f n := by
  refine Finset.sum_nbij' (fun a => a.val) (fun n => (n : ZMod N)) ?_ ?_ ?_ ?_ ?_
  · intro a _; simp [ZMod.val_lt]
  · intro n _; simp
  · intro a _; simp
  · intro n hn
    simp [ZMod.val_natCast, Nat.mod_eq_of_lt (mem_range.1 hn)]
  · intro a _; simp

theorem sum_range_two_mul_add_one_parity {M : Type*} [AddCommMonoid M] (g : ℕ → M) (m : ℕ) :
    ∑ n ∈ range (2 * m + 1), g n = g 0 + ∑ i ∈ range m, (g (2 * i + 1) + g (2 * i + 2)) := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [show 2 * (m + 1) + 1 = 2 * m + 1 + 1 + 1 by ring, sum_range_succ, sum_range_succ, ih,
      sum_range_succ, add_assoc, add_assoc]

theorem sum_range_two_mul_add_one_halves {M : Type*} [AddCommMonoid M] (g : ℕ → M) (m : ℕ) :
    ∑ n ∈ range (2 * m + 1), g n
      = ∑ n ∈ range m, g n + g m + ∑ j ∈ range m, g (m + 1 + j) := by
  rw [show 2 * m + 1 = (m + 1) + m by ring, sum_range_add, sum_range_succ]

/-! ### (c) the reindexing `a = 2k+1` -/

theorem cast_add_self {m r : ℕ} :
    (((2 * m + 1) + r : ℕ) : ZMod (2 * m + 1)) = (r : ZMod (2 * m + 1)) := by
  rw [Nat.cast_add, ZMod.natCast_self, zero_add]

theorem sum_mul_char_two_mul_add_one {K : Type*} [Field K] {m : ℕ} (hm : 1 ≤ m)
    (χ : MulChar (ZMod (2 * m + 1)) K) (hχ : χ ≠ 1) :
    ∑ k ∈ range (2 * m + 1), (k : K) * χ ((2 * k + 1 : ℕ) : ZMod (2 * m + 1))
      = (1 - χ 2) * ∑ a : ZMod (2 * m + 1), χ a * (a.val : K) := by
  haveI : Fact (1 < 2 * m + 1) := ⟨by omega⟩
  have hχ0 : χ 0 = 0 := MulChar.map_nonunit χ not_isUnit_zero
  -- the sum of `χ` vanishes
  have hZ : ∑ i ∈ range m, (χ ((2 * i + 1 : ℕ) : ZMod (2 * m + 1))
      + χ ((2 * i + 2 : ℕ) : ZMod (2 * m + 1))) = 0 := by
    have h := MulChar.sum_eq_zero_of_ne_one hχ
    rw [sum_zmod_eq_sum_range, sum_range_two_mul_add_one_parity] at h
    simpa [hχ0] using h
  -- `T` split into low and high halves
  have hT : ∑ k ∈ range (2 * m + 1), (k : K) * χ ((2 * k + 1 : ℕ) : ZMod (2 * m + 1))
      = ∑ i ∈ range m, ((i : K) * χ ((2 * i + 1 : ℕ) : ZMod (2 * m + 1))
          + ((m + 1 + i : ℕ) : K) * χ ((2 * i + 2 : ℕ) : ZMod (2 * m + 1))) := by
    rw [sum_range_two_mul_add_one_halves, ZMod.natCast_self, hχ0, mul_zero, add_zero,
      ← sum_add_distrib]
    refine sum_congr rfl fun j _ => ?_
    rw [show 2 * (m + 1 + j) + 1 = (2 * m + 1) + (2 * j + 2) by ring, cast_add_self]
  -- `S1`
  have hS : ∑ a : ZMod (2 * m + 1), χ a * (a.val : K)
      = ∑ i ∈ range m, (χ ((2 * i + 1 : ℕ) : ZMod (2 * m + 1)) * ((2 * i + 1 : ℕ) : K)
          + χ ((2 * i + 2 : ℕ) : ZMod (2 * m + 1)) * ((2 * i + 2 : ℕ) : K)) := by
    rw [sum_zmod_eq_sum_range]
    rw [sum_congr rfl fun n hn => by
      rw [ZMod.val_natCast, Nat.mod_eq_of_lt (mem_range.1 hn)]]
    rw [sum_range_two_mul_add_one_parity]
    simp
  -- `χ(2) S1`
  have hX : χ 2 * ∑ a : ZMod (2 * m + 1), χ a * (a.val : K)
      = ∑ i ∈ range m, (χ ((2 * i + 2 : ℕ) : ZMod (2 * m + 1)) * ((i + 1 : ℕ) : K)
          + χ ((2 * i + 1 : ℕ) : ZMod (2 * m + 1)) * ((m + 1 + i : ℕ) : K)) := by
    rw [mul_sum]
    rw [sum_congr rfl fun a _ => (by rw [← mul_assoc, ← map_mul] :
      χ 2 * (χ a * (a.val : K)) = χ (2 * a) * (a.val : K))]
    rw [sum_zmod_eq_sum_range]
    rw [sum_congr rfl fun (n : ℕ) hn => (by
      rw [ZMod.val_natCast, Nat.mod_eq_of_lt (mem_range.1 hn)]; push_cast; ring :
      χ (2 * (n : ZMod (2 * m + 1))) * (((n : ZMod (2 * m + 1)).val : ℕ) : K)
        = χ ((2 * n : ℕ) : ZMod (2 * m + 1)) * (n : K))]
    rw [sum_range_two_mul_add_one_halves, ← sum_range_succ, sum_range_succ']
    simp only [mul_zero, Nat.cast_zero, add_zero]
    rw [← sum_add_distrib]
    refine sum_congr rfl fun j _ => ?_
    rw [show 2 * (m + 1 + j) = (2 * m + 1) + (2 * j + 1) by ring, cast_add_self,
      show 2 * (j + 1) = 2 * j + 2 by ring]
  rw [hT, sub_mul, one_mul, hX, hS, ← sub_eq_zero, ← sum_sub_distrib, ← sum_sub_distrib]
  calc _ = ∑ i ∈ range m, (m : K) * (χ ((2 * i + 1 : ℕ) : ZMod (2 * m + 1))
            + χ ((2 * i + 2 : ℕ) : ZMod (2 * m + 1))) :=
        sum_congr rfl fun i _ => by push_cast; ring
    _ = 0 := by rw [← mul_sum, hZ, mul_zero]

/-! ### (b) `U(ζ^c)` as a finite Fourier sum -/

theorem U_pow_eq {K : Type*} [Field K] [CharZero K] {m : ℕ} {ζ : K}
    (hζ : IsPrimitiveRoot ζ (2 * m + 1)) (c : ZMod (2 * m + 1)) (hc : c ≠ 0) :
    U (ζ ^ c.val) = (2 / ((2 * m + 1 : ℕ) : K))
      * ∑ k ∈ range (2 * m + 1), (k : K) * ζ ^ (c.val * (2 * k + 1)) := by
  have hw : (ζ ^ (2 * c.val)) ^ (2 * m + 1) = 1 := by
    rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]
  have hw1 : ζ ^ (2 * c.val) ≠ 1 := by
    intro h
    rw [hζ.pow_eq_one_iff_dvd] at h
    have hcop : Nat.Coprime (2 * m + 1) 2 := Nat.coprime_two_right.2 (odd_two_mul_add_one m)
    have h2 : (2 * m + 1) ∣ c.val := (Nat.Coprime.dvd_mul_left hcop).1 h
    have h3 : c.val = 0 := Nat.eq_zero_of_dvd_of_lt h2 (ZMod.val_lt c)
    exact hc ((ZMod.val_eq_zero c).1 h3)
  have hsq : (ζ ^ c.val) ^ 2 = ζ ^ (2 * c.val) := by rw [← pow_mul, mul_comm]
  rw [U, div_eq_mul_inv, hsq, sub_one_mul_inv_eq (by omega) hw hw1, mul_sum, mul_sum, mul_sum]
  refine sum_congr rfl fun k _ => ?_
  rw [show ζ ^ (c.val * (2 * k + 1)) = ζ ^ c.val * (ζ ^ (2 * c.val)) ^ k by
    rw [← pow_mul, ← pow_add, show c.val * (2 * k + 1) = c.val + 2 * c.val * k by ring]]
  ring

/-! ### (d) Abel's lemma at `j = 0`, Gauss-sum property as hypothesis -/

theorem abel_zero_of_gauss {K : Type*} [Field K] [CharZero K] {m : ℕ} (hm : 1 ≤ m) {ζ : K}
    (hζ : IsPrimitiveRoot ζ (2 * m + 1)) (χ ψ : MulChar (ZMod (2 * m + 1)) K) (hχ : χ ≠ 1)
    (τ : K) (hgauss : ∀ a : ZMod (2 * m + 1), ∑ c, ψ c * ζ ^ (a * c).val = χ a * τ) :
    ∑ c : ZMod (2 * m + 1), ψ c * U (ζ ^ c.val) = 2 * τ * (1 - χ 2) * B1 χ := by
  haveI : Fact (1 < 2 * m + 1) := ⟨by omega⟩
  have hψ0 : ψ 0 = 0 := MulChar.map_nonunit ψ not_isUnit_zero
  have hζN : ζ ^ (2 * m + 1) = 1 := hζ.pow_eq_one
  -- step 1: expand each `U(ζ^c)`
  have h1 : ∑ c : ZMod (2 * m + 1), ψ c * U (ζ ^ c.val)
      = ∑ c : ZMod (2 * m + 1), ψ c * ((2 / ((2 * m + 1 : ℕ) : K))
          * ∑ k ∈ range (2 * m + 1), (k : K) * ζ ^ (c.val * (2 * k + 1))) := by
    refine sum_congr rfl fun c _ => ?_
    by_cases hc : c = 0
    · rw [hc, hψ0, zero_mul, zero_mul]
    · rw [U_pow_eq hζ c hc]
  -- step 2: swap sums and apply the Gauss property
  have h2 : ∀ k : ℕ, ∑ c : ZMod (2 * m + 1), ψ c * ζ ^ (c.val * (2 * k + 1))
      = χ ((2 * k + 1 : ℕ) : ZMod (2 * m + 1)) * τ := by
    intro k
    rw [← hgauss]
    refine sum_congr rfl fun c _ => ?_
    congr 1
    rw [pow_eq_pow_mod _ hζN, pow_eq_pow_mod ((((2 * k + 1 : ℕ) : ZMod (2 * m + 1)) * c).val) hζN]
    congr 1
    rw [← ZMod.natCast_eq_natCast_iff', ZMod.natCast_zmod_val]
    push_cast
    simp only [ZMod.natCast_zmod_val]
    ring
  rw [h1]
  have h3 : ∑ c : ZMod (2 * m + 1), ψ c * ((2 / ((2 * m + 1 : ℕ) : K))
          * ∑ k ∈ range (2 * m + 1), (k : K) * ζ ^ (c.val * (2 * k + 1)))
      = (2 / ((2 * m + 1 : ℕ) : K)) * ∑ k ∈ range (2 * m + 1),
          (k : K) * ∑ c : ZMod (2 * m + 1), ψ c * ζ ^ (c.val * (2 * k + 1)) := by
    simp only [mul_sum]
    rw [sum_comm]
    refine sum_congr rfl fun k _ => sum_congr rfl fun c _ => by ring
  rw [h3]
  simp only [h2]
  rw [show ∑ k ∈ range (2 * m + 1), (k : K) * (χ ((2 * k + 1 : ℕ) : ZMod (2 * m + 1)) * τ)
      = τ * ∑ k ∈ range (2 * m + 1), (k : K) * χ ((2 * k + 1 : ℕ) : ZMod (2 * m + 1)) by
    rw [mul_sum]; exact sum_congr rfl fun k _ => by ring]
  rw [sum_mul_char_two_mul_add_one hm χ hχ, B1]
  field_simp

/-! ### (e) primitive Dirichlet characters: the Gauss property from Mathlib -/

theorem abel_zero {K : Type*} [Field K] [CharZero K] {m : ℕ} (hm : 1 ≤ m) {ζ : K}
    (hζ : IsPrimitiveRoot ζ (2 * m + 1)) (χ : DirichletCharacter K (2 * m + 1))
    (hχ : χ.IsPrimitive) :
    ∑ c : ZMod (2 * m + 1), χ⁻¹ c * U (ζ ^ c.val)
      = 2 * gaussSum χ⁻¹ (AddChar.zmodChar (2 * m + 1) hζ.pow_eq_one) * (1 - χ 2) * B1 χ := by
  have hχ1 : χ ≠ 1 := by
    rintro rfl
    have := hχ
    rw [DirichletCharacter.isPrimitive_def, DirichletCharacter.conductor_one] at this
    omega
  have hχinv : (χ⁻¹).IsPrimitive := by
    rw [DirichletCharacter.isPrimitive_def, DirichletCharacter.conductor_inv]; exact hχ
  refine abel_zero_of_gauss hm hζ χ χ⁻¹ hχ1 _ fun a => ?_
  have h := gaussSum_mulShift_of_isPrimitive (AddChar.zmodChar (2 * m + 1) hζ.pow_eq_one)
    hχinv a
  rw [inv_inv] at h
  rw [← h, gaussSum]
  rfl

end ConjeturaC

#print axioms ConjeturaC.sub_one_mul_inv_eq
#print axioms ConjeturaC.sum_mul_char_two_mul_add_one
#print axioms ConjeturaC.U_pow_eq
#print axioms ConjeturaC.abel_zero_of_gauss
#print axioms ConjeturaC.abel_zero
