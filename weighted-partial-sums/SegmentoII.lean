/- SegmentoII.lean — digit recursion for weighted sums of a completely multiplicative,
   q-periodic-off-multiples function, zeros from base-q digits, and the complete-sequence
   (subset-sum) lemma with its pairs version.
   Author: Carles Marín  <karlesmarin@gmail.com>   (with Claude, Anthropic, as assistant)

   Formalizes three results of paper II (orbit-pair, segmento II), over ℤ-valued χ : ℕ → ℤ with
     (mult) χ (a*b) = χ a * χ b,  (per) χ (j*q + r) = χ r for 1 ≤ r < q,  (s) χ q = s,
     (sum0) ∑_{n=1}^{q-1} χ n = 0.
   * `digit_recursion`   : S (q*k + r) = s*q*S k + k*(Tq + q*C r) + T r      (r < q)
   * `zero_step`, `zeros_from_digits` : zeros of S propagate along base-q digits in R.
   * `complete_list`, `complete_range`, `pairs_subset_sum` : complete-sequence lemma.

   Fast-loop build (specific imports):
     lake env lean <abs path>/SegmentoII.lean     (from E:\proyectos\godsil-gutman-lean) -/
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Lemmas
import Mathlib.Algebra.BigOperators.Group.Finset.Piecewise
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Group.Nat.Even
import Mathlib.Data.Nat.Digits.Defs
import Mathlib.Data.List.Induction
import Mathlib.Tactic.Ring
import Mathlib.Tactic.LinearCombination

open Finset

namespace SegmentoII

/-- `S m = ∑_{n=1}^{m} n χ(n)`. -/
def S (χ : ℕ → ℤ) (m : ℕ) : ℤ := ∑ n ∈ Icc 1 m, (n : ℤ) * χ n

/-- `C r = ∑_{n=1}^{r} χ(n)`. -/
def C (χ : ℕ → ℤ) (r : ℕ) : ℤ := ∑ n ∈ Icc 1 r, χ n

/-- `T r = ∑_{n=1}^{r} n χ(n)` (same as `S`; used for `r < q`). -/
def T (χ : ℕ → ℤ) (r : ℕ) : ℤ := ∑ n ∈ Icc 1 r, (n : ℤ) * χ n

/-- `Tq = T (q-1)`. -/
def Tq (χ : ℕ → ℤ) (q : ℕ) : ℤ := T χ (q - 1)

lemma S_zero (χ : ℕ → ℤ) : S χ 0 = 0 := by simp [S]
lemma T_zero (χ : ℕ → ℤ) : T χ 0 = 0 := by simp [T]
lemma C_zero (χ : ℕ → ℤ) : C χ 0 = 0 := by simp [C]

lemma S_succ (χ : ℕ → ℤ) (n : ℕ) : S χ (n + 1) = S χ n + ((n : ℤ) + 1) * χ (n + 1) := by
  unfold S; rw [Finset.sum_Icc_succ_top (by omega)]; push_cast; ring

lemma T_succ (χ : ℕ → ℤ) (n : ℕ) : T χ (n + 1) = T χ n + ((n : ℤ) + 1) * χ (n + 1) :=
  S_succ χ n

lemma C_succ (χ : ℕ → ℤ) (n : ℕ) : C χ (n + 1) = C χ n + χ (n + 1) := by
  unfold C; rw [Finset.sum_Icc_succ_top (by omega)]

/-- Inner step: `r ↦ r+1` inside a block (uses periodicity). -/
private lemma step_r (q : ℕ) (χ : ℕ → ℤ) (s : ℤ)
    (hper : ∀ j r, 1 ≤ r → r < q → χ (j * q + r) = χ r) (k r : ℕ) (hr : r + 1 < q)
    (ih : S χ (q * k + r) = s * q * S χ k + k * (Tq χ q + q * C χ r) + T χ r) :
    S χ (q * k + (r + 1)) = s * q * S χ k + k * (Tq χ q + q * C χ (r + 1)) + T χ (r + 1) := by
  have hχ : χ (q * k + r + 1) = χ (r + 1) := by
    rw [mul_comm q k, add_assoc]; exact hper k (r + 1) (by omega) hr
  have e : q * k + (r + 1) = (q * k + r) + 1 := by omega
  rw [e, S_succ, C_succ, T_succ, hχ, ih]
  push_cast; ring

/-- Block step: `(k, q-1) ↦ (k+1, 0)` (uses multiplicativity, `χ q = s`, `sum0`). -/
private lemma step_k (q : ℕ) (hq : 1 ≤ q) (χ : ℕ → ℤ) (s : ℤ)
    (hmul : ∀ a b, χ (a * b) = χ a * χ b) (hs : χ q = s) (hsum0 : C χ (q - 1) = 0) (k : ℕ)
    (ih : S χ (q * k + (q - 1)) = s * q * S χ k + k * (Tq χ q + q * C χ (q - 1)) + T χ (q - 1)) :
    S χ (q * (k + 1) + 0) = s * q * S χ (k + 1) + ((k + 1 : ℕ) : ℤ) * (Tq χ q + q * C χ 0)
      + T χ 0 := by
  have e : q * (k + 1) + 0 = (q * k + (q - 1)) + 1 := by
    rw [Nat.mul_succ]; omega
  have hχ : χ (q * k + (q - 1) + 1) = s * χ (k + 1) := by
    rw [← e, add_zero, hmul, hs]
  have hcast : ((q * k + (q - 1) : ℕ) : ℤ) + 1 = (q : ℤ) * ((k : ℤ) + 1) := by
    have h2 : q * k + (q - 1) + 1 = q * (k + 1) := by rw [← e, add_zero]
    exact_mod_cast h2
  have hTq : Tq χ q = T χ (q - 1) := rfl
  rw [e, S_succ, hχ, ih, hcast, S_succ, C_zero, T_zero, hsum0, hTq]
  push_cast; ring

/-- **Theorem 1 (digit recursion).** For `r < q`:
`S (q*k + r) = s*q*S k + k*(Tq + q*C r) + T r`. -/
theorem digit_recursion (q : ℕ) (hq : 2 ≤ q) (χ : ℕ → ℤ) (s : ℤ)
    (hmul : ∀ a b, χ (a * b) = χ a * χ b)
    (hper : ∀ j r, 1 ≤ r → r < q → χ (j * q + r) = χ r)
    (hs : χ q = s) (hsum0 : C χ (q - 1) = 0) (k r : ℕ) (hr : r < q) :
    S χ (q * k + r) = s * q * S χ k + k * (Tq χ q + q * C χ r) + T χ r := by
  -- from `P k 0` get `P k r` for every `r < q`
  have fill : ∀ k, (S χ (q * k + 0) = s * q * S χ k + k * (Tq χ q + q * C χ 0) + T χ 0) →
      ∀ r, r < q → S χ (q * k + r) = s * q * S χ k + k * (Tq χ q + q * C χ r) + T χ r := by
    intro k h0 r
    induction r with
    | zero => intro _; exact h0
    | succ r ihr => intro hr; exact step_r q χ s hper k r hr (ihr (by omega))
  revert r
  induction k with
  | zero =>
    apply fill 0
    simp [S_zero, T_zero, C_zero]
  | succ k ih =>
    apply fill (k + 1)
    exact step_k q (by omega) χ s hmul hs hsum0 k (ih (q - 1) (by omega))

/-- The digit set `R = {r | r < q ∧ Tq + q*C r = 0 ∧ T r = 0}`. -/
def InR (χ : ℕ → ℤ) (q r : ℕ) : Prop := r < q ∧ Tq χ q + (q : ℤ) * C χ r = 0 ∧ T χ r = 0

/-- **Theorem 2 (zeros from digits), one step.** -/
theorem zero_step (q : ℕ) (hq : 2 ≤ q) (χ : ℕ → ℤ) (s : ℤ)
    (hmul : ∀ a b, χ (a * b) = χ a * χ b)
    (hper : ∀ j r, 1 ≤ r → r < q → χ (j * q + r) = χ r)
    (hs : χ q = s) (hsum0 : C χ (q - 1) = 0) (k r : ℕ)
    (hk : S χ k = 0) (hR : InR χ q r) : S χ (q * k + r) = 0 := by
  obtain ⟨hr, h1, h2⟩ := hR
  rw [digit_recursion q hq χ s hmul hper hs hsum0 k r hr, hk, h1, h2]
  ring

/-- **Theorem 2 (zeros from digits), corollary.** If every base-`q` digit of `m` lies in `R`,
then `S m = 0`. -/
theorem zeros_from_digits (q : ℕ) (hq : 2 ≤ q) (χ : ℕ → ℤ) (s : ℤ)
    (hmul : ∀ a b, χ (a * b) = χ a * χ b)
    (hper : ∀ j r, 1 ≤ r → r < q → χ (j * q + r) = χ r)
    (hs : χ q = s) (hsum0 : C χ (q - 1) = 0) (m : ℕ)
    (hm : ∀ d ∈ Nat.digits q m, InR χ q d) : S χ m = 0 := by
  induction m using Nat.strong_induction_on with
  | _ m ih =>
    rcases Nat.eq_zero_or_pos m with rfl | hpos
    · exact S_zero χ
    · rw [Nat.digits_def' (by omega) hpos] at hm
      have hlt : m / q < m := Nat.div_lt_self hpos (by omega)
      have h0 := ih (m / q) hlt (fun d hd => hm d (List.mem_cons_of_mem _ hd))
      have hR := hm (m % q) (List.mem_cons_self ..)
      have := zero_step q hq χ s hmul hper hs hsum0 (m / q) (m % q) h0 hR
      rwa [Nat.div_add_mod] at this

/-- **Theorem 3 (complete sequence), list form.** If each entry is at most `1 +` the sum of the
earlier entries, every `x ≤ d.sum` is the sum of a sublist of `d`.
(The hypothesis "every entry ≥ 1" is not needed and has been dropped.) -/
theorem complete_list (d : List ℕ)
    (hd : ∀ l₁ a l₂, d = l₁ ++ a :: l₂ → a ≤ 1 + l₁.sum) :
    ∀ x ≤ d.sum, ∃ l : List ℕ, l.Sublist d ∧ l.sum = x := by
  induction d using List.reverseRecOn with
  | nil =>
    intro x hx
    exact ⟨[], List.Sublist.slnil, by simp at hx ⊢; omega⟩
  | append_singleton l a ih =>
    have hl : ∀ l₁ b l₂, l = l₁ ++ b :: l₂ → b ≤ 1 + l₁.sum :=
      fun l₁ b l₂ h => hd l₁ b (l₂ ++ [a]) (by rw [h]; simp)
    have ha : a ≤ 1 + l.sum := hd l a [] rfl
    intro x hx
    simp only [List.sum_append, List.sum_cons, List.sum_nil, add_zero] at hx
    by_cases hx' : x ≤ l.sum
    · obtain ⟨t, ht, hsum⟩ := ih hl x hx'
      exact ⟨t, ht.trans (List.sublist_append_left l [a]), hsum⟩
    · obtain ⟨t, ht, hsum⟩ := ih hl (x - a) (by omega)
      refine ⟨t ++ [a], ht.append (List.Sublist.refl _), ?_⟩
      simp only [List.sum_append, List.sum_cons, List.sum_nil, add_zero, hsum]
      omega

/-- **Theorem 3 (complete sequence), indexed form** (`d j`, `j < t`, subsets of `range t`). -/
theorem complete_range (d : ℕ → ℕ) (t : ℕ)
    (hd : ∀ j < t, d j ≤ 1 + ∑ i ∈ range j, d i) :
    ∀ y ≤ ∑ i ∈ range t, d i, ∃ A ⊆ range t, ∑ i ∈ A, d i = y := by
  induction t with
  | zero =>
    intro y hy
    exact ⟨∅, empty_subset _, by simp at hy ⊢; omega⟩
  | succ t ih =>
    have hd' : ∀ j < t, d j ≤ 1 + ∑ i ∈ range j, d i := fun j hj => hd j (by omega)
    have ht := hd t (by omega)
    intro y hy
    rw [sum_range_succ] at hy
    by_cases hy' : y ≤ ∑ i ∈ range t, d i
    · obtain ⟨A, hA, hsum⟩ := ih hd' y hy'
      exact ⟨A, hA.trans (range_subset_range.mpr (by omega)), hsum⟩
    · obtain ⟨A, hA, hsum⟩ := ih hd' (y - d t) (by omega)
      have hnot : t ∉ A := fun h => by simpa using hA h
      refine ⟨insert t A, ?_, ?_⟩
      · intro i hi
        rcases mem_insert.mp hi with rfl | hi
        · simp
        · exact range_subset_range.mpr (by omega) (hA hi)
      · rw [sum_insert hnot, hsum]; omega

/-- Choosing `v j` on `A` and `u j` off `A` adds `2 * ∑_A d` to `U = ∑ u`. -/
lemma choice_sum (t : ℕ) (d u v : ℕ → ℕ) (hv : ∀ j, v j = u j + 2 * d j)
    (A : Finset ℕ) (hA : A ⊆ range t) :
    ∑ j ∈ range t, (if j ∈ A then v j else u j)
      = ∑ j ∈ range t, u j + 2 * ∑ j ∈ A, d j := by
  have : ∀ j, (if j ∈ A then v j else u j) = u j + (if j ∈ A then 2 * d j else 0) := by
    intro j; split_ifs <;> simp [hv]
  simp only [this, sum_add_distrib, sum_ite_mem, inter_eq_right.mpr hA, mul_sum]

/-- **Theorem 3, pairs version.** Pairs `(u j, v j)` with `v j = u j + 2 d j`, `d` complete,
`q0` odd, `U = ∑ u`, `H = ∑ d`. Every `x ∈ [U + q0 - 1, U + 2H + 1]` is a sum of a selection:
optionally `q0`, plus one element from each pair (`v j` for `j ∈ A`, `u j` otherwise).
(For distinct `q0, u j, v j` this is a subset sum of `{q0, u j, v j}`; the hypothesis
`q0 ≤ 2H + 1` is not needed and has been dropped.) -/
theorem pairs_subset_sum (t : ℕ) (d u v : ℕ → ℕ) (q0 : ℕ) (hq0 : Odd q0)
    (hv : ∀ j, v j = u j + 2 * d j)
    (hd : ∀ j < t, d j ≤ 1 + ∑ i ∈ range j, d i) (x : ℕ)
    (hx1 : ∑ j ∈ range t, u j + q0 - 1 ≤ x)
    (hx2 : x ≤ ∑ j ∈ range t, u j + 2 * ∑ j ∈ range t, d j + 1) :
    ∃ b : Bool, ∃ A ⊆ range t,
      x = (if b then q0 else 0) + ∑ j ∈ range t, (if j ∈ A then v j else u j) := by
  set U := ∑ j ∈ range t, u j with hU
  set H := ∑ j ∈ range t, d j with hH
  obtain ⟨c, hc⟩ := hq0
  have hUx : U ≤ x := by omega
  rcases Nat.even_or_odd (x - U) with ⟨y, hy⟩ | ⟨y, hy⟩
  · -- x = U + 2y, y ≤ H
    obtain ⟨A, hA, hsum⟩ := complete_range d t hd y (by omega)
    refine ⟨false, A, hA, ?_⟩
    rw [choice_sum t d u v hv A hA, hsum]; simp; omega
  · -- x = U + q0 + 2y', y' ≤ H
    obtain ⟨A, hA, hsum⟩ := complete_range d t hd (y - c) (by omega)
    refine ⟨true, A, hA, ?_⟩
    rw [choice_sum t d u v hv A hA, hsum]; simp; omega

end SegmentoII

#print axioms SegmentoII.digit_recursion
#print axioms SegmentoII.zero_step
#print axioms SegmentoII.zeros_from_digits
#print axioms SegmentoII.complete_list
#print axioms SegmentoII.complete_range
#print axioms SegmentoII.pairs_subset_sum
