# -*- coding: utf-8 -*-
"""The proof of the Question on the short form of the sign, checked step by step.

sign_ayyer_idiom.py verifies the STATEMENT. This script verifies the three steps of the PROOF, so
that a disagreement can be localised:

  STEP 1  the parity lemma.  For a word w and a position p carrying the letter c,
          let I(p) be the number of inversions of w involving p.  Then
                I(p) = (p-1) + #{letters of w strictly below c} + #{occurrences of c before p}   (mod 2).
          Checked on random words, where it can fail.

  STEP 2  the consequence for our word.  With the blocks ordered so that r_A <= r_B,
                inv(w) - inv(b_S) = j_A1 + j_B1 + r_A + r_B + 1 + [a1 < b1]   (mod 2),
          one formula covering the two-class and the size-three profiles.

  STEP 3  the arithmetic fact that closes it,  floor(t/2) + t + C(t+3,2) = 1  (mod 2).

  CONTROL the same STEP 2 with the ordering hypothesis dropped.  It must fail, because
          #{letters below r_A} picks up the extra occurrence of r_B when r_B < r_A.

Authors: Carles Marin, Claude (AI assistant)."""
import random

from law_control import partitions
from theorem_full import setup

TMAX, NMAX = 9, 15          # t = 2..8, |lambda| <= 14


def inv_count(word):
    return sum(1 for i in range(len(word)) for j in range(i + 1, len(word)) if word[i] > word[j])


def inversions_at(word, p):
    """number of inversions of `word` in which position p takes part."""
    c = word[p]
    return (sum(1 for j in range(p) if word[j] > c)
            + sum(1 for j in range(p + 1, len(word)) if word[j] < c))


# ----------------------------------------------------------------- STEP 1
def step1(trials=4000, seed=20260806):
    rng = random.Random(seed)
    ok = bad = 0
    for _ in range(trials):
        n = rng.randint(2, 12)
        alpha = rng.randint(2, 7)
        word = [rng.randrange(alpha) for _ in range(n)]
        p = rng.randrange(n)
        c = word[p]
        lhs = inversions_at(word, p) % 2
        rhs = (p + sum(1 for x in word if x < c) + sum(1 for j in range(p) if word[j] == c)) % 2
        if lhs == rhs:
            ok += 1
        else:
            bad += 1
    return ok, bad


# ----------------------------------------------------------------- STEP 0
def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def step0():
    """sgn(sigma) = (-1)^inv(w): sigma is the stable sort of the residue word, because beta
    decreases with the column index, so 'values decreasing within a class' is 'columns increasing'."""
    ok = bad = 0
    for t in range(2, TMAX):
        N = t + 2
        for n in range(0, NMAX):
            for lam in partitions(n, N):
                st = setup(lam, t)
                if st is None:
                    continue
                beta = st[0]
                order = []
                for a in range(t):
                    order += [j for j, b in enumerate(beta) if b % t == a]
                w = [b % t for b in beta]
                if perm_sign(order) == (-1) ** inv_count(w):
                    ok += 1
                else:
                    bad += 1
    return ok, bad


# ----------------------------------------------------------------- STEP 2
def step2(order_by_residue):
    ok = bad = 0
    for t in range(2, TMAX):
        N = t + 2
        for n in range(0, NMAX):
            for lam in partitions(n, N):
                st = setup(lam, t)
                if st is None:
                    continue
                beta, Ac, Bc = st
                rA, rB = beta[Ac[0]] % t, beta[Bc[0]] % t
                if order_by_residue and rA > rB:
                    Ac, Bc = Bc, Ac
                    rA, rB = rB, rA
                jA, jB = Ac[0], Bc[0]
                if jA == jB:
                    continue
                a1, a2 = beta[Ac[0]], beta[Ac[1]]
                b1, b2 = beta[Bc[0]], beta[Bc[1]]
                if a1 + a2 - b1 - b2 == 0:
                    continue
                w = [b % t for b in beta]
                bS = [w[j] for j in range(N) if j not in (jA, jB)]
                D = (inv_count(w) - inv_count(bS)) % 2
                pred = (jA + jB + rA + rB + 1 + (1 if a1 < b1 else 0)) % 2
                if D == pred:
                    ok += 1
                else:
                    bad += 1
    return ok, bad


# ----------------------------------------------------------------- STEP 3
def step3(tmax=200):
    bad = [t for t in range(2, tmax)
           if (t // 2 + t + (t + 3) * (t + 2) // 2) % 2 != 1]
    return tmax - 2 - len(bad), bad


if __name__ == "__main__":
    ok, bad = step0()
    print("STEP 0  sgn(sigma) = (-1)^inv(w):          %4d ok, %d fail" % (ok, bad))

    ok, bad = step1()
    print("STEP 1  parity lemma on random words:      %4d ok, %d fail" % (ok, bad))

    ok, bad = step2(order_by_residue=True)
    print("STEP 2  inv(w) - inv(b_S), blocks ordered: %4d ok, %d fail" % (ok, bad))

    ok, bad = step3()
    print("STEP 3  floor(t/2)+t+C(t+3,2) = 1 mod 2:   %4d ok, %s" % (ok, bad or "no exception"))

    ok, bad = step2(order_by_residue=False)
    print("CONTROL same, ordering dropped:            %4d ok, %d fail  <-- must fail" % (ok, bad))
