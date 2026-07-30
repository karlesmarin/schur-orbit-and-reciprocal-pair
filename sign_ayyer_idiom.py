# -*- coding: utf-8 -*-
"""Is the derived sign expressible in Ayyer-Kumari's own idiom (their sgn(sigma), eq. 2.1)?

This is the script behind the Question on the short form of eps_lambda, and it prints three things:

  CONTROL   the same comparison with the two blocks taken in whatever order `setup` returns.
            It FAILS, and it is supposed to: lambda_11 is antisymmetric under A <-> B (through its
            sgn(a1-b1) factor) and so is sgn(a1+a2-b1-b2), so their product eps_lambda is
            swap-invariant, while the candidate formula is not -- (-1)^(r_A+r_B) is symmetric, so only
            its last factor changes sign. Without a rule fixing which block is A the right-hand side
            is not a function of lambda at all. This control is what makes the ordering hypothesis
            load-bearing rather than decorative.

  TEST 1    the displayed closed form with A and B ordered so that r_A <= r_B.

  TEST 2    the weaker structural claim: with the blocks so ordered, is the residual
            eps_lambda * sgn(sigma) * sgn(a1+a2-b1-b2) a function of (t, r_A, r_B) alone?
            Equivalently, does any cell (t,r_A,r_B) carry both signs?

Authors: Carles Marin, Claude (AI assistant)."""
from mpmath import mp
from law_control import partitions
from theorem_full import setup, lambda11
mp.dps = 30

TMAX, NMAX = 9, 15          # t = 2..8, |lambda| <= 14


def perm_sign(seq):
    seq = list(seq); s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def sgn_sigma(beta, t):
    """their (2.1): regroup beta by residue class, classes increasing, values decreasing."""
    order = []
    for a in range(t):
        order += [j for j, b in enumerate(beta) if b % t == a]
    return perm_sign(order)


def run(order_by_residue):
    """compare eps_lambda with the short form; return (agree, disagree, cells)."""
    agree = disagree = 0
    cells = {}
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
                l11 = lambda11(beta, Ac, Bc, N)     # recomputed AFTER any swap: it is antisymmetric
                if l11 is None:
                    continue
                a1, a2 = beta[Ac[0]], beta[Ac[1]]
                b1, b2 = beta[Bc[0]], beta[Bc[1]]
                d = a1 + a2 - b1 - b2
                if d == 0:
                    continue
                sgn_d = 1 if d > 0 else -1
                eps = (-1) ** (t + (t + 2) * (t + 3) // 2) * l11 * sgn_d
                theirs = ((-1) ** (t // 2)) * sgn_sigma(beta, t) * ((-1) ** (rA + rB)) * sgn_d
                if eps == theirs:
                    agree += 1
                else:
                    disagree += 1
                cells.setdefault((t, rA, rB), set()).add(eps * sgn_sigma(beta, t) * sgn_d)
    return agree, disagree, cells


a0, d0, _ = run(order_by_residue=False)
print("CONTROL -- blocks in arbitrary order (the ordering hypothesis dropped):")
print("   agree %d   disagree %d   <-- must fail, and does" % (a0, d0))

a1, d1, cells = run(order_by_residue=True)
print("\nTEST 1 -- the displayed closed form, A and B ordered so that r_A <= r_B:")
print("   agree %d   disagree %d" % (a1, d1))

both = [k for k, v in cells.items() if len(v) == 2]
print("\nTEST 2 -- is the residual a function of (t,r_A,r_B)?")
print("   cells occupied: %d ; cells carrying BOTH signs: %d" % (len(cells), len(both)))
if both:
    print("   offending cells:", sorted(both)[:12])
else:
    print("   => every cell is sign-homogeneous, so once the blocks are ordered the sign")
    print("      depends on t, r_A and r_B alone.")

print("\nrange: t = 2..%d, |lambda| <= %d" % (TMAX - 1, NMAX - 1))
