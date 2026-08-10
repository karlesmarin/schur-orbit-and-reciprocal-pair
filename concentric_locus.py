# -*- coding: utf-8 -*-
"""The concentric branch of the vanishing criterion, in quotient coordinates (the Question on the
vanishing locus as an arrangement).

Corollary 3.2 says Phi_t vanishes when a residue class is empty or when d3 = 0.  Writing
a1 + a2 = t*(S_A) + 2*rA with S_A = |lam^(rA)| + 1, the second condition reads

        t*(s_A - s_B) = 2*(rB - rA),        s_r = |lam^(r)| ,

so it is a single hyperplane in the QUOTIENT-SIZE vector, and it has integer points only when
t divides 2*(rB - rA).  With 0 <= rA < rB <= t-1 that forces rB - rA = t/2, hence t even; and
then s_A = s_B + 1.  Three things are checked:

  CLAIM 1   for t ODD the concentric branch is empty -- Phi_t vanishes only by an empty class.
            The control is the even-t count on the same range: it must be nonzero, or CLAIM 1
            is vacuous rather than true.

  CLAIM 2   for t EVEN, d3 = 0  <==>  rB = rA + t/2  and  |lam^(rA)| = |lam^(rB)| + 1,
            tested in both directions.

  CLAIM 3   d3 is NOT a function of the residue profile: two shapes with the same profile (hence
            the same t-core, hence the same point of the root lattice under Garvan-Kim-Stanton)
            can carry different d3.  So the locus is an arrangement in the quotient sizes, not in
            the core lattice.

Authors: Carles Marin, Claude (AI assistant)."""
import contextlib
import io
from collections import defaultdict

from law_control import partitions
from theorem_full import setup

with contextlib.redirect_stdout(io.StringIO()):   # extra_structure runs its own sweep on import
    from extra_structure import quot

TMAX, NMAX = 10, 17


def data(lam, t):
    """(rA, rB, sA, sB, d3, overlap) with the blocks ordered by residue, or None."""
    st = setup(lam, t)
    if st is None:
        return None
    beta, Ac, Bc = st
    rA, rB = beta[Ac[0]] % t, beta[Bc[0]] % t
    if rA > rB:
        Ac, Bc = Bc, Ac
        rA, rB = rB, rA
    a1, a2 = beta[Ac[0]], beta[Ac[1]]
    b1, b2 = beta[Bc[0]], beta[Bc[1]]
    N = t + 2
    q = quot(lam, t, N)
    return rA, rB, sum(q[rA]), sum(q[rB]), abs(a1 + a2 - b1 - b2), len(set(Ac) & set(Bc)) > 0


def sweep():
    odd_conc = even_conc = 0
    c2_ok = c2_bad = 0
    profiles = defaultdict(set)
    shapes = 0
    for t in range(2, TMAX):
        N = t + 2
        for n in range(0, NMAX):
            for lam in partitions(n, N):
                d = data(lam, t)
                if d is None:
                    continue
                rA, rB, sA, sB, d3, overlap = d
                shapes += 1
                if d3 == 0:
                    (even_conc, odd_conc) = (even_conc + 1, odd_conc) if t % 2 == 0 \
                        else (even_conc, odd_conc + 1)
                if t % 2 == 0:
                    lhs = (d3 == 0)
                    rhs = (not overlap) and (rB - rA == t // 2) and (sA == sB + 1)
                    if lhs == rhs:
                        c2_ok += 1
                    else:
                        c2_bad += 1
                st = setup(lam, t)
                prof = tuple(sorted((st[0][j] % t) for j in range(N)))
                profiles[(t, prof)].add(d3)
    return shapes, odd_conc, even_conc, c2_ok, c2_bad, profiles


if __name__ == "__main__":
    shapes, odd_c, even_c, ok, bad, profiles = sweep()
    print("range: t = 2..%d, |lambda| <= %d, %d shapes with a non-degenerate profile"
          % (TMAX - 1, NMAX - 1, shapes))
    print()
    print("CLAIM 1  concentric (d3 = 0) at ODD t:      %d      <-- must be 0" % odd_c)
    print("CONTROL  concentric (d3 = 0) at EVEN t:     %d      <-- must NOT be 0" % even_c)
    print()
    print("CLAIM 2  t even, d3=0 <=> rB-rA=t/2 and sA=sB+1:  %d ok, %d fail" % (ok, bad))
    print()
    split = {k: v for k, v in profiles.items() if len(v) > 1}
    print("CLAIM 3  residue profiles carrying more than one value of d3: %d of %d"
          % (len(split), len(profiles)))
    if split:
        k = sorted(split, key=lambda k: (k[0], len(k[1])))[0]
        print("         smallest witness: t = %d, profile %s carries d3 in %s"
              % (k[0], list(k[1]), sorted(split[k])))
        print("         => d3 is not a function of the residue vector, so the locus is an")
        print("            arrangement in the quotient sizes and not in the core lattice.")
