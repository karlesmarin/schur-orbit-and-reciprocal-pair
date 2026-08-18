"""Problem 10.7 -- the fibres of the evaluation invariant.

The paper asks: "How many partitions share a given (d1,d2,d3), and which?  ...  A generating
function for the fibre sizes would say how much of a partition this alphabet is blind to."

It also says where to look: by Proposition prop:quotient the fibre is the set of t-quotients with
prescribed sl2 contents and a prescribed difference of sizes.  This file takes that literally and
asks whether the fibre can be BUILT from the quotient side rather than searched for.

THE CLAIM UNDER TEST (two-class profile, N = t+2 beta numbers over t residue classes):

  * the profile is: two classes rA != rB carry two beta numbers each, the other t-2 carry one each;
  * core_t(lambda) is then determined by the profile alone (Remark 8.7), so it is CONSTANT on the
    part of a fibre with a given profile;
  * lambda^(rA) and lambda^(rB) have two parts; d1, d2 fix their sl2 contents
    kA = d1/t - 1, kB = d2/t - 1, so each component is determined by its SIZE alone:
        lambda^(rA) = (v + kA, v),   sA = 2v + kA,      v >= 0
        lambda^(rB) = (u + kB, u),   sB = 2u + kB,      u >= 0
  * the other t-2 components carry one beta number each, so each is a single part m_i >= 0, FREE;
  * d3 = |t(sA - sB) + 2(rA - rB)| ties v - u to a constant;
  * |lambda| = |core| + t*(sA + sB + sum m_i).

If that is right then a fibre is a lattice of the shape (one free parameter v, with u = v - const)
times (t-2 free nonnegative integers), every point of it sitting at a size in an explicit
progression -- and its generating function is forced:

        sum over (profile, branch) of    q^{c} / [ (1 - q^{4t}) * (1 - q^{t})^{t-2} ]

with 4t from v moving both sA and sB by 2 each, and one factor 1 - q^t per free component.

THE TEST is set equality, not a count: build every fibre from the quotient side inside a size
bound, build it again by brute force from the beta-sets, and compare AS SETS.  A count could agree
by accident; the sets cannot.

CONTROLS, each able to fail:
  C1  a WRONG denominator -- predict with (1-q^{2t}) in place of (1-q^{4t}) -- must disagree.
  C2  the constructed fibre must be checked against the brute-force one for EVERY triple, not a
      sample, and the script prints the number of triples compared.
  C3  the size-three and degenerate profiles are reported separately and NOT counted as agreement;
      the claim above is about the two-class profile only.

Authors: Carles Marin, Claude (AI assistant).
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "paper", "anc"))

from law_control import partitions          # noqa: E402
from theorem_full import setup              # noqa: E402
from extra_structure import quot            # noqa: E402


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - (j + 1) for j in range(N)]


def profile(lam, t, N):
    """how many beta numbers in each residue class"""
    b = beta_of(lam, N)
    p = [0] * t
    for v in b:
        p[v % t] += 1
    return tuple(p)


def lam_from_beta(b, N):
    """the partition whose beta-set (N numbers) is b"""
    b = sorted(b, reverse=True)
    lam = [b[j] - (N - (j + 1)) for j in range(N)]
    assert all(x >= 0 for x in lam), b
    while lam and lam[-1] == 0:
        lam.pop()
    return tuple(lam)


def triple(lam, t):
    """(d1,d2,d3) via the paper's own setup(); None if the profile is not two-class"""
    st = setup(lam, t)
    if st is None:
        return None
    beta, Ac, Bc = st
    a1, a2 = beta[Ac[0]], beta[Ac[1]]
    b1, b2 = beta[Bc[0]], beta[Bc[1]]
    if len(Ac) != 2 or len(Bc) != 2:
        return None
    return (a1 - a2, b1 - b2, abs(a1 + a2 - b1 - b2))


def brute_fibres(t, nmax):
    """triple -> set of partitions, by brute force over |lambda| <= nmax"""
    N = t + 2
    fib = defaultdict(set)
    other = 0
    for n in range(0, nmax + 1):
        for lam in partitions(n, N):
            d = triple(lam, t)
            if d is None:
                other += 1
                continue
            fib[d].add(tuple(lam))
    return fib, other


def built_fibre(t, d, nmax):
    """the fibre over d = (d1,d2,d3), BUILT from the quotient description, inside |lambda| <= nmax.

    Constructed by walking the free parameters directly: choose the two distinguished residues,
    choose v and u with the d3 constraint, choose the t-2 free single parts, reassemble the
    beta-set, and read off lambda.  Nothing here looks at the brute-force answer."""
    N = t + 2
    d1, d2, d3 = d
    out = set()
    if d1 % t or d2 % t:
        return out
    kA, kB = d1 // t - 1, d2 // t - 1
    if kA < 0 or kB < 0:
        return out
    for rA in range(t):
        for rB in range(t):
            if rA == rB:
                continue
            for v in range(0, nmax // max(t, 1) + 2):
                sA = 2 * v + kA
                for sB_sign in (+1, -1):
                    # t(sA - sB) + 2(rA - rB) = sB_sign * d3
                    num = sB_sign * d3 - 2 * (rA - rB)
                    if num % t:
                        continue
                    sB = sA - num // t
                    if sB < kB or (sB - kB) % 2:
                        continue
                    u = (sB - kB) // 2
                    # the two distinguished components, as two beta numbers each
                    # lambda^(rA) = (v+kA, v) -> its two beta numbers within the class
                    aa = [t * (v + kA + 1) + rA, t * v + rA]
                    bb = [t * (u + kB + 1) + rB, t * u + rB]
                    frees = [r for r in range(t) if r not in (rA, rB)]
                    # the t-2 free single parts m_i >= 0, one per remaining class
                    def rec(i, acc, budget):
                        if i == len(frees):
                            b = aa + bb + acc
                            if len(set(b)) != N:
                                return
                            lam = lam_from_beta(b, N)
                            if sum(lam) <= nmax:
                                out.add(lam)
                            return
                        for m in range(0, budget + 1):
                            rec(i + 1, acc + [t * m + frees[i]], budget - m)
                    rec(0, [], nmax // t + 2)
    return out


def main():
    print("=" * 96)
    print("Problem 10.7: is the fibre BUILDABLE from the quotient description?")
    print("=" * 96)
    print()
    for t, nmax in [(2, 16), (3, 14), (4, 12)]:
        fib, other = brute_fibres(t, nmax)
        agree = disagree = 0
        firstbad = []
        for d in sorted(fib):
            built = built_fibre(t, d, nmax)
            if built == fib[d]:
                agree += 1
            else:
                disagree += 1
                if len(firstbad) < 3:
                    firstbad.append((d, sorted(fib[d] - built)[:3],
                                     sorted(built - fib[d])[:3]))
        npart = sum(len(v) for v in fib.values())
        print("t = %d,  |lambda| <= %d" % (t, nmax))
        print("    %d partitions in two-class profile, on %d triples "
              "(%d partitions in other profiles, not claimed)" % (npart, len(fib), other))
        print("    fibres rebuilt exactly: %d / %d" % (agree, agree + disagree))
        for d, missing, extra in firstbad:
            print("        d = %s   brute-not-built %s   built-not-brute %s" % (d, missing, extra))
        print()


if __name__ == "__main__":
    main()
