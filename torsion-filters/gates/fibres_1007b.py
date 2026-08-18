"""Problem 10.7, second pass: the fibre as a LATTICE, with the paper's own (A,B) convention.

The first pass rebuilt 72/110, 36/49, 24/30 fibres and missed the rest.  Every mismatch was a
labelling one: theorem_full.setup() calls A the class holding the LARGEST beta number (it sorts the
big classes by column index, and column 0 carries the largest beta), so a parametrisation that loops
over ordered pairs (rA,rB) generates each configuration twice with the two labellings, and buckets
half of them under the transposed triple.

So this pass does not build one fibre at a time.  It generates EVERY two-class configuration from
the free parameters, labels each one by the paper's rule, buckets it, and compares the whole
bucketing against brute force -- dictionary against dictionary, not fibre by fibre.  That is a
stronger test: it also catches a parametrisation that is complete but produces spurious extras.

THE PARAMETRISATION (N = t+2 beta numbers, t residue classes, two-class profile):
    residues       rA != rB carry two beta numbers each; the other t-2 carry one each
    component A    (v + kA, v)   ->  beta  t(v+kA+1)+rA,  t*v+rA          v >= 0, kA >= 0
    component B    (u + kB, u)   ->  beta  t(u+kB+1)+rB,  t*u+rB          u >= 0, kB >= 0
    free classes   one part m_i >= 0 each  ->  beta  t*m_i + r_i
and |lambda| is then determined by the beta-set.  Four integer parameters (kA, kB, v, u) plus t-2
free ones, with no constraint tying them: the whole two-class stratum is that lattice.

If the bucketing agrees, then the fibre over a triple is an explicit lattice section and its
generating function follows by inspection, which is what Problem 10.7 asks for.

CONTROLS:
  C1  dictionary equality both ways -- triples present, and the set behind each triple.
  C2  a WRONG labelling (A = smallest beta instead of largest) must disagree.  If both labellings
      agree, the test is not seeing the convention and proves nothing about it.
  C3  the count of partitions covered must equal the brute-force count of two-class partitions.

Authors: Carles Marin, Claude (AI assistant).
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "paper", "anc"))

from law_control import partitions          # noqa: E402
from theorem_full import setup              # noqa: E402


def lam_from_beta(b, N):
    b = sorted(b, reverse=True)
    lam = [b[j] - (N - 1 - j) for j in range(N)]
    if any(x < 0 for x in lam):
        return None
    while lam and lam[-1] == 0:
        lam.pop()
    return tuple(lam)


def brute(t, nmax):
    """triple -> set of partitions, from the paper's own setup()"""
    N = t + 2
    fib = defaultdict(set)
    other = 0
    for n in range(0, nmax + 1):
        for lam in partitions(n, N):
            st = setup(lam, t)
            if st is None:
                other += 1
                continue
            beta, Ac, Bc = st
            if len(set(Ac) | set(Bc)) != 4:      # size-three profile shares a column
                other += 1
                continue
            a1, a2 = beta[Ac[0]], beta[Ac[1]]
            b1, b2 = beta[Bc[0]], beta[Bc[1]]
            fib[(a1 - a2, b1 - b2, abs(a1 + a2 - b1 - b2))].add(tuple(lam))
    return fib, other


def built(t, nmax, wrong_label=False):
    """every two-class configuration, generated from the lattice and labelled by the paper's rule"""
    N = t + 2
    fib = defaultdict(set)
    K = nmax // t + 3                                    # a generous cap on every free parameter
    for rA in range(t):
        for rB in range(rA + 1, t):                      # unordered; the label is decided below
            frees = [r for r in range(t) if r not in (rA, rB)]

            def rec(i, acc, budget):
                if i == len(frees):
                    yield acc
                    return
                for m in range(0, budget + 1):
                    yield from rec(i + 1, acc + [t * m + frees[i]], budget - m)

            for kA in range(0, K):
                for v in range(0, K):
                    pA = [t * (v + kA + 1) + rA, t * v + rA]
                    if pA[0] > t * K + N:
                        break
                    for kB in range(0, K):
                        for u in range(0, K):
                            pB = [t * (u + kB + 1) + rB, t * u + rB]
                            if pB[0] > t * K + N:
                                break
                            for rest in rec(0, [], K):
                                b = pA + pB + rest
                                if len(set(b)) != N:
                                    continue
                                lam = lam_from_beta(b, N)
                                if lam is None or sum(lam) > nmax:
                                    continue
                                # the paper's rule: A is the class of the LARGEST beta number
                                first, second = (pA, pB) if pA[0] > pB[0] else (pB, pA)
                                if wrong_label:
                                    first, second = second, first
                                d = (first[0] - first[1], second[0] - second[1],
                                     abs(sum(first) - sum(second)))
                                fib[d].add(lam)
    return fib


def compare(t, nmax):
    B, other = brute(t, nmax)
    for wrong in (False, True):
        G = built(t, nmax, wrong_label=wrong)
        same = (dict(B) == dict(G))
        nb = sum(len(v) for v in B.values())
        ng = len(set().union(*G.values())) if G else 0
        tag = "  [CONTROL: WRONG LABELLING]" if wrong else ""
        print("t = %d,  |lambda| <= %d%s" % (t, nmax, tag))
        print("    brute: %d partitions on %d triples   (%d in other profiles)"
              % (nb, len(B), other))
        print("    built: %d partitions on %d triples" % (ng, len(G)))
        print("    dictionaries equal: %s" % ("YES" if same else "no"))
        if not same:
            missing = [d for d in B if d not in G]
            extra = [d for d in G if d not in B]
            diff = [d for d in B if d in G and B[d] != G[d]]
            print("        triples only in brute: %d %s" % (len(missing), missing[:4]))
            print("        triples only in built: %d %s" % (len(extra), extra[:4]))
            print("        triples present in both but different sets: %d %s"
                  % (len(diff), diff[:4]))
            for d in diff[:2]:
                print("            d = %s  brute-not-built %s  built-not-brute %s"
                      % (d, sorted(B[d] - G[d])[:3], sorted(G[d] - B[d])[:3]))
        print()


if __name__ == "__main__":
    print("=" * 96)
    print("Problem 10.7: the two-class stratum as an explicit lattice")
    print("=" * 96)
    print()
    for t, nmax in [(2, 14), (3, 12), (4, 10)]:
        compare(t, nmax)
