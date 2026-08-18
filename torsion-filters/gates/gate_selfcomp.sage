# -*- coding: utf-8 -*-
# THE GATE THAT CAN KILL IT: is the criterion just branch (b) in disguise?
#
# The criterion measured today is
#     Phi_t == 0  <=>  (i) the excess values S are concentric, C - S = S, C = min S + max S,
#                      (ii) both residues k with 2k = C (mod t) are excess classes.
# The five r=3 zeros printed by hand are ALL self-complementary of odd width -- which is branch (b)
# of Theorem 9.1 of our own arXiv:2608.09619, sufficient by a two-line argument that needs only an
# inversion-closed alphabet with product -1 (true for every even t).  If (i)+(ii) is equivalent to
# branch (b) then nothing new has been found in the interior, and the honest verdict is that the
# locus there is the classical one and only the CONVERSE is open.
#
# The distinction is structural, not cosmetic.  If every residue class has size >= 2 then the
# excess values ARE the whole beta set, concentricity becomes full self-complementarity, and
# C even <=> w odd, so (i)+(ii) collapses to branch (b).  The criterion can only be strictly
# stronger through shapes with a SINGLETON class, where it constrains the excess block alone.  At
# r = 1 those are the generic case -- Corollary 3.2's d3 = 0 is a condition on 4 of the t+2 beta
# numbers, not on lambda -- so at r = 1 the two are certainly different.  The question is r >= 2.
#
# WHAT IS MEASURED, per (t, r): the 2 x 2 table  criterion x self-complementary-of-odd-width, over
# the shapes with every residue class occupied (branch (a) is the empty-class branch and is out of
# range by construction).  Every zero that is NOT self-complementary is printed by hand, and so is
# every self-complementary shape that does not vanish.
#
# ANCHOR, fatal: at t = 2 the criterion must reproduce Conjecture 9.4 of the paper (branch (b)),
# which was verified there over its own ranges.  A disagreement at t = 2 means the criterion is
# wrong, not new.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

PRIME = 2013265921
F = GF(PRIME)
G = F.multiplicative_generator()


def roots_of(t):
    z = G ** ((PRIME - 1) // t)
    return [z ** k for k in range(t)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def classes(beta, t):
    cls = {}
    for k, b in enumerate(beta):
        cls.setdefault(b % t, []).append(k)
    return cls


def criterion(beta, cls, t):
    S = set(beta[i] for v in cls.values() if len(v) >= 2 for i in v)
    if not S:
        return False, None
    C = min(S) + max(S)
    if set(C - b for b in S) != S:
        return False, None
    exc = [k for k in cls if len(cls[k]) >= 2]
    s = sum(1 for k in exc if (2 * k - C) % t == 0)
    return s == 2, C


def selfcomp_odd(lam, N):
    """lambda_i + lambda_{N+1-i} = w for all i, with w odd."""
    l = list(lam) + [0] * (N - len(lam))
    w = l[0] + l[N - 1]
    if w % 2 == 0:
        return False
    return all(l[i] + l[N - 1 - i] == w for i in range(N))


print("=" * 94)
print("IS THE CRITERION BRANCH (b) IN DISGUISE?   table: criterion x self-complementary odd width")
print("=" * 94)
print("")
print("   t   r   N |lam|<=  shapes  ZEROS   crit&SC  crit&!SC  !crit&SC  singleton-class zeros")
print("  " + "-" * 90)

CASES = ((2, 1, 30), (2, 2, 26), (2, 3, 22), (2, 4, 20),
         (4, 1, 34), (4, 2, 30), (4, 3, 26), (4, 4, 22),
         (6, 2, 22), (6, 3, 18), (8, 2, 20))
ODD = []
for t, r, MAX in CASES:
    N = t + 2 * r
    EMAX = MAX + N
    RT = roots_of(t)
    set_random_seed(313 + 10 * t + r)
    pts = []
    while len(pts) < 3:
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        al = RT + [y for x in zz for y in (x, 1 / x)]
        if len(set(al)) == N:
            pts.append(al)
    POW = [[[a ** e for e in range(EMAX + 1)] for a in al] for al in pts]

    nsh = zeros = c_sc = c_nsc = nc_sc = sing = bad = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            lam = list(l)
            beta = beta_of(lam, N)
            cls = classes(beta, t)
            if len(cls) < t:
                continue
            nsh += 1
            v = True
            for p in range(3):
                if matrix(F, [[POW[p][b][beta[a]] for b in range(N)]
                              for a in range(N)]).det() != 0:
                    v = False
                    break
            crit, C = criterion(beta, cls, t)
            if crit != v:
                bad += 1
            sc = selfcomp_odd(lam, N)
            if v:
                zeros += 1
                if min(len(vv) for vv in cls.values()) == 1:
                    sing += 1
                if not sc:
                    ODD.append((t, r, lam, beta, sorted(len(vv) for vv in cls.values()), C))
            if crit and sc:
                c_sc += 1
            elif crit and not sc:
                c_nsc += 1
            elif sc and not crit:
                nc_sc += 1
    print("  %2d %3d %3d %6d %8d %6d %8d %9d %9d %10d%s"
          % (t, r, N, MAX, nsh, zeros, c_sc, c_nsc, nc_sc, sing,
             "   <-- CRIT!=DET on %d" % bad if bad else ""))

print("")
print("  'crit&!SC' is the column that decides.  If it is 0 everywhere for r >= 2, the interior")
print("  locus is exactly branch (b) and the criterion adds nothing there but the converse.")
print("  '!crit&SC' must be 0: branch (b) is PROVED sufficient, so a self-complementary shape of")
print("  odd width that the criterion misses would refute the criterion, not the branch.")
print("")
print("  every zero that is NOT self-complementary of odd width, by hand (%d of them):" % len(ODD))
for t, r, lam, beta, sizes, C in ODD[:25]:
    print("    t=%d r=%d  lam=%-26s beta=%s  sizes=%s  C=%s" % (t, r, lam, beta, sizes, C))
if len(ODD) > 25:
    print("    ... and %d more" % (len(ODD) - 25))
print("")
print("DONE")
