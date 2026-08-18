# -*- coding: utf-8 -*-
# *** THE CLOSING PROSE OF L3 IS WRONG AND IS LEFT STANDING.  It predicts that at r = 1 the two
# *** last columns coincide, because "concentricity forces s = 2 there".  The run says 2212 of
# *** 2436.  The reason is the case Corollary 3.2 excludes by hypothesis: at r = 1 the excess can
# *** also be ONE class of three, which is often concentric and never vanishes, and it is exactly
# *** condition (ii) that kills it.  So (ii) is doing work already at r = 1.  The line was written
# *** before the run; it is a refuted prediction and it stays visible instead of being edited away.
#
# The r = 1 anchor: does the new criterion reproduce the PUBLISHED Corollary 3.2?
#
# The criterion measured today is
#
#     Phi_t(lambda; z) == 0   <=>   (i) the excess values S are concentric, C - S = S with
#                                       C = min S + max S, and
#                                  (ii) BOTH residues k with 2k = C (mod t) are excess classes.
#
# At r = 1 the answer is already in print: Corollary 3.2 of arXiv:2608.09619 says Phi_t vanishes
# iff a residue class is empty or d3 = 0, and d3 = |a1 + a2 - b1 - b2| on the two excess classes.
# If the new criterion is right it must agree with that on the nose -- and if it does not, every
# number produced today is wrong, because Corollary 3.2 is refereed and this is not.
#
# Three predicates are compared on the same shapes:
#   DET    the bialternant over GF(p), the ground truth;
#   NEW    (i) and (ii) above;
#   PAPER  d3 = 0 on the two excess classes (shapes with an empty class are excluded from the
#          sweep, so that half of Corollary 3.2 is not at issue).
#
# CONTROLS able to fail:
#   L1  DET vs NEW, DET vs PAPER, NEW vs PAPER: every disagreement printed by hand.
#   L2  the degenerate profiles -- one excess class of size 3, which Corollary 3.2 does not cover
#       because setup() returns None on them -- are counted and tested separately.  NEW predicts
#       no zeros there (s = 1), and that is a prediction the published statement does not make.
#   L3  a control that can fail: the same tally with condition (ii) DROPPED, i.e. concentricity
#       alone.  At r = 1 it must agree too (there (ii) is automatic); at r >= 2 it must NOT, or
#       (ii) is idle and the criterion is just concentricity.
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


def centre(beta, cls):
    S = set(beta[i] for v in cls.values() if len(v) >= 2 for i in v)
    if not S:
        return None
    C = min(S) + max(S)
    return C if set(C - b for b in S) == S else None


def new_criterion(beta, cls, t):
    """((i) and (ii),  (i) alone)"""
    C = centre(beta, cls)
    if C is None:
        return False, False
    exc = [k for k in cls if len(cls[k]) >= 2]
    s = sum(1 for k in exc if (2 * k - C) % t == 0)
    return s == 2, True


def paper_d3(beta, cls):
    """Corollary 3.2's d3, or None on a degenerate profile."""
    exc = [k for k in cls if len(cls[k]) >= 2]
    if len(exc) != 2 or any(len(cls[k]) != 2 for k in exc):
        return None
    (a1, a2) = [beta[i] for i in cls[exc[0]]]
    (b1, b2) = [beta[i] for i in cls[exc[1]]]
    return abs(a1 + a2 - b1 - b2)


print("=" * 92)
print("L1  the r = 1 anchor:  determinant  vs  the new criterion  vs  Corollary 3.2")
print("=" * 92)
print("")
print("   t   N |lam|<=  shapes  ZEROS   DET!=NEW   DET!=PAPER   degenerate  zeros among them")
print("  " + "-" * 88)

for t, MAX in ((4, 40), (6, 30), (8, 26), (10, 22), (3, 30), (5, 24)):
    r = 1
    N = t + 2
    EMAX = MAX + N
    RT = roots_of(t)
    set_random_seed(900 + t)
    pts = []
    while len(pts) < 3:
        z = F.random_element()
        if z == 0:
            continue
        al = RT + [z, 1 / z]
        if len(set(al)) == N:
            pts.append(al)
    POW = [[[a ** e for e in range(EMAX + 1)] for a in al] for al in pts]

    nsh = zeros = bad_new = bad_paper = ndeg = zdeg = bad_conc = 0
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
            if v:
                zeros += 1
            nw, conc = new_criterion(beta, cls, t)
            if nw != v:
                bad_new += 1
                if bad_new <= 3:
                    print("     DET!=NEW  t=%d lam=%s beta=%s det0=%s new=%s" % (t, lam, beta, v, nw))
            if conc != v:
                bad_conc += 1
            d3 = paper_d3(beta, cls)
            if d3 is None:
                ndeg += 1
                if v:
                    zdeg += 1
            else:
                if (d3 == 0) != v:
                    bad_paper += 1
                    if bad_paper <= 3:
                        print("     DET!=PAPER t=%d lam=%s beta=%s det0=%s d3=%s"
                              % (t, lam, beta, v, d3))
    print("  %2d %3d %6d %8d %6d %10d %12d %12d %17d"
          % (t, N, MAX, nsh, zeros, bad_new, bad_paper, ndeg, zdeg))

print("")
print("  L2  'degenerate' = one excess class of size 3, which Corollary 3.2 does not cover.")
print("      The new criterion predicts NO zeros there (s = 1); the last column is the test.")

print("")
print("=" * 92)
print("L3  the control: is condition (ii) doing any work?  Same tally with (i) ALONE")
print("=" * 92)
print("")
print("   t   r   N |lam|<=  concentric  ZEROS  (i)alone wrong on   (i)+(ii) wrong on")
print("  " + "-" * 84)
for t, r, MAX in ((4, 1, 34), (4, 2, 30), (4, 3, 26), (6, 2, 20), (6, 3, 18)):
    N = t + 2 * r
    EMAX = MAX + N
    RT = roots_of(t)
    set_random_seed(700 + 10 * t + r)
    pts = []
    while len(pts) < 3:
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        al = RT + [y for x in zz for y in (x, 1 / x)]
        if len(set(al)) == N:
            pts.append(al)
    POW = [[[a ** e for e in range(EMAX + 1)] for a in al] for al in pts]
    ncon = zeros = w1 = w2 = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            cls = classes(beta, t)
            if len(cls) < t:
                continue
            v = True
            for p in range(3):
                if matrix(F, [[POW[p][b][beta[a]] for b in range(N)]
                              for a in range(N)]).det() != 0:
                    v = False
                    break
            zeros += 1 if v else 0
            nw, conc = new_criterion(beta, cls, t)
            if conc:
                ncon += 1
            if conc != v:
                w1 += 1
            if nw != v:
                w2 += 1
    print("  %2d %3d %3d %6d %11d %6d %19d %19d" % (t, r, N, MAX, ncon, zeros, w1, w2))

print("")
print("  If (ii) were idle the two last columns would be equal.  At r = 1 they are (concentricity")
print("  forces s = 2 there, which is why Corollary 3.2 needs only d3); from r = 2 they separate.")
print("")
print("DONE")
