# -*- coding: utf-8 -*-
# The reflection sign q, in closed form, and the criterion on a big sample.
#
# parity_terms.sage measured, with no sign convention anywhere, the ratio q = term(C-T)/term(T) of
# the Laplace expansion.  It is +-1, constant over the terms, independent of the point, and on
# 2138 concentric shapes at r = 1,2,3,4 it satisfied
#
#       Phi vanishes identically   <=>   beta concentric  and  q = -1
#
# with no exception either way.  There is NO r-parity in it: the criterion is the same sign for
# every r.  The alternation reported earlier came from an eps that omitted the sign of the
# root-of-unity minor (it is +-V, and which sign depends on the order of the picked residues).
#
# This script does two things.
#
# (1) ROUTE B.  From the diagnosis, the ratio must be
#         q  =  (-1)^{sum(pick) + sum(pick2)} * sigma(pick) * sigma(pick2),
#     sigma(P) = sign of the residue word (beta_i mod t, i in P sorted by index) as a permutation
#     of {0..t-1}.  This uses no field arithmetic.  It is checked term by term against route A,
#     the measured field ratio, on every concentric shape of the acceptance range.  Disagreement
#     is fatal.
#
# (2) THE WIDE RUN.  With route B the criterion costs nothing, so r = 3 is pushed to |lambda| <= 50
#     -- from 6 zeros this morning to several hundred -- and the biconditional is re-tested there.
#     Then q is tabulated against coarse invariants of beta, to see what a closed form could
#     possibly depend on.
#
# CONTROLS able to fail:
#   C1  route A vs route B, term by term.  Must be 0 disagreements.
#   C2  q must be constant over the terms of a shape.  Non-constant shapes are counted, not hidden.
#   C3  a zero that is not concentric.  Must be 0.
#   C4  q = -1 and not vanishing, or vanishing and q != -1.  Either one kills the criterion.
#   C5  a control that CAN fail: the same tally for a deliberately wrong sign, q' = q * (-1)^r,
#       which is what the retracted reading predicted.  It must do badly, or this run proves
#       nothing about which of the two is right.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools
import time

t = 4
PRIME = 998244353
F = GF(PRIME)
I4 = F(-1).sqrt()
ROOTS = [F(1), I4, F(-1), -I4]
SGN_T = (-1) ** (t * (t - 1) // 2)


def make_points(r, N, seed, howmany):
    set_random_seed(seed)
    pts = []
    while len(pts) < howmany:
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        alph = ROOTS + [y for x in zz for y in (x, 1 / x)]
        if len(set(alph)) == N:
            pts.append(alph)
    return pts


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def classes(beta):
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


def word_sign(seq):
    """sign of a sequence of distinct values, as a permutation."""
    sg = 1
    s = list(seq)
    for a in range(len(s)):
        for b in range(a + 1, len(s)):
            if s[a] > s[b]:
                sg = -sg
    return sg


def q_route_B(beta, cls, C, N, all_terms=True):
    """(q, constant?) with no field arithmetic."""
    pos = dict((b, i) for i, b in enumerate(beta))
    iota = {}
    for i, b in enumerate(beta):
        j = pos.get(C - b)
        if j is not None:
            iota[i] = j
    vals = set()
    keys = sorted(cls)
    picks = itertools.product(*[cls[k] for k in keys]) if all_terms \
        else [tuple(cls[k][0] for k in keys)]
    for pick in picks:
        pick = tuple(sorted(pick))
        T = tuple(sorted(set(range(N)) - set(pick)))
        if any(i not in iota for i in T):
            return None, False
        T2 = tuple(sorted(iota[i] for i in T))
        pick2 = tuple(sorted(set(range(N)) - set(T2)))
        q = (-1) ** (sum(pick) + sum(pick2)) \
            * word_sign([beta[i] % t for i in pick]) \
            * word_sign([beta[i] % t for i in pick2])
        vals.add(q)
        if len(vals) > 1:
            return 0, False
    return vals.pop(), True


def q_route_A(beta, cls, alph, C, N):
    """the measured field ratio, term by term."""
    pos = dict((b, i) for i, b in enumerate(beta))
    iota = dict((i, pos[C - b]) for i, b in enumerate(beta) if (C - b) in pos)
    M = matrix(F, [[alph[b] ** beta[a] for b in range(N)] for a in range(N)])
    tm = {}
    for pick in itertools.product(*[cls[k] for k in sorted(cls)]):
        pick = tuple(sorted(pick))
        T = tuple(sorted(set(range(N)) - set(pick)))
        tm[T] = (-1) ** sum(pick) * SGN_T \
            * M[list(pick), list(range(t))].det() * M[list(T), list(range(t, N))].det()
    if sum(tm.values()) != M.det():
        return "SUMFAIL"
    vals = set()
    for T, v in tm.items():
        T2 = tuple(sorted(iota[i] for i in T))
        if v == 0:
            continue
        r = tm[T2] / v
        vals.add(1 if r == F(1) else (-1 if r == F(-1) else 99))
    if len(vals) != 1:
        return 0
    return vals.pop()


# =============================================================== C1: route A vs route B, fatal
print("=" * 90)
print("C1  route A (measured field ratio) vs route B (closed form), term by term")
print("=" * 90)
agree = dis = 0
for r, MAX in ((1, 18), (2, 18), (3, 20), (4, 16)):
    N = t + 2 * r
    pts = make_points(r, N, 555 + r, 2)
    a = d = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            cls = classes(beta)
            if len(cls) < t:
                continue
            C = centre(beta, cls)
            if C is None:
                continue
            qB, _ = q_route_B(beta, cls, C, N)
            for al in pts:
                qA = q_route_A(beta, cls, al, C, N)
                if qA == "SUMFAIL":
                    print("   LAPLACE SUM FAILED at lam=%s" % list(l))
                    raise SystemExit(1)
                if qA == qB:
                    a += 1
                else:
                    d += 1
                    if d <= 3:
                        print("   MISMATCH r=%d lam=%s  A=%s  B=%s" % (r, list(l), qA, qB))
    print("   r=%d |lam|<=%d :  %6d agree,  %d disagree" % (r, MAX, a, d))
    agree += a
    dis += d
print("")
print("   TOTAL: %d agree, %d disagree  ->  %s" % (agree, dis, "PASS" if dis == 0 else "FAIL"))
if dis:
    raise SystemExit(1)


# ============================================================================== the wide run
def wide(r, MAX, seed):
    N = t + 2 * r
    EMAX = MAX + N
    pts = make_points(r, N, seed, 3)
    POW = [[[a ** e for e in range(EMAX + 1)] for a in al] for al in pts]

    def vanishes(beta):
        for p in range(3):
            if matrix(F, [[POW[p][b][beta[a]] for b in range(N)]
                          for a in range(N)]).det() != 0:
                return False
        return True

    nsh = zeros = nonconc_zero = nonconst = 0
    tab = {}
    rows = []
    t0 = time.time()
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            lam = list(l)
            beta = beta_of(lam, N)
            cls = classes(beta)
            if len(cls) < t:
                continue
            nsh += 1
            v = vanishes(beta)
            zeros += 1 if v else 0
            C = centre(beta, cls)
            if C is None:
                if v:
                    nonconc_zero += 1
                    if nonconc_zero <= 3:
                        print("   C3 VIOLATION: non-concentric zero lam=%s beta=%s" % (lam, beta))
                continue
            q, const = q_route_B(beta, cls, C, N)
            if not const:
                nonconst += 1
            tab[(q, v)] = tab.get((q, v), 0) + 1
            rows.append((lam, beta, cls, C, q, v))
    return dict(r=r, N=N, MAX=MAX, nsh=nsh, zeros=zeros, tab=tab, rows=rows,
                nonconc_zero=nonconc_zero, nonconst=nonconst, el=time.time() - t0)


print("")
print("=" * 90)
print("THE CRITERION ON A BIG SAMPLE:  Phi = 0  <=>  concentric and q = -1")
print("=" * 90)
print("")
print("   r   N |lam|<=   shapes  concentric  ZEROS   q=-1&0   q=-1&!0   q=+1&0   q=+1&!0   time")
print("  " + "-" * 88)
RUNS = []
for r, MAX, seed in ((3, 50, 4245), (1, 40, 7001), (2, 36, 7002), (4, 24, 7004)):
    S = wide(r, MAX, seed)
    RUNS.append(S)
    T = S['tab']
    print("  %2d %3d %6d %9d %11d %6d %8d %9d %8d %9d %6.1fs"
          % (S['r'], S['N'], S['MAX'], S['nsh'], sum(T.values()), S['zeros'],
             T.get((-1, True), 0), T.get((-1, False), 0),
             T.get((1, True), 0), T.get((1, False), 0), S['el']))

print("")
print("  C3  non-concentric zeros (must be 0) :", [S['nonconc_zero'] for S in RUNS])
print("  C2  q not constant over terms        :", [S['nonconst'] for S in RUNS])
print("  C4  the two failing cells are  'q=-1 & does not vanish'  and  'q=+1 & vanishes'.")
tot_bad = sum(S['tab'].get((-1, False), 0) + S['tab'].get((1, True), 0) for S in RUNS)
tot_zero = sum(S['zeros'] for S in RUNS)
print("      total exceptions over the four r : %d      total zeros: %d" % (tot_bad, tot_zero))

# ---- C5: the control that CAN fail -- the retracted reading, q' = q*(-1)^r
print("")
print("  C5  the SAME tally scored against the retracted reading q' = q*(-1)^r, which is what")
print("      'the criterion alternates with the parity of r' asserts.  It must do badly:")
for S in RUNS:
    sgn = (-1) ** S['r']
    bad = sum(1 for lam, beta, cls, C, q, v in S['rows'] if ((q * sgn) == -1) != v)
    print("      r=%d  exceptions under the alternating reading: %5d of %d concentric shapes"
          % (S['r'], bad, len(S['rows'])))

# ---- what q depends on
print("")
print("  q against coarse invariants of beta, on the r = 3 run (looking for a closed form):")
S = RUNS[0]
cells = {}
for lam, beta, cls, C, q, v in S['rows']:
    exc = sorted(k for k in cls if len(cls[k]) >= 2)
    m = sum(len(cls[k]) for k in exc)
    k = len(exc)
    selfp = sum(1 for kk in exc if (C - kk) % t == kk % t)
    singles = [cls[kk][0] for kk in cls if len(cls[kk]) == 1]
    lo = min(beta[i] for kk in exc for i in cls[kk])
    hi = max(beta[i] for kk in exc for i in cls[kk])
    below = sum(1 for i in singles if lo < beta[i] < hi)
    key = (m, k, selfp, below, tuple(sorted(len(cls[kk]) for kk in cls)))
    cells.setdefault(key, set()).add(q)
mixed = [c for c, s in cells.items() if len(s) > 1]
print("    cells (m,k,selfpaired,below,sizes) : %d      of them carrying BOTH signs of q : %d"
      % (len(cells), len(mixed)))
for c in sorted(cells)[:12]:
    print("      m=%d k=%d selfp=%d below=%d sizes=%s  ->  q in %s"
          % (c[0], c[1], c[2], c[3], c[4], sorted(cells[c])))

print("")
print("  five zeros of the wide r = 3 run, by hand:")
n = 0
for lam, beta, cls, C, q, v in S['rows']:
    if v and n < 5:
        print("    lam=%-28s beta=%s  C=%d  q=%+d" % (lam, beta, C, q))
        n += 1

print("")
print("DONE")
