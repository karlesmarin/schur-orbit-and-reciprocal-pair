# -*- coding: utf-8 -*-
# The term signs of the Laplace expansion, computed instead of derived.
#
# WHY THIS SCRIPT EXISTS.  parity_r3_wide.sage widened r = 3 from 6 zeros to 29 and its control C1
# fired: 125 concentric shapes have eps = +1 and DO NOT VANISH.  Under the reading of
# parity_alternates.sage -- r odd cancels pairwise iff eps = +1 -- each of those 125 would have to
# vanish.  So the quantity called eps there is not the ratio of the two term signs.  Two candidate
# omissions, both about conventions:
#
#   (i)  the t x t root-of-unity minor is not a sign, it is +- V (Vandermonde), and which sign
#        depends on the ORDER in which the picked residues appear -- that varies with the pick;
#   (ii) reflecting beta -> C - beta reverses the row order of the 2r x 2r block, and the reversal
#        sign (-1)^{r} exactly cancels the (-1)^{r} that column inversion produces, so with rows
#        SORTED -- which is the convention the Laplace sign assumes -- A(T2) = +A(T) for every r.
#
# If (ii) is right the r-parity alternation is an artefact of mixing two conventions.  Rather than
# argue conventions, this computes the terms.  For each shape and each point it forms
#
#     term(T)  =  (-1)^{sum(pick) + t(t-1)/2} * det M[pick, roots] * det M[T, zs],
#
# checks sum_T term(T) == det M  (fatal control: if that fails nothing below means anything), and
# then asks the reflection directly:  what IS term(T2)/term(T)?  No sign convention is used
# anywhere; the ratio is a field element that is measured.
#
# ACCEPTANCE, fatal:
#   A1  the Laplace terms must sum to the determinant, on every shape tested (the control is that
#       it is checked on all of them, not on one).
#   A2  r=3, |lambda| <= 20 reproduces 2387 shapes / 6 zeros.
#
# CONTROLS able to fail:
#   C1  ratio(T) must be +-1 and must not depend on the point.  If it does, the pairing is not a
#       sign at all and the whole picture dies.
#   C2  "every pair cancels" must imply "vanishes", 100%, with no exception.  This is now a
#       theorem about the code as much as about the algebra.
#   C3  shapes that vanish while some pair does NOT cancel must be reported by hand.
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


def laplace_terms(beta, cls, alph, N):
    """{T (sorted tuple of the 2r non-picked rows) : term value}, plus det M."""
    M = matrix(F, [[alph[b] ** beta[a] for b in range(N)] for a in range(N)])
    out = {}
    for pick in itertools.product(*[cls[k] for k in sorted(cls)]):
        pick = tuple(sorted(pick))
        T = tuple(sorted(set(range(N)) - set(pick)))
        s = (-1) ** (sum(pick)) * SGN_T
        out[T] = s * M[list(pick), list(range(t))].det() * M[list(T), list(range(t, N))].det()
    return out, M.det()


def reflect(beta, C):
    pos = dict((b, i) for i, b in enumerate(beta))
    iota = {}
    for i, b in enumerate(beta):
        j = pos.get(C - b)
        if j is not None:
            iota[i] = j
    return iota


# ------------------------------------------------------------------------------------- the sweep
def run(r, MAX, seed):
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

    nsh = 0
    conc = []           # (lam, beta, cls, C, vanishes)
    zeros = 0
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
            if v:
                zeros += 1
            C = centre(beta, cls)
            if C is not None:
                conc.append((lam, beta, cls, C, v))

    # --- the term analysis, on the concentric shapes only (all zeros are concentric)
    nsum_bad = nratio_bad = npoint_bad = 0
    rows = []
    for lam, beta, cls, C, v in conc:
        iota = reflect(beta, C)
        cancels = True
        ratios = set()
        for p in range(3):
            tm, dd = laplace_terms(beta, cls, pts[p], N)
            if sum(tm.values()) != dd:
                nsum_bad += 1
            per = {}
            for T, val in tm.items():
                T2 = tuple(sorted(iota[i] for i in T))
                if T2 not in tm:
                    cancels = False
                    continue
                if tm[T2] != -val:
                    cancels = False
                if val != 0:
                    q = tm[T2] / val
                    if q not in (F(1), F(-1)):
                        nratio_bad += 1
                    per[T] = 1 if q == F(1) else (-1 if q == F(-1) else 0)
            if p == 0:
                ref = per
            elif per != ref:
                npoint_bad += 1
            ratios |= set(per.values())
        rows.append(dict(lam=lam, beta=beta, cls=cls, C=C, v=v, cancels=cancels,
                         ratios=tuple(sorted(ratios))))
    return dict(r=r, N=N, MAX=MAX, nsh=nsh, zeros=zeros, conc=rows, el=time.time() - t0,
                nsum_bad=nsum_bad, nratio_bad=nratio_bad, npoint_bad=npoint_bad)


# ------------------------------------------------------------------------------------ acceptance
A2 = run(3, 20, 4245)
ok = (A2['nsh'] == 2387 and A2['zeros'] == 6 and A2['nsum_bad'] == 0)
print("A2  r=3 |lam|<=20 : shapes=%d(want 2387)  zeros=%d(want 6)  "
      "Laplace-sum failures=%d(want 0)  ->  %s"
      % (A2['nsh'], A2['zeros'], A2['nsum_bad'], "PASS" if ok else "FAIL"))
if not ok:
    raise SystemExit(1)

print("")
print("=" * 90)
print("WHAT THE REFLECTION ACTUALLY DOES TO A TERM,  t = 4")
print("=" * 90)
print("")
print("   r   N  |lam|<=  shapes  concentric  zeros   all pairs cancel   of those, vanish")
print("  " + "-" * 86)

RUNS = []
for r, MAX, seed in ((1, 30, 7001), (2, 26, 7002), (3, 34, 4245), (4, 18, 7004)):
    S = run(r, MAX, seed)
    RUNS.append(S)
    canc = [d for d in S['conc'] if d['cancels']]
    print("  %2d %3d %7d %8d %11d %6d %18d %17d"
          % (S['r'], S['N'], S['MAX'], S['nsh'], len(S['conc']), S['zeros'],
             len(canc), sum(1 for d in canc if d['v'])))

print("")
print("  CONTROLS (all must be 0):")
for S in RUNS:
    print("    r=%d   Laplace sum != det : %d    ratio not +-1 : %d    ratio depends on point : %d"
          % (S['r'], S['nsum_bad'], S['nratio_bad'], S['npoint_bad']))

print("")
print("  C2  'all pairs cancel' => 'vanishes'  (any exception kills the picture):")
for S in RUNS:
    bad = [d for d in S['conc'] if d['cancels'] and not d['v']]
    print("    r=%d   exceptions: %d" % (S['r'], len(bad)))
    for d in bad[:3]:
        print("        lam=%s beta=%s C=%d" % (d['lam'], d['beta'], d['C']))

print("")
print("  C3  zeros the reflection does NOT explain, per r:")
for S in RUNS:
    un = [d for d in S['conc'] if d['v'] and not d['cancels']]
    print("    r=%d   %d of %d zeros unexplained" % (S['r'], len(un), S['zeros']))
    for d in un[:3]:
        print("        lam=%-24s beta=%s  C=%d  ratios seen=%s"
              % (d['lam'], d['beta'], d['C'], d['ratios']))

print("")
print("  the ratio multiset per r, over ALL concentric shapes -- if the pairing were governed by")
print("  the parity of r, the r odd rows and the r even rows would look different:")
for S in RUNS:
    tab = {}
    for d in S['conc']:
        tab[d['ratios']] = tab.get(d['ratios'], 0) + 1
    print("    r=%d  %s" % (S['r'], sorted(tab.items(), key=lambda kv: -kv[1])))

print("")
print("DONE")
