# -*- coding: utf-8 -*-
# A closed form for the reflection sign q, and the criterion at a second t.
#
# WHERE THIS STANDS.  parity_terms.sage measured the ratio q = term(C-T)/term(T) of the Laplace
# expansion with no sign convention; parity_closed.sage gave it a closed form checked term by term
# against the measurement (876 agreements, 0 disagreements) and then established, on 516763 shapes
# at t = 4 and r = 1,2,3,4,
#
#       Phi_t(lambda; z) == 0    <=>    beta concentric  and  q(beta) = -1,
#
# with 821 zeros, 7311 concentric shapes and NOT ONE exception in either direction.  There is no
# r-parity: the retracted alternation was a (-1)^r picked up by mixing two row conventions, and the
# control shows its signature exactly -- it is right on every even r and wrong on every odd r.
#
# TWO THINGS ARE STILL MISSING and this script goes after both.
#
# (1) q is defined by a sum over transversals.  The r=3 table hinted it depends only on
#     s = #{excess classes k with C - k = k mod t}, through q = (-1)^{binom(s,2)}.  That is fitted
#     here over EVERY concentric shape, not only the vanishing ones -- fitting on the zeros alone
#     proves nothing, since q = -1 is what vanishing means.  Rival candidates are scored beside it.
#
# (2) everything so far is t = 4.  The criterion is re-run whole at t = 6, determinants included,
#     which is the first test that it is not an accident of one alphabet.
#
# CONTROLS able to fail:
#   K1  a candidate must agree on ALL concentric shapes; the tally is printed for each.
#   K2  cells of the invariant carrying both signs of q are printed; any is fatal to that fit.
#   K3  at t = 6: non-concentric zeros must be 0, and the two failing cells of the biconditional
#       ('q=-1 and does not vanish', 'q=+1 and vanishes') must be empty.
#   K4  a deliberately wrong grading, q'' = (-1)^s, is scored beside the candidate.  If it also
#       fits, the sample does not discriminate and nothing is claimed.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools
import time

PRIME = 2013265921            # 15*2^27 + 1 : has both 4th and 6th roots of unity
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


def word_sign(seq):
    sg = 1
    s = list(seq)
    for a in range(len(s)):
        for b in range(a + 1, len(s)):
            if s[a] > s[b]:
                sg = -sg
    return sg


def q_of(beta, cls, C, N, t):
    pos = dict((b, i) for i, b in enumerate(beta))
    iota = dict((i, pos[C - b]) for i, b in enumerate(beta) if (C - b) in pos)
    vals = set()
    for pick in itertools.product(*[cls[k] for k in sorted(cls)]):
        pick = tuple(sorted(pick))
        T = tuple(sorted(set(range(N)) - set(pick)))
        if any(i not in iota for i in T):
            return None
        T2 = tuple(sorted(iota[i] for i in T))
        pick2 = tuple(sorted(set(range(N)) - set(T2)))
        vals.add((-1) ** (sum(pick) + sum(pick2))
                 * word_sign([beta[i] % t for i in pick])
                 * word_sign([beta[i] % t for i in pick2]))
        if len(vals) > 1:
            return 0
    return vals.pop()


def invariants(beta, cls, C, t):
    exc = sorted(k for k in cls if len(cls[k]) >= 2)
    s = sum(1 for k in exc if (C - k) % t == k % t)
    m = sum(len(cls[k]) for k in exc)
    return dict(s=s, k=len(exc), m=m, pairs=(len(exc) - s) // 2,
                sizes=tuple(sorted(len(cls[k]) for k in cls)))


CAND = {
    "(-1)^binom(s,2)": lambda d: (-1) ** (d['s'] * (d['s'] - 1) // 2),
    "(-1)^s            [K4 decoy]": lambda d: (-1) ** d['s'],
    "(-1)^(s//2)": lambda d: (-1) ** (d['s'] // 2),
    "(-1)^k": lambda d: (-1) ** d['k'],
    "(-1)^binom(m,2)": lambda d: (-1) ** (d['m'] * (d['m'] - 1) // 2),
    "(-1)^pairs": lambda d: (-1) ** d['pairs'],
    "(-1)^(binom(s,2)+pairs)": lambda d: (-1) ** (d['s'] * (d['s'] - 1) // 2 + d['pairs']),
}

# ==================================================== (1) the fit, t = 4, no determinants needed
print("=" * 92)
print("THE CLOSED FORM FOR q,  fitted over EVERY concentric shape (not only the zeros),  t = 4")
print("=" * 92)
print("")
t = 4
DATA = []
for r, MAX in ((1, 44), (2, 40), (3, 34), (4, 26), (5, 22)):
    N = t + 2 * r
    n = 0
    t0 = time.time()
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            cls = classes(beta, t)
            if len(cls) < t:
                continue
            C = centre(beta, cls)
            if C is None:
                continue
            q = q_of(beta, cls, C, N, t)
            if q is None or q == 0:
                continue
            d = invariants(beta, cls, C, t)
            d['q'] = q
            d['r'] = r
            DATA.append(d)
            n += 1
    print("  r=%d  N=%2d  |lam|<=%2d :  %6d concentric shapes   [%.1fs]" % (r, N, MAX, n, time.time() - t0))

print("")
print("  %-32s %s" % ("candidate", "agrees with the measured q"))
print("  " + "-" * 66)
for name, f in CAND.items():
    ok = sum(1 for d in DATA if f(d) == d['q'])
    print("  %-32s %7d / %-7d %s" % (name, ok, len(DATA),
                                     "  <-- EXACT" if ok == len(DATA) else ""))

print("")
print("  K2  cells of s carrying both signs of q (any is fatal to the fit):")
cell = {}
for d in DATA:
    cell.setdefault(d['s'], set()).add(d['q'])
for s in sorted(cell):
    n = sum(1 for d in DATA if d['s'] == s)
    print("      s=%d  ->  q in %-9s  (%d shapes)" % (s, sorted(cell[s]), n))
print("      mixed cells: %d" % sum(1 for s in cell if len(cell[s]) > 1))

print("")
print("  q by r, to show the closed form does not depend on r:")
for r in sorted(set(d['r'] for d in DATA)):
    sub = [d for d in DATA if d['r'] == r]
    ok = sum(1 for d in sub if (-1) ** (d['s'] * (d['s'] - 1) // 2) == d['q'])
    print("      r=%d  (-1)^binom(s,2) agrees %5d / %-5d" % (r, ok, len(sub)))


# ============================================================ (2) the whole criterion at t = 6
print("")
print("=" * 92)
print("THE CRITERION AT t = 6,  determinants included -- the first alphabet that is not t = 4")
print("=" * 92)
print("")
print("   t   r   N |lam|<=  shapes  concentric  ZEROS   q=-1&0  q=-1&!0  q=+1&0  q=+1&!0  fit")
print("  " + "-" * 88)

for t6, r, MAX, seed in ((6, 1, 26, 11), (6, 2, 22, 12), (6, 3, 18, 13)):
    N = t6 + 2 * r
    EMAX = MAX + N
    RT = roots_of(t6)
    set_random_seed(seed)
    pts = []
    while len(pts) < 3:
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        al = RT + [y for x in zz for y in (x, 1 / x)]
        if len(set(al)) == N:
            pts.append(al)
    POW = [[[a ** e for e in range(EMAX + 1)] for a in al] for al in pts]

    nsh = zeros = ncon = nonconc = nfit = 0
    tab = {}
    t0 = time.time()
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            cls = classes(beta, t6)
            if len(cls) < t6:
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
            C = centre(beta, cls)
            if C is None:
                if v:
                    nonconc += 1
                continue
            q = q_of(beta, cls, C, N, t6)
            if q is None or q == 0:
                continue
            ncon += 1
            tab[(q, v)] = tab.get((q, v), 0) + 1
            d = invariants(beta, cls, C, t6)
            if (-1) ** (d['s'] * (d['s'] - 1) // 2) == q:
                nfit += 1
    print("  %2d %3d %3d %6d %8d %11d %6d %7d %8d %7d %8d  %d/%d  [%.0fs]"
          % (t6, r, N, MAX, nsh, ncon, zeros, tab.get((-1, True), 0), tab.get((-1, False), 0),
             tab.get((1, True), 0), tab.get((1, False), 0), nfit, ncon, time.time() - t0))
    if nonconc:
        print("      K3 VIOLATION: %d non-concentric zeros" % nonconc)

print("")
print("  the last two columns before 'fit' are the failing cells: 'q=-1 and does not vanish'")
print("  and 'q=+1 and vanishes'.  Both must be 0 for the biconditional to hold at t = 6.")
print("")
print("DONE")
