# -*- coding: utf-8 -*-
# THE TOP STRATUM HAS ITS OWN (i)+(ii): closure is the concentricity, and the missing half is a SIGN.
#
# topstratum.sage measured that closure of the degree-maximiser set G under T -> C - T is NECESSARY
# for [Phi_t]_top = 0 (0 shapes in 12613 vanish without it) and nowhere near sufficient (4438 closed
# shapes do not vanish).  That is exactly the shape of the main theorem, where concentricity alone
# fails in 2212 of 2436 shapes and the second condition does the rest.  So the stratum should have
# its own second condition, and it should be a sign.
#
# ------------------------------------------------------------------------------------------------
# LEMMA (proved, three lines, machine-checked below).  Write P(T) = a_H(z) a_L(1/z) for the
# top-degree part of A(T), with H the top r values of T and L the bottom r.  Then
#       P(C - T)  =  P(T)      for every T and every C.
# Proof.  The top half of C - T is C - L listed decreasingly, i.e. L listed increasingly, so
#       a_{C-L}(z) = (prod_j z_j^C) * (-1)^{r(r-1)/2} * a_L(1/z),
# and symmetrically a_{C-H}(1/z) = (prod_j z_j^{-C}) * (-1)^{r(r-1)/2} * a_H(z).  The two powers of
# prod z_j^C cancel and (-1)^{r(r-1)} = +1 because r(r-1) is even.  QED
#
# CONSEQUENCE.  [Phi_t]_top = sum_{T in G} w(T) P(T) with P constant on the orbits of T -> C - T.
# If T -> C - T has no fixed point on G and P separates the orbits, the sum vanishes if and only if
#       G is closed,  and  w(C - T) = -w(T)  for every T in G.
# That is a criterion with no polynomial in it: a closure test and |G|/2 sign comparisons.
# ------------------------------------------------------------------------------------------------
#
# TESTED, each able to fail:
#   H4  the biconditional:  [Phi]_top = 0  <=>  G closed, fixed-point-free, and w(C-T) = -w(T).
#       Both directions counted separately; either can fail.
#   D1  DECOY, closure alone -- known to over-fire, printed so the gain is visible.
#   D2  DECOY, closure plus fixed-point-freeness but WITHOUT the sign.  If D2 scores as well as H4
#       the sign is not doing the work.
#   L1  the lemma P(C-T) = P(T), checked on every maximiser of every shape.
#   S1  does P separate the orbits?  #distinct P against |G|/2.  If it does not, the criterion
#       above is only sufficient and the biconditional has to survive on its own.
#
# CONTROLS:
#   K1  non-vacuity: each cell of the H4 contingency table that is meant to be non-empty must be.
#   K2  forced: on criterion-holding shapes H4's right-hand side must hold.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

CONF = [(4, 1, 30), (4, 2, 24), (4, 3, 18), (6, 2, 18), (6, 3, 14), (8, 2, 16)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def perm_sign(q):
    n = len(q)
    seen = [False] * n
    s = 1
    for i in range(n):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = q[j]
            ln += 1
        if ln % 2 == 0:
            s = -s
    return s


def topdeg_dict(T, r):
    D = {}
    n = 2 * r
    for a in itertools.permutations(range(r)):
        for b in itertools.permutations(range(r)):
            q = [0] * n
            e = [0] * r
            for i in range(r):
                q[i] = 2 * a[i]
                e[a[i]] += T[i]
            for i in range(r):
                q[r + i] = 2 * b[i] + 1
                e[b[i]] -= T[r + i]
            k = tuple(e)
            D[k] = D.get(k, 0) + perm_sign(q)
    return tuple(sorted((k, v) for k, v in D.items() if v != 0))


def analyse(beta, t, r):
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    S = sorted((beta[i] for k in E for i in cl[k]), reverse=True)
    C = S[0] + S[-1]
    conc = set(C - b for b in S) == set(S)
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    cond_ii = (len(fixed) == 2 and all(k in E for k in fixed))
    tm = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        Ps = set(P)
        tm.append((w, tuple(beta[i] for i in range(N) if i not in Ps)))
    degs = [sum(T[:r]) - sum(T[r:]) for _, T in tm]
    Dmax = max(degs)
    G = [tm[i] for i, d in enumerate(degs) if d == Dmax]
    W = dict((T, w) for w, T in G)
    Ts = list(W)
    P = dict((T, topdeg_dict(list(T), r)) for T in Ts)

    def rho(T):
        return tuple(sorted((C - x for x in T), reverse=True))

    closed = all(rho(T) in W for T in Ts)
    lemma_bad = sum(1 for T in Ts if rho(T) in P and P[rho(T)] != P[T])
    fpf = all(rho(T) != T for T in Ts)
    signs = closed and fpf and all(W[rho(T)] == -W[T] for T in Ts)
    # the actual top-degree part
    tot = {}
    for T in Ts:
        for k, v in P[T]:
            tot[k] = tot.get(k, 0) + W[T] * v
    top0 = not any(v != 0 for v in tot.values())
    nP = len(set(P.values()))
    return dict(crit=conc and cond_ii, conc=conc, cond_ii=cond_ii, nG=len(Ts), C=C,
                closed=closed, fpf=fpf, signs=signs, top0=top0,
                lemma_bad=lemma_bad, nP=nP)


print("=" * 104)
print("L1  the lemma  P(C - T) = P(T)      S1  does P separate the reflection orbits of G?")
print("=" * 104)
print("")
print("     t   r   shapes   maximisers seen   lemma violations   #G > 2   |G| vs #distinct P")
print("  " + "-" * 100)

l1bad = 0
for t, r, MAX in CONF:
    N = t + 2 * r
    nsh = nmax = lb = big = 0
    pairs = {}
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            nsh += 1
            nmax += a['nG']
            lb += a['lemma_bad']
            if a['nG'] > 2:
                big += 1
            pairs[(a['nG'], a['nP'])] = pairs.get((a['nG'], a['nP']), 0) + 1
    l1bad += lb
    top5 = sorted(pairs.items(), key=lambda kv: -kv[1])[:4]
    print("  %4d %3d %8d %17d %18d %8d   %s"
          % (t, r, nsh, nmax, lb, big,
             ", ".join("|G|=%d:#P=%d (%d)" % (k[0], k[1], v) for k, v in top5)))

print("")
if l1bad:
    print("  L1 FAILED -- the lemma is false, stop.")
    raise SystemExit(1)
print("  L1 PASS: P(C-T) = P(T) on every maximiser of every shape.")
print("  S1: read the last column.  |G| = 2k with #P = k means P separates the orbits exactly.")

# ---------------------------------------------------------------- H4, D1, D2 ---------------------
print("")
print("=" * 104)
print("H4  [Phi]_top = 0  <=>  G closed, fixed-point-free, and w(C-T) = -w(T)")
print("=" * 104)
print("")
print("     t   r   shapes | H4 says 0  really 0 | H4 wrong -> | says0 not0   not-say0 is0"
      " | D1 closed  D2 +fpf")
print("  " + "-" * 100)

tot = dict(a=0, b=0, e1=0, e2=0, d1=0, d2=0, sh=0)
K2bad = 0
for t, r, MAX in CONF:
    N = t + 2 * r
    nsh = sa = sb = e1 = e2 = d1 = d2 = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            nsh += 1
            sa += 1 if a['signs'] else 0
            sb += 1 if a['top0'] else 0
            if a['signs'] and not a['top0']:
                e1 += 1
            if a['top0'] and not a['signs']:
                e2 += 1
            d1 += 1 if a['closed'] else 0
            d2 += 1 if (a['closed'] and a['fpf']) else 0
            if a['crit'] and not a['signs']:
                K2bad += 1
    print("  %4d %3d %8d | %9d %9d | %12s %11d %13d | %8d %8d"
          % (t, r, nsh, sa, sb, "", e1, e2, d1, d2))
    tot['sh'] += nsh
    tot['a'] += sa
    tot['b'] += sb
    tot['e1'] += e1
    tot['e2'] += e2
    tot['d1'] += d1
    tot['d2'] += d2

print("")
print("  totals over %d shapes: H4 predicts vanishing %d times, it really vanishes %d times."
      % (tot['sh'], tot['a'], tot['b']))
print("  H4 exceptions: predicted 0 but nonzero %d;  vanished but not predicted %d."
      % (tot['e1'], tot['e2']))
print("  D1 decoy (closure alone) fires %d times -- %d false positives."
      % (tot['d1'], tot['d1'] - tot['b']))
print("  D2 decoy (closure + fixed-point-free, no sign) fires %d times -- %d false positives."
      % (tot['d2'], tot['d2'] - tot['b']))
print("  K2 criterion-holding shapes where H4's right-hand side fails: %d (must be 0)" % K2bad)
print("  K1 non-vacuity: H4 fires (%s) and does not fire (%s)"
      % ("yes" if tot['a'] > 0 else "NO", "yes" if tot['a'] < tot['sh'] else "NO"))
print("")
print("DONE")
