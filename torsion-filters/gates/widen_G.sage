# -*- coding: utf-8 -*-
# WIDENING THE ONE SIZE THAT MATTERS: is |G| <= 2 real, or an artefact of the range?
#
# stratum_sign.sage found H4 exact on 12613 shapes -- [Phi_t]_top = 0 iff the degree-maximiser set G
# is closed under T -> C - T with no fixed point -- and its own decoy D2 showed the SIGN condition
# is redundant: closed and fixed-point-free already fires exactly 209 times with 0 false positives.
# But every one of those tables was computed in a range where |G| never exceeded 2, and with |G| <= 2
# the biconditional is nearly a tautology: closed and fixed-point-free just means "there are two
# maximisers and the reflection swaps them".  A criterion that is only tested where it cannot be
# stressed is not tested.  So:
#
#   W1  push |lambda| far past the earlier range on the SAME configurations, and add r = 4 and
#       t = 10, 12, and report the full distribution of |G|.  If |G| > 2 never occurs, that is a
#       statement about the geometry of transversals worth proving.  If it does occur, H4 has to
#       survive on those shapes or it dies.
#   W2  re-run H4 and the two decoys on the widened range, reported separately on the shapes with
#       |G| > 2 if there are any -- that is the only place the sign can earn its keep.
#   W3  the same for the proved corollary: how often is the maximiser unique, far out?
#
# CONTROLS:
#   K1  acceptance, fatal: on the shapes shared with the earlier range the counts must reproduce
#       (t=4 r=1 |lam|<=30 must again give 141 criterion shapes and 141 H4 hits).
#   K2  forced: criterion-holding shapes must always be closed and fixed-point-free.
#   K3  non-vacuity: the widened range must actually contain shapes the old one did not.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

WIDE = [(4, 1, 44), (4, 2, 34), (4, 3, 26), (4, 4, 20),
        (6, 2, 26), (6, 3, 20), (8, 2, 22), (8, 3, 16), (10, 2, 20), (12, 2, 18)]
OLD = {(4, 1): 30, (4, 2): 24, (4, 3): 18, (6, 2): 18, (6, 3): 14, (8, 2): 16}


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
    return {k: v for k, v in D.items() if v != 0}


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
    best = None
    G = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        d = sum(T[:r]) - sum(T[r:])
        if best is None or d > best:
            best, G = d, []
        if d == best:
            w = perm_sign([beta[i] % t for i in P])
            if sum(P) % 2:
                w = -w
            G.append((w, T))
    W = {}
    for w, T in G:
        W[T] = w
    Ts = list(W)

    def rho(T):
        return tuple(sorted((C - x for x in T), reverse=True))

    closed = all(rho(T) in W for T in Ts)
    fpf = all(rho(T) != T for T in Ts)
    signs = closed and fpf and all(W[rho(T)] == -W[T] for T in Ts)
    tot = {}
    for T in Ts:
        for k, v in topdeg_dict(list(T), r).items():
            tot[k] = tot.get(k, 0) + W[T] * v
    top0 = not any(v != 0 for v in tot.values())
    return dict(crit=conc and cond_ii, cond_ii=cond_ii, nG=len(Ts),
                closed=closed, fpf=fpf, signs=signs, top0=top0)


print("=" * 106)
print("W1/W2/W3  the widened range")
print("=" * 106)
print("")
print("     t   r  |lam|<=  shapes   crit | max |G|   |G| distribution        | H4 hits  really 0"
      "  H4 bad | D2 bad  uniq |G|=1")
print("  " + "-" * 102)

BIG = []
K1 = {}
tot = dict(sh=0, cr=0, h=0, z=0, bad=0, d2=0, u=0)
for t, r, MAX in WIDE:
    N = t + 2 * r
    nsh = ncr = h = z = bad = d2bad = uniq = 0
    dist = {}
    k1c = k1h = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            nsh += 1
            dist[a['nG']] = dist.get(a['nG'], 0) + 1
            if a['crit']:
                ncr += 1
                if not (a['closed'] and a['fpf']):
                    bad += 1
            if a['nG'] == 1:
                uniq += 1
            pred = a['closed'] and a['fpf'] and a['signs']
            d2 = a['closed'] and a['fpf']
            if pred:
                h += 1
            if a['top0']:
                z += 1
            if pred != a['top0']:
                bad += 1
            if d2 != a['top0']:
                d2bad += 1
            if a['nG'] > 2:
                BIG.append((t, r, list(l), a))
            if (t, r) in OLD and size <= OLD[(t, r)]:
                k1c += 1 if a['crit'] else 0
                k1h += 1 if pred else 0
    K1[(t, r)] = (k1c, k1h)
    tot['sh'] += nsh
    tot['cr'] += ncr
    tot['h'] += h
    tot['z'] += z
    tot['bad'] += bad
    tot['d2'] += d2bad
    tot['u'] += uniq
    dd = ", ".join("%d:%d" % kv for kv in sorted(dist.items()))
    print("  %4d %3d %8d %7d %6d | %7d   %-22s | %7d %9d %7d | %6d %8d"
          % (t, r, MAX, nsh, ncr, max(dist), dd[:22], h, z, bad, d2bad, uniq))

print("")
print("  totals: %d shapes, %d criterion.  H4 fires %d, really zero %d, H4 wrong %d."
      % (tot['sh'], tot['cr'], tot['h'], tot['z'], tot['bad']))
print("  D2 (closure + fixed-point-free, NO sign) wrong %d times." % tot['d2'])
print("  unique maximiser (the proved corollary applies) on %d of %d shapes (%.1f%%)."
      % (tot['u'], tot['sh'], 100.0 * tot['u'] / tot['sh']))
print("")
print("  K1 reproduction of the earlier range (criterion count, H4 hits):")
for k in sorted(K1):
    if k in OLD:
        print("     t=%d r=%d |lam|<=%d -> %s" % (k[0], k[1], OLD[k], str(K1[k])))
print("     stratum_sign_OUT.txt gave, on those same six rows, criterion 141/30/6/7/2/3 and")
print("     H4 hits 141/44/8/10/3/3.  Those are the pairs that must be reproduced above.")
print("")
print("  W1  shapes with |G| > 2: %d" % len(BIG))
if BIG:
    print("")
    print("     t   r  lambda                      |G|  closed  fpf  signs  [Phi]top=0")
    print("  " + "-" * 102)
    nsg = 0
    for t, r, lam, a in BIG[:25]:
        print("  %4d %3d  %-27s %3d %7s %5s %6s %10s"
              % (t, r, str(lam), a['nG'], a['closed'], a['fpf'], a['signs'], a['top0']))
    for t, r, lam, a in BIG:
        if (a['closed'] and a['fpf']) != (a['closed'] and a['fpf'] and a['signs']):
            nsg += 1
    print("")
    print("     of the %d shapes with |G| > 2, the sign condition changes the verdict on %d."
          % (len(BIG), nsg))
else:
    print("  |G| never exceeds 2 in this whole range.  Then 'closed and fixed-point-free' just")
    print("  says |G| = 2 and the reflection swaps the two, and the sign condition cannot be")
    print("  tested at all -- which is exactly why D2 scores like H4.")
print("")
print("DONE")
