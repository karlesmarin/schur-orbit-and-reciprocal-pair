# -*- coding: utf-8 -*-
# Is "missing" at t = 2 a block, or just a bigger partition than I swept?
#
# The last run left me with 495 triples missing a permutation at |lambda| <= 34, and I was about to
# call that an obstruction of unknown kind.  One of its examples can be settled by hand, which means
# it must be: (8,20,26) is claimed missing.  Parametrise the two-class profile directly --
#     a_2 = 2u, a_1 = 2u + d_1, b_2 = 2v+1, b_1 = 2v+1+d_2,   d_3 = |4(u-v) + d_1 - d_2 - 2|
# -- and solve for u, v >= 0.  That is arithmetic, not a sweep, so it answers for ALL sizes at once.
#
# For every triple reported missing, this computes the MINIMAL |lambda| that realises each of its
# absent permutations, or says no solution exists.  Two outcomes, and they are opposite:
#
#   if every absent permutation turns out realisable at some larger size, then t = 2 has no
#   obstruction at all and the 57% was the edge of my sweep -- the question I announced is void;
#
#   if some are genuinely unsolvable, those are the real obstruction and this prints them.
#
# Authors: Carles Marin, Claude (AI assistant).

t, N = 2, 4


def min_size_two_class(d1, d2, d3):
    """Least |lambda| realising (d1,d2,d3) in the two-class profile, or None."""
    if d1 <= 0 or d2 <= 0 or d1 % 2 or d2 % 2 or d3 % 2:
        return None
    best = None
    for eps in (1, -1):
        # 4(u-v) = eps*d3 - d1 + d2 + 2
        rhs = eps * d3 - d1 + d2 + 2
        if rhs % 4:
            continue
        k = rhs // 4                      # k = u - v, any integer; u,v >= 0
        for base in range(0, 40):         # v = base, u = base + k  (or the mirror)
            u, v = base + k, base
            if u < 0 or v < 0:
                continue
            beta = sorted([2 * u + d1, 2 * u, 2 * v + 1 + d2, 2 * v + 1], reverse=True)
            if len(set(beta)) < 4:
                continue
            lam = [beta[i] - (N - 1 - i) for i in range(N)]
            if any(x < 0 for x in lam) or any(lam[i] < lam[i + 1] for i in range(N - 1)):
                continue
            sz = sum(lam)
            if best is None or sz < best:
                best = sz
            break
    return best


def min_size_size_three(d1, d2, d3):
    """The size-three profile forces d3 = d1 + d2; least |lambda|, or None."""
    if d3 != d1 + d2 or d1 % 2 or d2 % 2:
        return None
    best = None
    for res in (0, 1):
        for w in range(0, 60):
            P = [w + d1 + d2, w + d2, w]
            P = [2 * x + res for x in P]
            other = [y for y in range(0, 200) if y % 2 != res]
            for o in other[:20]:
                beta = sorted(P + [o], reverse=True)
                if len(set(beta)) < 4:
                    continue
                lam = [beta[i] - (N - 1 - i) for i in range(N)]
                if any(x < 0 for x in lam) or any(lam[i] < lam[i + 1] for i in range(N - 1)):
                    continue
                sz = sum(lam)
                if best is None or sz < best:
                    best = sz
    return best


def min_size(d):
    a = min_size_two_class(*d)
    b = min_size_size_three(*d)
    xs = [x for x in (a, b) if x is not None]
    return min(xs) if xs else None


MAX = 34


def datum(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for b in beta:
        cls.setdefault(b % t, []).append(b)
    if any(len(cls.get(i, [])) == 0 for i in range(t)):
        return None
    big = sorted([i for i in range(t) if len(cls[i]) >= 2])
    if len(big) == 2:
        A = sorted(cls[big[0]], reverse=True)[:2]
        B = sorted(cls[big[1]], reverse=True)[:2]
    else:
        P = sorted(cls[big[0]], reverse=True)
        A, B = [P[0], P[1]], [P[1], P[2]]
    return (A[0] - A[1], B[0] - B[1], abs(A[0] + A[1] - B[0] - B[1]))


att = set()
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        d = datum(list(l))
        if d is not None and 0 not in d:
            att.add(d)

print("=" * 76)
print("t = 2 : are the absent permutations blocked, or just larger than |lambda| <= %d ?" % MAX)
print("=" * 76)

realisable = unsolvable = 0
sizes = []
worst = []
for d in sorted(att):
    perms = set(tuple(x) for x in Permutations(list(d)).list())
    for q in sorted(perms - att):
        m = min_size(q)
        if m is None:
            unsolvable += 1
            if len(worst) < 6:
                worst.append((d, q))
        else:
            realisable += 1
            sizes.append(m)

print("")
print("  absent permutations examined      : %d" % (realisable + unsolvable))
print("  realisable at a larger |lambda|   : %d" % realisable)
print("  with no solution at any size      : %d" % unsolvable)
if sizes:
    sizes.sort()
    print("  their minimal |lambda|: min %d, median %d, max %d  (my sweep stopped at %d)"
          % (sizes[0], sizes[len(sizes) // 2], sizes[-1], MAX))
for d, q in worst:
    print("    genuinely absent: %s is attained, %s has no solution" % (str(d), str(q)))

print("")
print("  CONTROL: the same routine on triples I DID attain must find them all, at a size <= %d." % MAX)
ok = sum(1 for d in att if (min_size(d) or 10 ** 9) <= MAX)
print("           %d of %d attained triples reproduced by the parametrisation" % (ok, len(att)))
print("")
print("DONE")
