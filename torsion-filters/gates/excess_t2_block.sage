# -*- coding: utf-8 -*-
# t = 2: what blocks a permutation, when the arithmetic mod t allows it?
#
# At t=2 every d_i is even, so the mod-t obstruction that settles t=3 and t=4 is absent, and yet
# 57% of attained triples still miss a permutation up to |lambda| <= 38.  Something else blocks it.
#
# The candidate, derived rather than guessed.  Write d_i = 2 e_i.  In the TWO-CLASS profile put
# a_2 = 2u, a_1 = 2u + d_1, b_2 = 2v+1, b_1 = 2v+1+d_2; then
#     d_3 = |4(u-v) + d_1 - d_2 - 2|,   so   e_3 = e_1 + e_2 + 1  (mod 2),
# i.e. e_1 + e_2 + e_3 is ODD, a condition symmetric in the three.  In the SIZE-THREE profile
# d_3 = d_1 + d_2 exactly, so e_1 + e_2 + e_3 = 2(e_1+e_2) is EVEN.  Two disjoint families by the
# parity of one symmetric quantity -- and the size-three one is cut out by a relation, d_3 = d_1+d_2,
# that is NOT symmetric.  So its permutations would need even parity, which only size-three can
# supply, and they are not of that form.
#
# Predictions, each falsifiable here:
#   P1  parity of e_1+e_2+e_3 separates the two profiles exactly, 0 exceptions;
#   P2  the triples missing a permutation are exactly the size-three ones;
#   P3  among two-class triples, every permutation is attained.
#
# P2 is the one I expect to fail, because 252 of 420 is far more than a one-parameter family should
# give, and if it fails the counts below say by how much and in which direction.
#
# Authors: Carles Marin, Claude (AI assistant).

t, N = 2, 4


def datum(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for b in beta:
        cls.setdefault(b % t, []).append(b)
    prof = [len(cls.get(i, [])) for i in range(t)]
    if 0 in prof:
        return None
    big = sorted([i for i in range(t) if prof[i] >= 2])
    if len(big) == 2:
        A = sorted(cls[big[0]], reverse=True)[:2]
        B = sorted(cls[big[1]], reverse=True)[:2]
        kind = "two-class"
    else:
        P = sorted(cls[big[0]], reverse=True)
        A, B = [P[0], P[1]], [P[1], P[2]]
        kind = "size-three"
    return (A[0] - A[1], B[0] - B[1], abs(A[0] + A[1] - B[0] - B[1])), kind


MAX = 34
kinds = {}
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        r = datum(list(l))
        if r is None or 0 in r[0]:
            continue
        kinds.setdefault(r[0], set()).add(r[1])

att = set(kinds)
print("=" * 74)
print("t = 2, |lambda| <= %d : %d attained triples" % (MAX, len(att)))
print("=" * 74)

# P1 -- parity of e1+e2+e3 against the profile
bad1 = 0
for d, ks in kinds.items():
    par = sum(x // 2 for x in d) % 2
    want = set(["two-class"]) if par == 1 else set(["size-three"])
    if ks != want:
        bad1 += 1
print("")
print("  P1  parity of e1+e2+e3 separates the profiles : %d exceptions in %d triples"
      % (bad1, len(att)))
print("      (odd -> two-class, even -> size-three)")

# P2, P3 -- who misses a permutation
tab = {}
for d, ks in kinds.items():
    perms = set(tuple(x) for x in Permutations(list(d)).list())
    full = perms <= att
    key = (tuple(sorted(ks)), full)
    tab.setdefault(key, []).append(d)

print("")
print("  %-14s %-14s %7s" % ("profile", "orbit", "count"))
print("  " + "-" * 40)
for (ks, full), ds in sorted(tab.items()):
    print("  %-14s %-14s %7d" % ("/".join(ks), "complete" if full else "MISSING", len(ds)))
    if not full:
        d = sorted(ds)[len(ds) // 2]
        perms = set(tuple(x) for x in Permutations(list(d)).list())
        print("      example %s ; not attained: %s"
              % (str(d), sorted(perms - att)[:2]))

# and the direct question: of the missing, how many are missing ONLY because a permutation
# has even parity while the family that could supply it is the size-three one
odd_only = 0
for d, ks in kinds.items():
    perms = set(tuple(x) for x in Permutations(list(d)).list())
    if perms <= att:
        continue
    miss = perms - att
    if all(sum(x // 2 for x in q) % 2 == sum(x // 2 for x in d) % 2 for q in miss):
        odd_only += 1
print("")
print("  every missing permutation has the SAME parity as its triple: %d of %d missing"
      % (odd_only, sum(len(v) for (k, f), v in tab.items() if not f)))
print("  -- so parity is NOT what blocks them, and the block is elsewhere.")
print("")
print("DONE")
