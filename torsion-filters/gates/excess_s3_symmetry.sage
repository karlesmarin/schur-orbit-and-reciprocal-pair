# -*- coding: utf-8 -*-
# The formula is symmetric in d1, d2, d3.  Is anything else?
#
# Theorem 3.1 is symmetric in the three arguments, but the construction that produces them is not:
# d1 and d2 are INTERNAL gaps of the two distinguished quotient components, while d3 is a
# DIFFERENCE OF SIZES between them (eq. after Prop. quotient, d3 = t(|lam^(rA)|-|lam^(rB)|) +
# 2(rA-rB) in the two-class profile).  Two species of quantity, and the formula treats them as one.
#
# If they really are one species, the S_3 that permutes them should be realised by something done to
# the partition.  That is a falsifiable prediction, and this measures it:
#
#   attained      how many ordered triples (d1,d2,d3) occur
#   full orbits   of the attained triples, how many have ALL their distinct permutations attained
#   |lambda|      whether a permutation can be realised at constant size -- if not, whatever
#                 realises it is not a size-preserving operation on partitions
#
# Control: the same three questions asked of a deliberately wrong triple, (d1, d2, d1+d2), which is
# the size-three relation imposed everywhere.  It must behave differently, or the measurement is not
# reading the actual image.
#
# Authors: Carles Marin, Claude (AI assistant).

Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()


def dtriple(lam, t):
    """(d1, d2, d3) ordered as the theorem builds them: r_A <= r_B.  None if degenerate."""
    N = t + 2
    lam = list(lam) + [0] * (N - len(lam))
    if len(lam) > N:
        return None
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for b in beta:
        cls.setdefault(b % t, []).append(b)
    prof = [len(cls.get(i, [])) for i in range(t)]
    if 0 in prof:
        return None
    big = sorted([i for i in range(t) if prof[i] >= 2])
    if len(big) == 2:
        rA, rB = big
        A = sorted(cls[rA], reverse=True)[:2]
        B = sorted(cls[rB], reverse=True)[:2]
    else:
        P = sorted(cls[big[0]], reverse=True)
        A, B = [P[0], P[1]], [P[1], P[2]]
    return (A[0] - A[1], B[0] - B[1], abs(A[0] + A[1] - B[0] - B[1]))


print("=" * 76)
print("Is the S_3 of the formula realised on partitions?")
print("=" * 76)

for t, MAX in ((2, 26), (3, 22), (4, 20)):
    N = t + 2
    seen = {}
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            d = dtriple(list(l), t)
            if d is None or 0 in d:
                continue
            seen.setdefault(d, []).append(size)

    attained = set(seen)
    full, partial, sizeok = 0, 0, 0
    missing_example = None
    for d in attained:
        perms = set(Permutations(list(d)).list())
        perms = set(tuple(x) for x in perms)
        if perms <= attained:
            full += 1
            # can a permutation be reached at the SAME |lambda|?
            for q in perms:
                if q != d and set(seen[d]) & set(seen[q]):
                    sizeok += 1
                    break
        else:
            partial += 1
            if missing_example is None:
                missing_example = (d, sorted(perms - attained)[:2])

    print("")
    print("  t = %d,  N = %d,  |lambda| <= %d" % (t, N, MAX))
    print("    attained triples          : %d" % len(attained))
    print("    with every permutation too: %d" % full)
    print("    missing some permutation  : %d" % partial)
    print("    a permutation reachable at the SAME |lambda|: %d of %d" % (sizeok, full))
    if missing_example:
        print("    e.g. %s attained, %s not" % (missing_example[0], missing_example[1]))

    # control: impose the size-three relation everywhere
    ctrl = set((a, b, a + b) for (a, b, _) in attained)
    print("    CONTROL, triples forced to d3 = d1+d2: %d, of which attained: %d"
          % (len(ctrl), len(ctrl & attained)))

print("")
print("DONE")
