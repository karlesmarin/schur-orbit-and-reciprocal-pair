# -*- coding: utf-8 -*-
# Trying to break the previous conclusion, before believing it.
#
# excess_s3_symmetry.sage reported that most attained triples are missing some permutation, and I
# was about to read that as "the S_3 of the formula is not realised on partitions".  There is an
# obvious way for that to be wrong: a permuted triple may exist and simply need a larger |lambda|
# than the cutoff.  Then the number would be a fact about my sweep and not about the object.
#
# Two attacks, and they have different characters.
#
#   A -- RANGE.  Raise the cutoff and watch "missing" as a fraction.  If it falls steadily the
#        conclusion is an artefact.  If it settles, it is not.
#
#   B -- ARITHMETIC, which no range can rescue.  Both a_1,a_2 lie in one residue class mod t and
#        b_1,b_2 in another, so d_1 = a_1-a_2 and d_2 = b_1-b_2 are MULTIPLES OF t, while
#        d_3 = |a_1+a_2-b_1-b_2| = 2(r_A-r_B) mod t.  So whenever t does not divide d_3, no
#        partition whatsoever has that value in slot 1 or 2, at any size.  This counts how often
#        that happens, and checks the two congruences on every attained triple rather than assuming
#        them.
#
# Authors: Carles Marin, Claude (AI assistant).

def dtriple(lam, t):
    N = t + 2
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
        rA, rB = big
        A = sorted(cls[rA], reverse=True)[:2]
        B = sorted(cls[rB], reverse=True)[:2]
    else:
        rA = rB = big[0]
        P = sorted(cls[big[0]], reverse=True)
        A, B = [P[0], P[1]], [P[1], P[2]]
    return (A[0] - A[1], B[0] - B[1], abs(A[0] + A[1] - B[0] - B[1])), rA, rB


print("=" * 74)
print("A -- is 'missing a permutation' an artefact of the cutoff?")
print("=" * 74)

for t, cutoffs in ((2, (14, 18, 22, 26, 30, 34, 38)), (4, (12, 16, 20, 24, 28))):
    N = t + 2
    print("")
    print("  t = %d" % t)
    print("    |lam|<=   triples   full orbits   missing   missing %")
    print("    " + "-" * 56)
    seen = set()
    done = 0
    for hi in cutoffs:
        for size in range(done, hi + 1):
            for l in Partitions(size, max_length=N):
                r = dtriple(list(l), t)
                if r is None or 0 in r[0]:
                    continue
                seen.add(r[0])
        done = hi + 1
        full = 0
        for d in seen:
            perms = set(tuple(x) for x in Permutations(list(d)).list())
            if perms <= seen:
                full += 1
        miss = len(seen) - full
        print("    %6d   %7d   %11d   %7d   %8.1f%%"
              % (hi, len(seen), full, miss, 100.0 * miss / len(seen)))

print("")
print("=" * 74)
print("B -- the arithmetic obstruction, checked and not assumed")
print("=" * 74)

for t, MAX in ((2, 24), (3, 22), (4, 20)):
    N = t + 2
    tot = bad12 = bad3 = blocked = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            r = dtriple(list(l), t)
            if r is None or 0 in r[0]:
                continue
            (d1, d2, d3), rA, rB = r
            tot += 1
            bad12 += (d1 % t != 0 or d2 % t != 0)
            bad3 += ((d3 - 2 * (rA - rB)) % t != 0)
            blocked += (d3 % t != 0)
    print("  t=%d : %5d triples | d1,d2 not multiples of t: %d | d3 off 2(rA-rB) mod t: %d"
          % (t, tot, bad12, bad3))
    print("        of those, %d (%.0f%%) have t not dividing d3, so a permutation moving d3 into"
          % (blocked, 100.0 * blocked / tot))
    print("        slot 1 or 2 is impossible at ANY size, not merely absent from this range.")

print("")
print("DONE")
