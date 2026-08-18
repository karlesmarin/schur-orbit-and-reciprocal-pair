# -*- coding: utf-8 -*-
# The criterion, stated on beta alone.
#
# The Laplace decomposition says Psi ~ sum_T sgn(T) A(T) with T running over the complements of
# transversals, and A(C-T) = (-1)^r A(T) because reflecting a row b -> C-b factors z_j^{C} out of
# one column and z_j^{-C} out of the other and swaps them, r swaps in all.  So for r even the sum
# collapses when the term set is closed under T -> C-T with opposite signs.
#
# But the T's are built only from the EXCESS classes: a class of size one gives its single element
# to every transversal and never appears in T.  So closure of the term set should be exactly
#
#     the UNION OF THE EXCESS CLASSES is invariant under b -> C - b
#
# for some C, which is then a common centre.  That is a condition on beta alone, with no
# determinant, no alternant and no evaluation in it.  It contains both known branches: two classes
# of size two invariant about C is concentric, and beta invariant about C is self-complementary.
#
# This checks that statement against the value, on the full concentric family and on the shapes
# with no concentric pair as well, so the claim is tested where it should FAIL too.
#
# Authors: Carles Marin, Claude (AI assistant).

t, r = 4, 2
N = t + 2 * r
MAX = 26
PRIME = 998244353
F = GF(PRIME)
I = F(-1).sqrt()
ROOTS = [F(1), I, F(-1), -I]
set_random_seed(2024)
POINTS = []
while len(POINTS) < 3:
    z1, z2 = F.random_element(), F.random_element()
    if z1 and z2:
        alph = ROOTS + [z1, 1 / z1, z2, 1 / z2]
        if len(set(alph)) == N:
            POINTS.append(alph)


def vanishes(beta):
    for alph in POINTS:
        if matrix(F, N, N, lambda a, b: alph[b] ** beta[a]).det() != 0:
            return False
    return True


def excess_union(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for b in beta:
        cls.setdefault(b % t, []).append(b)
    if len(cls) < t:
        return None, None
    U = sorted([b for k, v in cls.items() if len(v) >= 2 for b in v])
    return beta, U


def reflective(U):
    """Is the union of the excess classes invariant under b -> C - b for some C?"""
    S = set(U)
    for C in range(2 * max(U) + 1):
        if set(C - b for b in U) == S:
            return C
    return None


print("=" * 78)
print("The criterion on beta alone,  t = 4, r = 2,  |lambda| <= %d" % MAX)
print("=" * 78)

a = b = c = d = 0
false_pos, false_neg = [], []
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        lam = list(l)
        beta, U = excess_union(lam)
        if beta is None:
            continue
        pred = reflective(U) is not None
        z = vanishes(beta)
        if pred and z:
            a += 1
        elif pred and not z:
            b += 1
            if len(false_pos) < 6:
                false_pos.append((lam, beta, U))
        elif z:
            c += 1
            if len(false_neg) < 6:
                false_neg.append((lam, beta, U))
        else:
            d += 1

print("")
print("  shapes with no empty class: %d" % (a + b + c + d))
print("")
print("                              vanishes    does not")
print("  excess union reflective    %9d %11d" % (a, b))
print("  not reflective             %9d %11d" % (c, d))
print("")
print("  false positives (predicted zero, is not): %d" % b)
for lam, beta, U in false_pos:
    print("    lam=%-24s beta=%-30s excess=%s" % (str(lam), str(beta), str(U)))
print("  false negatives (vanishes, not predicted): %d" % c)
for lam, beta, U in false_neg:
    print("    lam=%-24s beta=%-30s excess=%s" % (str(lam), str(beta), str(U)))

print("")
print("DONE")
