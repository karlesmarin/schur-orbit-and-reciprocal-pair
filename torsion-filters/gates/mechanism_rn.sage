# -*- coding: utf-8 -*-
# The mechanism at general r, checked before it is believed.
#
# CLAIM.  The 2r x 2r alternant A(T) = det[ z_j^{b}, z_j^{-b} ]_{b in T} has its columns in PAIRS,
# one pair per variable.  A determinant whose columns come in blocks expands by Laplace over the
# blocks, each block taking two rows, and the 2x2 minor of the pair (z_j^b, z_j^{-b}) on rows
# {b, b'} is z_j^{b-b'} - z_j^{b'-b} = 2 sinh((b-b') u_j), the same single sine that carries the
# whole r = 1 case.  Hence
#
#     A(T)  =  sum over the ways to split T into r ORDERED pairs, one per variable,
#              of  +-  prod_j  2 sinh( (b_{j,1} - b_{j,2}) u_j ).
#
# Two consequences, if it holds:
#   * r = 1 has a single matching, so the four transversal terms telescope to the three factors of
#     Theorem 3.1.  For r >= 2 there is a SUM over matchings and nothing telescopes -- which is
#     Conjecture 10.4 seen as a mechanism rather than as evidence.
#   * under b -> C - b every difference negates and sinh is odd, so every matching term picks up
#     (-1)^r.  That is A(C-T) = (-1)^r A(T) with a two-line combinatorial reason instead of a
#     column-factoring one, and it is what makes the zero locus alternate with the parity of r.
#
# Checked here symbolically for r = 1, 2, 3 on random row sets, against the honest determinant.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools


def check(r, trials, seed):
    N = 2 * r
    R = LaurentPolynomialRing(QQ, ['z%d' % j for j in range(r)])
    zs = list(R.gens())
    set_random_seed(seed)
    bad = 0
    for _ in range(trials):
        T = []
        while len(T) < N:
            x = ZZ.random_element(0, 14)
            if x not in T:
                T.append(x)
        T = sorted(T, reverse=True)
        M = matrix(R, N, N,
                   lambda a, b: (zs[b // 2] ** T[a]) if b % 2 == 0 else (zs[b // 2] ** (-T[a])))
        lhs = M.det()

        # the matching expansion: split the row positions into r ordered pairs, pair j to variable j
        rhs = R(0)
        idx = list(range(N))
        for perm in itertools.permutations(idx):
            # take pairs (perm[0],perm[1]) -> var 0, (perm[2],perm[3]) -> var 1, ...
            pairs = [(perm[2 * j], perm[2 * j + 1]) for j in range(r)]
            # canonical: first element of each pair smaller index, pairs sorted by first element,
            # so each unordered matching is generated exactly once
            if any(p[0] > p[1] for p in pairs):
                continue
            if any(pairs[j][0] > pairs[j + 1][0] for j in range(r - 1)):
                pass  # the assignment of pairs to VARIABLES matters, so do not sort them
            sg = Permutation([x + 1 for x in perm]).signature()
            term = R(sg)
            for j, (u, v) in enumerate(pairs):
                term *= zs[j] ** (T[u] - T[v]) - zs[j] ** (T[v] - T[u])
            rhs += term
        if lhs != rhs:
            bad += 1
            if bad == 1:
                print("    MISMATCH at r=%d, rows %s" % (r, T))
    return bad


print("=" * 74)
print("The matching expansion of the alternant, checked against the determinant")
print("=" * 74)
print("")
for r, trials in ((1, 6), (2, 5), (3, 3)):
    bad = check(r, trials, 100 + r)
    print("  r = %d : %d mismatches in %d random row sets%s"
          % (r, bad, trials, "   <-- holds" if bad == 0 else ""))

print("")
print("  CONTROL: the same expansion with one sinh replaced by a cosh must FAIL, or the check")
print("           above is not testing anything.")
r = 2
N = 4
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
zs = list(R.gens())
T = [9, 5, 2, 0]
M = matrix(R, N, N, lambda a, b: (zs[b // 2] ** T[a]) if b % 2 == 0 else (zs[b // 2] ** (-T[a])))
lhs = M.det()
rhs = R(0)
for perm in itertools.permutations(range(N)):
    pairs = [(perm[0], perm[1]), (perm[2], perm[3])]
    if any(p[0] > p[1] for p in pairs):
        continue
    sg = Permutation([x + 1 for x in perm]).signature()
    term = R(sg)
    u, v = pairs[0]
    term *= zs[0] ** (T[u] - T[v]) + zs[0] ** (T[v] - T[u])          # cosh instead of sinh
    u, v = pairs[1]
    term *= zs[1] ** (T[u] - T[v]) - zs[1] ** (T[v] - T[u])
    rhs += term
print("           deliberately wrong version agrees with the determinant: %s  (want False)"
      % (lhs == rhs))

print("")
print("DONE")
