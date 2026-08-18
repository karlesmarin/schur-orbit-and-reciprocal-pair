# -*- coding: utf-8 -*-
# Two of my scripts disagree.  One case, printed term by term, decides which.
#
# quadrant_mechanism.sage reported that for the vanishing shapes the term set is closed under
# T -> C-T with every matched pair carrying opposite signs, 30 of 30 and 0 of 635.
# quadrant_parity.sage computes the same sign ratio term by term and finds it NOT constant on 93
# shapes and equal to -1 on only 2.  Both cannot be right.
#
# lambda = (5,4,3), beta = [12,10,8,4,3,2,1,0], C = 12 is one of the 30.  Print every transversal,
# its complement T, its Laplace sign, the reflected C-T, and that one's sign, and add the signed
# alternants up.  Whatever the table says is the answer.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

t, r = 4, 2
N = 8
lam = [5, 4, 3]
lam = lam + [0] * (N - len(lam))
beta = [lam[i] + N - 1 - i for i in range(N)]
C = 12
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
zs = list(R.gens())

cls = {}
for k, b in enumerate(beta):
    cls.setdefault(b % t, []).append(k)
pos = dict((b, i) for i, b in enumerate(beta))


def shuffle_sign(pick, rest):
    perm = list(pick) + list(rest)
    sg = 1
    for a in range(len(perm)):
        for b in range(a + 1, len(perm)):
            if perm[a] > perm[b]:
                sg = -sg
    return sg


def alt(T):
    return matrix(R, len(T), len(T),
                  lambda a, b: (zs[b // 2] ** T[a]) if b % 2 == 0
                  else (zs[b // 2] ** (-T[a]))).det()


print("=" * 78)
print("lambda = %s   beta = %s   C = %d" % (lam, beta, C))
print("classes: %s" % {k: [beta[i] for i in v] for k, v in cls.items()})
print("=" * 78)
print("")
print("   T (values)             sgn   C-T (values)          sgn   in beta?   A(T)=A(C-T)?")
print("   " + "-" * 74)

total = R(0)
for pick in itertools.product(*[cls[k] for k in sorted(cls)]):
    rest = tuple(sorted(set(range(N)) - set(pick)))
    T = [beta[i] for i in rest]
    sg = shuffle_sign(sorted(pick), rest)
    total += sg * alt(T)
    TR = [C - x for x in T]
    inside = all(x in pos for x in TR)
    if inside:
        rest2 = tuple(sorted(pos[x] for x in TR))
        pick2 = tuple(sorted(set(range(N)) - set(rest2)))
        sg2 = shuffle_sign(pick2, rest2)
        same_alt = (alt(T) == alt(sorted(TR, reverse=True)))
    else:
        sg2, same_alt = None, None
    print("   %-22s %+3d   %-22s %s   %-8s   %s"
          % (str(T), sg, str(TR), ("%+3d" % sg2) if sg2 is not None else "  -",
             "yes" if inside else "NO", str(same_alt)))

print("")
print("   signed total is zero: %s" % (total == 0))
print("")
print("DONE")
