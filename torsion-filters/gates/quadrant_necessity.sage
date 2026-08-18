# -*- coding: utf-8 -*-
# The attack on my own sweep: I filtered by concentric BEFORE computing.
#
# quadrant_big.sage swept only the concentric shapes, on the strength of a necessity established
# on a much shorter range (|lambda| <= 15).  That makes the big sweep unable, by construction, to
# see a vanishing shape that is not concentric.  If one exists above 15, the necessity is false and
# every conclusion drawn from that run inherits the error.
#
# So: no filter.  Every shape with no empty class, tested with the same finite-field bialternant,
# and the full 2x2 table.  Also counts the shapes that are beta-symmetric about SOME centre without
# being concentric, since those are precisely where a counterexample would hide -- the sufficient
# condition found in the big sweep would then predict a non-concentric zero.
#
# Authors: Carles Marin, Claude (AI assistant).

t, r = 4, 2
N = t + 2 * r
MAX = 28
PRIME = 998244353
F = GF(PRIME)
I = F(-1).sqrt()
ROOTS = [F(1), I, F(-1), -I]

set_random_seed(999)
POINTS = []
while len(POINTS) < 3:
    z1, z2 = F.random_element(), F.random_element()
    if z1 == 0 or z2 == 0:
        continue
    alph = ROOTS + [z1, 1 / z1, z2, 1 / z2]
    if len(set(alph)) == N:
        POINTS.append(alph)


def vanishes(beta):
    for alph in POINTS:
        if matrix(F, N, N, lambda a, b: alph[b] ** beta[a]).det() != 0:
            return False
    return True


def info(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = [sorted([b for b in beta if b % t == i], reverse=True) for i in range(t)]
    return beta, cls


def concentric(cls):
    big = [c for c in cls if len(c) >= 2]
    for x in range(len(big)):
        for y in range(x + 1, len(big)):
            if big[x][0] + big[x][-1] == big[y][0] + big[y][-1]:
                return True
    return False


def beta_symmetric(beta):
    """Symmetric about ANY centre, not only about a concentric class's centre."""
    S = set(beta)
    for C in range(2 * max(beta) + 1):
        if set(C - b for b in beta) == S:
            return C
    return None


print("=" * 78)
print("Necessity of 'concentric', with no pre-filter,  t=4 r=2,  |lambda| <= %d" % MAX)
print("=" * 78)

a = b = c = d = 0
rogue = []
symnotconc = 0
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        lam = list(l)
        beta, cls = info(lam)
        if any(len(x) == 0 for x in cls):
            continue
        k = concentric(cls)
        z = vanishes(beta)
        if k and z:
            a += 1
        elif k and not z:
            b += 1
        elif z:
            c += 1
            if len(rogue) < 12:
                rogue.append((lam, beta))
        else:
            d += 1
        if not k and beta_symmetric(beta) is not None:
            symnotconc += 1

print("")
print("  shapes with no empty class : %d" % (a + b + c + d))
print("")
print("             vanishes   does not")
print("  conc      %8d %10d" % (a, b))
print("  not conc  %8d %10d   <- c must be 0 or the necessity is false" % (c, d))
print("")
print("  beta-symmetric but NOT concentric: %d shapes" % symnotconc)
print("  (that is where a counterexample would hide: the sufficient condition would predict")
print("   a vanishing shape outside the concentric family)")
for lam, beta in rogue:
    print("    ROGUE zero, not concentric: lam=%s beta=%s" % (lam, beta))

print("")
print("DONE")
