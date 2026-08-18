# -*- coding: utf-8 -*-
# The same quadrant, with the control the first run was missing, and a longer range.
#
# quadrant_zeros.sage counted how many EXTRA zeros are concentric and got 11 of 11 at (t,r)=(4,1)
# and 1 of 1 at (4,2) and (4,3).  That is only half a measurement: if being concentric is common
# among NON-vanishing shapes too, the condition predicts nothing and those ratios are an artefact of
# only ever looking at the zeros.  So this counts the full 2x2 table --
#
#            vanishes   does not
#   conc        a           b        <- b must be 0 for the condition to be sufficient
#   not conc    c           d        <- c must be 0 (given no empty class) for it to be necessary
#
# -- and pushes |lambda| further, because at r >= 2 the whole finding rested on one or two zeros.
#
# Authors: Carles Marin, Claude (AI assistant).

def psi(lam, t, R, zs, p, s):
    out = R(0)
    for rho, c in p(s[list(lam)] if lam else s.one()):
        term = R(c)
        for k in rho:
            term *= R((t if k % t == 0 else 0) + sum(zz ** k + zz ** (-k) for zz in zs))
        out += term
    return out


def classes(lam, t, N):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for b in beta:
        cls.setdefault(b % t, []).append(b)
    return [sorted(cls.get(i, []), reverse=True) for i in range(t)]


def concentric(cl):
    """Two classes CONCENTRIC as intervals: the same centre, min+max equal.

    Not nested.  The first version of this function tested nesting, and the control caught it at
    once: nesting left 116 concentric non-vanishing shapes at (t,r)=(4,1), where the paper's
    criterion is an equivalence.  Concentric means concentric -- (a_1+a_2)/2 = (b_1+b_2)/2, which
    on the excess-2 line is exactly d_3 = 0."""
    big = [c for c in cl if len(c) >= 2]
    for i in range(len(big)):
        for j in range(i + 1, len(big)):
            a, b = big[i], big[j]
            if a[0] + a[-1] == b[0] + b[-1]:
                return True
    return False


Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()

print("=" * 78)
print("Is 'concentric' a criterion, or only a property of the zeros I looked at?")
print("=" * 78)
print("")
print("  among shapes with NO empty class:")
print("   t  r  |lam|<=  shapes   conc&zero  conc&nonzero  notconc&zero  notconc&nonzero")
print("  " + "-" * 76)

for t, r, MAX in ((3, 1, 18), (4, 1, 18),
                  (3, 2, 15), (4, 2, 15),
                  (3, 3, 11), (4, 3, 11)):
    N = t + 2 * r
    R = LaurentPolynomialRing(QQ, ['z%d' % j for j in range(r)])
    zs = list(R.gens())
    a = b = c = d = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            lam = list(l)
            cl = classes(lam, t, N)
            if any(len(x) == 0 for x in cl):
                continue
            z = (psi(lam, t, R, zs, p, s) == 0)
            k = concentric(cl)
            if k and z:
                a += 1
            elif k and not z:
                b += 1
            elif z:
                c += 1
            else:
                d += 1
    print("  %2d %2d %7d %7d %11d %13d %13d %16d"
          % (t, r, MAX, a + b + c + d, a, b, c, d))

print("")
print("  b > 0 means concentric does NOT imply vanishing.")
print("  c > 0 means vanishing without an empty class and without concentric: a THIRD branch.")
print("")
print("DONE")
