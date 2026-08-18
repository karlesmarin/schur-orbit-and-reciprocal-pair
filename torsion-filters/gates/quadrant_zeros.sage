# -*- coding: utf-8 -*-
# The interior of the (t, r) map: where does s_lambda(mu_t, z_1^{+-1}, ..., z_r^{+-1}) vanish?
#
# The paper settles the line r = 1 for every t, and the line t = 2 for every r.  The interior
# t >= 3, r >= 2 it says nothing about.  Both settled lines share their FIRST branch -- the value
# vanishes when a residue class mod t is empty -- and differ in the second: concentric intervals at
# r = 1 (even t only), self-complementary of odd width at t = 2.  So the question with a shape is:
# does branch one survive into the interior, and what replaces branch two?
#
# The alphabet's power sums are rational, p_k = t*[t|k] + sum_j (z_j^k + z_j^{-k}), because the
# t-th roots of unity sum to t when t | k and to 0 otherwise.  No cyclotomic field is needed and
# Psi lives in Z[z^{+-1}].
#
# Measured per (t, r):
#   empty class -> zero        must be 100%, it is branch one
#   EXTRA zeros                vanishing with no empty class: this is where the new criterion lives
#   and for each extra zero, its profile, t-core size, t-quotient and two candidate predictors:
#     SC   self-complementary of odd width, which is the t = 2 branch
#     CONC the concentric-interval condition of the r = 1 branch, in the only form that makes sense
#          above excess 2: SOME pair of excess classes concentric as intervals
#
# CONTROL: at r = 1 the machinery must reproduce the paper's criterion exactly, 0 false positives
# and 0 false negatives.  Without that line nothing below is evidence of anything.
#
# Authors: Carles Marin, Claude (AI assistant).

def alphabet_power(k, t, zs):
    return (t if k % t == 0 else 0) + sum(zz ** k + zz ** (-k) for zz in zs)


def psi(lam, t, R, zs, p, s):
    out = R(0)
    for rho, c in p(s[list(lam)] if lam else s.one()):
        term = R(c)
        for k in rho:
            term *= R(alphabet_power(k, t, zs))
        out += term
    return out


def profile(lam, t, N):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for b in beta:
        cls.setdefault(b % t, []).append(b)
    return beta, [sorted(cls.get(i, []), reverse=True) for i in range(t)]


def core_quotient(cl, t):
    quot = []
    for i in range(t):
        e = cl[i]
        n = len(e)
        quot.append(tuple((b - i) // t - (n - 1 - j) for j, b in enumerate(e)))
    return tuple(quot)


def self_comp_odd_width(lam):
    lam = [x for x in lam if x > 0]
    if not lam:
        return False
    h, w = len(lam), lam[0]
    if w % 2 == 0:
        return False
    return all(lam[i] + lam[h - 1 - i] == w for i in range(h))


def concentric(cl, t):
    """Some pair of excess classes nested as intervals [min, max]."""
    big = [c for c in cl if len(c) >= 2]
    for i in range(len(big)):
        for j in range(i + 1, len(big)):
            a, b = big[i], big[j]
            lo1, hi1, lo2, hi2 = a[-1], a[0], b[-1], b[0]
            if (lo1 < lo2 and hi2 < hi1) or (lo2 < lo1 and hi1 < hi2):
                return True
    return False


Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()

print("=" * 78)
print("The (t, r) interior: zeros of s_lambda(mu_t, z_1^{+-1}, ..., z_r^{+-1})")
print("=" * 78)
print("")
print("  t  r   N  |lam|<=  shapes  zeros  empty-class  EXTRA  of extra: SC  CONC")
print("  " + "-" * 72)

extras = []
for t, r, MAX in ((3, 1, 14), (4, 1, 14),
                  (3, 2, 11), (4, 2, 11),
                  (3, 3, 9), (4, 3, 9)):
    N = t + 2 * r
    R = LaurentPolynomialRing(QQ, ['z%d' % j for j in range(r)])
    zs = list(R.gens())
    nsh = nz = nempty = nextra = nsc = nconc = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            lam = list(l)
            beta, cl = profile(lam, t, N)
            empty = any(len(c) == 0 for c in cl)
            v = psi(lam, t, R, zs, p, s)
            nsh += 1
            if v == 0:
                nz += 1
                if empty:
                    nempty += 1
                else:
                    nextra += 1
                    sc = self_comp_odd_width(lam)
                    co = concentric(cl, t)
                    nsc += sc
                    nconc += co
                    if len(extras) < 400:
                        extras.append((t, r, tuple(lam), tuple(len(c) for c in cl),
                                       core_quotient(cl, t), sc, co))
            elif empty:
                print("    !! empty class but NONZERO: t=%d r=%d lam=%s" % (t, r, lam))
    print("  %2d %2d %3d %7d %7d %6d %12d %6d %13d %5d"
          % (t, r, N, MAX, nsh, nz, nempty, nextra, nsc, nconc))

print("")
print("  the extra zeros, in full:")
print("    t  r  lambda                profile        SC  CONC  t-quotient")
for t, r, lam, prof, q, sc, co in extras[:40]:
    print("    %d  %d  %-20s %-14s %-3s %-4s %s"
          % (t, r, str(list(lam)), str(list(prof)), "yes" if sc else ".",
             "yes" if co else ".", str([list(x) for x in q])))
if len(extras) > 40:
    print("    ... and %d more" % (len(extras) - 40))

print("")
print("DONE")
