# -*- coding: utf-8 -*-
# Does the value collapse, and does the collapse survive above excess 2?
#
# This replaces a first attempt whose signature -- (profile, within-class gaps, cross offset, sign)
# -- turned out to be INJECTIVE on the shapes it was tested over: it separated every partition, so
# its "0 conflicts" measured nothing. 125 signature classes over 125 shapes is not a result.
# A datum that is a repackaging of beta always determines the value.  The question worth asking is
# the opposite one, and it needs no datum at all:
#
#     how much does lambda |-> Psi_r(lambda) collapse?
#
# At excess 2 the collapse is enormous and that is the content of the paper: Theorem 3.1 says the
# value factors through (d_1,d_2,d_3,eps), and Corollary "the fibre of the triple" says each fibre
# is INFINITE.  So if the graph reading extends above excess 2, the collapse must survive; if the
# map becomes essentially injective at excess 4, then compression is special to excess 2, which is
# Conjecture 10.4 seen from the other side.
#
# Measured, for each r, over every shape in the range:
#   shapes          how many partitions
#   values          how many distinct Psi_r
#   values up to +- how many distinct |Psi_r|, since the sign is a separate factor
#   biggest fibre   the most shapes sharing one value (zero excluded: it is its own story)
#
# Control that the instrument sees a collapse when there is one: at r = 1 the paper's own triple
# (d_1,d_2,d_3) must have 0 conflicts up to sign.  If that column is not 0 the measurement is
# broken and no other number here means anything.
#
# Authors: Carles Marin, Claude (AI assistant).

Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()


def psi(lam, r, R, zs):
    out = R(0)
    for rho, c in p(s[list(lam)] if lam else s.one()):
        term = R(c)
        for k in rho:
            term *= R(1 + (-1) ** k + sum(zz ** k + zz ** (-k) for zz in zs))
        out += term
    return out


def triple(lam, N):
    """The paper's interval triple, at t = 2.  None if the profile is degenerate."""
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for b in beta:
        cls.setdefault(b % 2, []).append(b)
    prof = [len(cls.get(i, [])) for i in range(2)]
    if 0 in prof:
        return None
    big = sorted([i for i in range(2) if prof[i] >= 2])
    if len(big) == 2:
        A = sorted(cls[big[0]], reverse=True)[:2]
        B = sorted(cls[big[1]], reverse=True)[:2]
    else:
        P = sorted(cls[big[0]], reverse=True)
        A, B = [P[0], P[1]], [P[1], P[2]]
    return (A[0] - A[1], B[0] - B[1], abs(A[0] + A[1] - B[0] - B[1]))


print("=" * 78)
print("Does the collapse survive above excess 2?")
print("=" * 78)
print("")
print("  r  N  exc  |lam|<=  shapes  values  up to +-  biggest fibre   zeros")
print("  " + "-" * 68)

for r, MAX in ((1, 14), (2, 12), (3, 9)):
    N = 2 * r + 2
    R = LaurentPolynomialRing(QQ, ['z%d' % j for j in range(r)])
    zs = list(R.gens())
    shapes = [(list(l), psi(list(l), r, R, zs))
              for size in range(0, MAX + 1)
              for l in Partitions(size, max_length=N)]
    nz = [(lam, v) for lam, v in shapes if v != 0]
    vals = set(v for _, v in nz)
    absv = {}
    for lam, v in nz:
        absv.setdefault(min(v, -v), []).append(lam)
    big = max(len(x) for x in absv.values())
    print("  %d  %d  %3d  %6d  %6d  %6d  %8d  %13d  %6d"
          % (r, N, N - 2, MAX, len(shapes), len(vals), len(absv), big,
             len(shapes) - len(nz)))

print("")
print("  CONTROL, r = 1: the paper's triple (d1,d2,d3) must give 0 conflicts up to sign.")
R = LaurentPolynomialRing(QQ, ['z0'])
z0 = R.gen()
bk = {}
for size in range(0, 15):
    for l in Partitions(size, max_length=4):
        d = triple(list(l), 4)
        if d is None:
            continue
        v = psi(list(l), 1, R, [z0])
        bk.setdefault(d, set()).add(min(v, -v))
conf = sum(1 for vs in bk.values() if len(vs) > 1)
print("      %d triples, %d conflicts, biggest fibre %d shapes"
      % (len(bk), conf, max(len(vs) for vs in bk.values())))

print("")
print("DONE")
