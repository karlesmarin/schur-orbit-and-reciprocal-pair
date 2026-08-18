# -*- coding: utf-8 -*-
# What IS the invariant above excess 2?  Measured against the true fibres, not guessed.
#
# The metric is the point of this script.  A candidate datum is graded on two numbers at once:
#
#   conflicts   shapes sharing the datum but not the value.  Must be 0, or the datum is too coarse.
#   classes     how many distinct data.  Compared with the number of distinct VALUES: any excess is
#               collapse the datum fails to see, so the datum is finer than the invariant.
#
# The ideal datum has conflicts = 0 AND classes = values.  My first attempt scored 0 conflicts and
# 125 classes against 76 values, which reads as success and is not: it was finer than the value and
# in that range injective.  Grading on one number is what let that pass.
#
# Positive control, r = 1: the paper's own datum is the MULTISET {d1,d2,d3} with the sign, unordered
# because the closed form is symmetric in the three arguments.  It must score conflicts 0 and
# classes = values.  If it does not, the metric is broken and nothing else here counts.
#
# Values are taken up to sign throughout and the zero value is excluded, so the sign is not part of
# any datum below.
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


def classes_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {0: [], 1: []}
    for b in beta:
        cls[b % 2].append(b)
    out = []
    for i in (0, 1):
        e = sorted(cls[i], reverse=True)
        n = len(e)
        gaps = tuple(e[j] - e[j + 1] for j in range(n - 1))
        quot = tuple((e[j] - i) // 2 - (n - 1 - j) for j in range(n))
        out.append(dict(res=i, n=n, elems=e, gaps=gaps, quot=quot, sum=sum(e)))
    return out


def triple(cl):
    """The paper's interval triple at t = 2, as an unordered multiset."""
    big = [c for c in cl if c['n'] >= 2]
    if any(c['n'] == 0 for c in cl):
        return None
    if len(big) == 2:
        A, B = big[0]['elems'][:2], big[1]['elems'][:2]
    else:
        P = big[0]['elems']
        A, B = [P[0], P[1]], [P[1], P[2]]
    return tuple(sorted([A[0] - A[1], B[0] - B[1], abs(A[0] + A[1] - B[0] - B[1])]))


# the ladder of candidates, all symmetric under swapping the two residue classes
CAND = {
    "multiset of gaps, all classes":
        lambda cl: tuple(sorted(c['gaps'] for c in cl)),
    "+ |cross|":
        lambda cl: (tuple(sorted(c['gaps'] for c in cl)),
                    abs(cl[0]['sum'] - cl[1]['sum'])),
    "multiset of gaps, EXCESS classes only, + |cross|":
        lambda cl: (tuple(sorted(c['gaps'] for c in cl if c['n'] >= 2)),
                    abs(cl[0]['sum'] - cl[1]['sum'])),
    "multiset of quotient components, excess classes":
        lambda cl: tuple(sorted(c['quot'] for c in cl if c['n'] >= 2)),
    "multiset of quotient components, all classes":
        lambda cl: tuple(sorted(c['quot'] for c in cl)),
}

print("=" * 78)
print("Grading candidate data against the true fibres")
print("=" * 78)

for r, MAX in ((1, 16), (2, 14)):
    N = 2 * r + 2
    R = LaurentPolynomialRing(QQ, ['z%d' % j for j in range(r)])
    zs = list(R.gens())
    rows = []
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            v = psi(list(l), r, R, zs)
            if v == 0:
                continue
            rows.append((list(l), min(v, -v), classes_of(list(l), N)))
    nvals = len(set(v for _, v, _ in rows))
    print("")
    print("  r = %d,  N = %d,  excess = %d,  |lambda| <= %d" % (r, N, N - 2, MAX))
    print("    %d nonzero shapes, %d distinct values up to sign" % (len(rows), nvals))
    print("    %-50s %9s %8s" % ("datum", "conflicts", "classes"))
    print("    " + "-" * 70)

    def grade(name, f):
        bk = {}
        for lam, v, cl in rows:
            k = f(cl)
            if k is None:
                continue
            bk.setdefault(k, set()).add(v)
        conf = sum(1 for vs in bk.values() if len(vs) > 1)
        flag = ""
        if conf == 0:
            flag = "  <-- exact" if len(bk) == nvals else "  (too fine by %d)" % (len(bk) - nvals)
        print("    %-50s %9d %8d%s" % (name, conf, len(bk), flag))

    if r == 1:
        grade("CONTROL: multiset {d1,d2,d3} (the paper's)", triple)
    for name, f in CAND.items():
        grade(name, f)

print("")
print("DONE")
