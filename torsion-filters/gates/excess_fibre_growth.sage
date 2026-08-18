# -*- coding: utf-8 -*-
# Does the fibre GROW?  The range-free form of the question.
#
# Counting distinct values at one cutoff cannot compare excess 2 with excess 4: the alphabets have
# different numbers of variables, so the higher one has more room to be distinct for reasons that
# have nothing to do with the mathematics.  What does compare is the GROWTH.
#
# At excess 2 the paper proves the fibre of a value is infinite (Corollary "the fibre of the
# triple": the free parameters are an integer v >= max(0,c) and t-2 further integers).  So as the
# cutoff on |lambda| rises, the largest fibre must rise with it, without bound.  If at excess 4 the
# largest fibre saturates instead, the collapse is a phenomenon of excess 2 and not of the graph --
# and that is Conjecture 10.4 measured rather than argued.
#
# Fibres are taken up to sign, since the sign is a separate factor of the value, and the zero value
# is excluded: it is the vanishing locus and it has its own criterion.
#
# Authors: Carles Marin, Claude (AI assistant).

Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()


def psi(lam, r, R, zs, cache={}):
    key = (tuple(lam), r)
    if key in cache:
        return cache[key]
    out = R(0)
    for rho, c in p(s[list(lam)] if lam else s.one()):
        term = R(c)
        for k in rho:
            term *= R(1 + (-1) ** k + sum(zz ** k + zz ** (-k) for zz in zs))
        out += term
    cache[key] = out
    return out


print("=" * 74)
print("Growth of the largest fibre with the cutoff")
print("=" * 74)

for r, cutoffs in ((1, (6, 8, 10, 12, 14, 16, 18)), (2, (6, 8, 10, 12, 14, 16))):
    N = 2 * r + 2
    R = LaurentPolynomialRing(QQ, ['z%d' % j for j in range(r)])
    zs = list(R.gens())
    print("")
    print("  r = %d,  N = %d,  excess = %d" % (r, N, N - 2))
    print("    |lam|<=   shapes   nonzero   values   largest fibre   mean fibre")
    print("    " + "-" * 62)
    fib = {}
    done = 0
    for hi in cutoffs:
        for size in range(done, hi + 1):
            for l in Partitions(size, max_length=N):
                v = psi(list(l), r, R, zs)
                if v != 0:
                    fib.setdefault(min(v, -v), []).append(tuple(l))
        done = hi + 1
        shapes = sum(Partitions(size, max_length=N).cardinality() for size in range(hi + 1))
        nz = sum(len(x) for x in fib.values())
        big = max(len(x) for x in fib.values())
        print("    %6d   %6d   %7d   %6d   %13d   %10.2f"
              % (hi, shapes, nz, len(fib), big, nz / len(fib)))

print("")
print("DONE")
