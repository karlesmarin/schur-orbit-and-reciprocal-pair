# -*- coding: utf-8 -*-
# What do the members of a large fibre at excess 4 share?
#
# The growth measurement (excess_fibre_growth.sage) killed the hypothesis that the collapse is
# special to excess 2: at excess 4 the largest fibre also grows without saturating, 3 5 7 7 9 11
# over cutoffs 6..16.  So there IS an invariant there; what is special to excess 2 is only that the
# invariant's value happens to be a PRODUCT of three characters.
#
# So look at the fibres and read off what their members have in common, instead of guessing a datum
# and testing it.  For the largest fibres this prints, per member: the beta set, the profile mod 2,
# the within-class gap tuples, the cross offset, the sign of sigma, and the 2-quotient.  Whatever is
# constant down a column is a candidate; whatever varies is not in the invariant.
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


def anatomy(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {0: [], 1: []}
    for b in beta:
        cls[b % 2].append(b)
    gaps, quot, sums = [], [], []
    for i in (0, 1):
        e = sorted(cls[i], reverse=True)
        gaps.append(tuple(e[j] - e[j + 1] for j in range(len(e) - 1)))
        n = len(e)
        quot.append(tuple((e[j] - i) // 2 - (n - 1 - j) for j in range(n)))
        sums.append(sum(e))
    order = sorted(range(N), key=lambda j: (beta[j] % 2, -beta[j]))
    sg = 1
    for a in range(N):
        for b in range(a + 1, N):
            if order[a] > order[b]:
                sg = -sg
    return dict(beta=tuple(beta), prof=(len(cls[0]), len(cls[1])),
                gaps=tuple(gaps), cross=sums[0] - sums[1], sgn=sg,
                quot=tuple(quot), size=sum(lam))


r, N, MAX = 2, 6, 16
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
zs = list(R.gens())

fib = {}
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        v = psi(list(l), r, R, zs)
        if v != 0:
            fib.setdefault(min(v, -v), []).append(tuple(l))

order = sorted(fib.items(), key=lambda kv: -len(kv[1]))
print("=" * 78)
print("Anatomy of the largest fibres at excess 4  (r = 2, N = 6, |lambda| <= %d)" % MAX)
print("=" * 78)

for val, members in order[:3]:
    print("")
    print("  fibre of %d shapes" % len(members))
    rows = [anatomy(list(m), N) for m in members]
    for f in ("prof", "gaps", "cross", "sgn"):
        vals = set(str(x[f]) for x in rows)
        mark = "CONSTANT" if len(vals) == 1 else "varies (%d)" % len(vals)
        print("    %-6s %-11s %s" % (f, mark, sorted(vals)[:4]))
    print("    %-22s %s" % ("sizes |lambda|", sorted(set(x['size'] for x in rows))))
    for m, f in list(zip(members, rows))[:6]:
        print("      lam=%-18s beta=%-22s quot=%s"
              % (str(list(m)), str(list(f['beta'])), str(list(f['quot']))))
    if len(members) > 6:
        print("      ... and %d more" % (len(members) - 6))

print("")
print("DONE")
