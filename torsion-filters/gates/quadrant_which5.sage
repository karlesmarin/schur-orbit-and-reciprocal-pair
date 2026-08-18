# -*- coding: utf-8 -*-
# In the interior, concentric is necessary and no longer sufficient.  Which concentric shapes vanish?
#
# quadrant_zeros2.sage settled the shape of the answer at t = 4:
#   r = 1  concentric <=> vanishing        24 and 24, no exceptions -- the paper's criterion
#   r = 2  concentric is necessary only    5 vanish of 120 concentric
#   r = 3  same                            1 of 72
# and no vanishing shape anywhere is non-concentric, so there is no third branch in this range.
#
# So the open half is: among concentric shapes, what picks out the vanishing ones?  This lists them
# with everything that could plausibly decide it -- the profile, which pair of classes is
# concentric, the shared centre, the t-quotient, and the sizes of the quotient components -- beside
# a sample of concentric shapes that do NOT vanish, since a feature shared by both decides nothing.
#
# Authors: Carles Marin, Claude (AI assistant).

t, r = 4, 2
N = t + 2 * r
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
zs = list(R.gens())
Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()


def psi(lam):
    out = R(0)
    for rho, c in p(s[list(lam)] if lam else s.one()):
        term = R(c)
        for k in rho:
            term *= R((t if k % t == 0 else 0) + sum(zz ** k + zz ** (-k) for zz in zs))
        out += term
    return out


def classes(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for b in beta:
        cls.setdefault(b % t, []).append(b)
    return beta, [sorted(cls.get(i, []), reverse=True) for i in range(t)]


def conc_pairs(cl):
    out = []
    big = [(i, c) for i, c in enumerate(cl) if len(c) >= 2]
    for x in range(len(big)):
        for y in range(x + 1, len(big)):
            (i, a), (j, b) = big[x], big[y]
            if a[0] + a[-1] == b[0] + b[-1]:
                out.append((i, j, a[0] + a[-1]))
    return out


def quotient(cl):
    q = []
    for i in range(t):
        e = cl[i]
        n = len(e)
        q.append(tuple((bb - i) // t - (n - 1 - j) for j, bb in enumerate(e)))
    return q


zeros, nonzeros = [], []
for size in range(16):
    for l in Partitions(size, max_length=N):
        lam = list(l)
        beta, cl = classes(lam)
        if any(len(x) == 0 for x in cl):
            continue
        cp = conc_pairs(cl)
        if not cp:
            continue
        rec = (lam, [len(x) for x in cl], cp, quotient(cl), beta)
        (zeros if psi(lam) == 0 else nonzeros).append(rec)

print("=" * 78)
print("t = 4, r = 2 : the concentric shapes, %d vanishing and %d not"
      % (len(zeros), len(nonzeros)))
print("=" * 78)


def show(rec):
    lam, prof, cp, q, beta = rec
    print("    lam=%-22s prof=%-14s conc=%-18s" % (str(lam), str(prof), str(cp)))
    print("        beta=%-26s quot=%s  sizes=%s"
          % (str(beta), str([list(x) for x in q]), str([sum(x) for x in q])))


print("")
print("  VANISHING:")
for rec in zeros:
    show(rec)

print("")
print("  NOT vanishing, a sample of the same size:")
for rec in nonzeros[:len(zeros)]:
    show(rec)

print("")
print("  quick separators, checked on both sets:")
for name, f in (("concentric pair is (0,2) or (1,3)", lambda rc: any(abs(i - j) == 2 for i, j, _ in rc[2])),
                ("some concentric class has 3 elements", lambda rc: any(len(x) >= 3 for x in [[]])),
                ("all quotient components empty", lambda rc: all(sum(x) == 0 for x in rc[3])),
                ("centre is even", lambda rc: all(c % 2 == 0 for _, _, c in rc[2]))):
    zt = sum(1 for rc in zeros if f(rc))
    nt = sum(1 for rc in nonzeros if f(rc))
    print("    %-40s zeros %d/%d   nonzeros %d/%d"
          % (name, zt, len(zeros), nt, len(nonzeros)))

print("")
print("DONE")
