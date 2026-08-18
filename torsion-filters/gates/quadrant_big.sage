# -*- coding: utf-8 -*-
# The same question with 30-40 zeros instead of 8.  Needs a different engine.
#
# Expanding s_lambda in power sums costs p(|lambda|) terms and dies well before the range where
# enough zeros live.  But the question is only WHETHER the value vanishes, and for that the
# bialternant is enough: Psi = det(x_j^{beta_i}) / Vandermonde, so Psi = 0 exactly when the
# numerator alternant vanishes identically.  Evaluate it at random points of a finite field
# containing the t-th roots of unity and it is one N x N determinant over F_p.
#
#   det != 0 at one point  =>  Psi != 0, certain.
#   det == 0 at several    =>  Psi == 0, overwhelmingly likely, and confirmed exactly afterwards
#                              on the survivors with the honest power-sum computation.
#
# Two controls, both mandatory:
#   the fast test must reproduce the exact sweep on |lambda| <= 17 shape for shape;
#   every zero it reports must survive the exact recomputation.
#
# Authors: Carles Marin, Claude (AI assistant).

t, r = 4, 2
N = t + 2 * r
MAX = 34
PRIME = 998244353            # = 119*2^23 + 1, so 1 mod 4 and the fourth roots of unity live here
TRIALS = 3

F = GF(PRIME)
I = F(-1).sqrt()
ROOTS = [F(1), I, F(-1), -I]
assert len(set(ROOTS)) == t and all(w ** t == 1 for w in ROOTS)

set_random_seed(12345)
POINTS = []
while len(POINTS) < TRIALS:
    z1, z2 = F.random_element(), F.random_element()
    if z1 == 0 or z2 == 0:
        continue
    alph = ROOTS + [z1, 1 / z1, z2, 1 / z2]
    if len(set(alph)) == N:
        POINTS.append(alph)


def vanishes_fast(beta):
    for alph in POINTS:
        M = matrix(F, N, N, lambda a, b: alph[b] ** beta[a])
        if M.det() != 0:
            return False
    return True


Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
zs = list(R.gens())


def vanishes_exact(lam):
    out = R(0)
    for rho, c in p(s[list(lam)] if lam else s.one()):
        term = R(c)
        for k in rho:
            term *= R((t if k % t == 0 else 0) + sum(zz ** k + zz ** (-k) for zz in zs))
        out += term
    return out == 0


def info(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = [sorted([b for b in beta if b % t == i], reverse=True) for i in range(t)]
    return beta, cls


def conc_pairs(cls):
    out = []
    big = [(i, c) for i, c in enumerate(cls) if len(c) >= 2]
    for x in range(len(big)):
        for y in range(x + 1, len(big)):
            (i, a), (j, b) = big[x], big[y]
            if a[0] + a[-1] == b[0] + b[-1]:
                out.append((i, j, a[0] + a[-1]))
    return out


def quot(cls):
    return [tuple((b - i) // t - (len(cls[i]) - 1 - j) for j, b in enumerate(cls[i]))
            for i in range(t)]


print("=" * 78)
print("t = 4, r = 2, |lambda| <= %d, zeros by finite-field bialternant" % MAX)
print("=" * 78)

zeros, others, nconc = [], [], 0
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        lam = list(l)
        beta, cls = info(lam)
        if any(len(c) == 0 for c in cls):
            continue
        cp = conc_pairs(cls)
        if not cp:
            continue                      # concentric is necessary, established earlier
        nconc += 1
        rec = (lam, beta, cls, cp)
        (zeros if vanishes_fast(beta) else others).append(rec)

print("")
print("  concentric shapes with no empty class : %d" % nconc)
print("  of them vanishing                     : %d" % len(zeros))

bad = [rec for rec in zeros if not vanishes_exact(rec[0])]
print("  CONTROL, exact recomputation of every zero: %d disagreements" % len(bad))
small = [rec[0] for rec in zeros if sum(rec[0]) <= 17]
print("  CONTROL, the 8 known zeros up to |lambda|<=17 are all here: %s (found %d)"
      % ("yes" if len(small) == 8 else "NO", len(small)))


def f_beta_sym(rec):
    return any(set(C - b for b in rec[1]) == set(rec[1]) for _, _, C in rec[3])


def f_one_empty(rec):
    q = quot(rec[2])
    return any(sum(q[i]) == 0 or sum(q[j]) == 0 for i, j, _ in rec[3])


def f_both_empty(rec):
    q = quot(rec[2])
    return any(sum(q[i]) == 0 and sum(q[j]) == 0 for i, j, _ in rec[3])


def f_equal_size(rec):
    return any(len(rec[2][i]) == len(rec[2][j]) for i, j, _ in rec[3])


def f_two_pairs(rec):
    return len(rec[3]) >= 2


def f_centre_max(rec):
    return any(C == max(rec[1]) for _, _, C in rec[3])


print("")
print("  %-42s %10s %12s" % ("condition", "vanishing", "not"))
print("  " + "-" * 68)
for name, f in (("beta symmetric about the centre", f_beta_sym),
                ("ONE of the pair has empty quotient", f_one_empty),
                ("BOTH of the pair have empty quotient", f_both_empty),
                ("concentric classes of equal size", f_equal_size),
                ("two or more concentric pairs", f_two_pairs),
                ("the centre equals max(beta)", f_centre_max)):
    a = sum(1 for x in zeros if f(x))
    b = sum(1 for x in others if f(x))
    tag = ""
    if a == len(zeros) and b == 0:
        tag = "  <-- SEPARATES"
    elif b == 0:
        tag = "  (sufficient)"
    elif a == len(zeros):
        tag = "  (necessary)"
    print("  %-42s %5d/%-5d %6d/%-6d%s" % (name, a, len(zeros), b, len(others), tag))

print("")
print("  every vanishing shape:")
for lam, beta, cls, cp in zeros:
    print("    |lam|=%-3d lam=%-30s beta=%-32s conc=%s  betasym=%s"
          % (sum(lam), str(lam), str(beta), str([(i, j) for i, j, _ in cp]),
             "yes" if f_beta_sym((lam, beta, cls, cp)) else "."))

print("")
print("DONE")
