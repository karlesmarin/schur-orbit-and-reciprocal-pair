# -*- coding: utf-8 -*-
# Excess 4: expand Psi_2 in symplectic characters and look at the SUPPORT.
#
# The determinant question has an answer that is not a conjecture.  Littlewood's restriction gives
# s_lambda = sum_mu (sum_beta c^lambda_{mu,beta}) o_mu, and Lemma 8.8 of the paper gives
# o_nu(A) = sp_nu(z_1,...,z_r) on this alphabet.  Composing,
#
#     Psi_r = sum_mu c_mu sp_mu(z_1,...,z_r),      c_mu in Z.
#
# So the value IS a combination of symplectic characters; what is unknown is which mu occur and
# whether the support can be read off the 2-quotient, the way the three arguments are at excess 2.
# At excess 2 the answer is one term after the sinh denominator is cleared, which is why it looks
# like a product; if the support stays small at excess 4, the graph programme has its object there.
#
# Method: peel greedily.  sp_mu at (z_1,z_2) is Weyl's quotient of alternants, and its leading
# monomial in dominance order is z_1^{mu_1} z_2^{mu_2} with coefficient 1, so the highest remaining
# monomial names the next mu.  Two controls, and both must hold or nothing here is worth reading:
#
#   the remainder after peeling must be EXACTLY zero -- otherwise the expansion is not in this span;
#   the coefficients must be INTEGERS -- otherwise the peeling picked the wrong basis.
#
# Authors: Carles Marin, Claude (AI assistant).

r, N = 2, 6
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
z0, z1 = R.gens()
zs = [z0, z1]
Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()


def psi(lam):
    out = R(0)
    for rho, c in p(s[list(lam)] if lam else s.one()):
        term = R(c)
        for k in rho:
            term *= R(1 + (-1) ** k + sum(zz ** k + zz ** (-k) for zz in zs))
        out += term
    return out


def alt(exps):
    """det[ z_j^{e_i} - z_j^{-e_i} ]."""
    M = matrix(R, r, r, lambda i, j: zs[j] ** exps[i] - zs[j] ** (-exps[i]))
    return M.det()


DEN = alt([r - i for i in range(r)])              # exponents r, r-1, ..., 1  -> here 2, 1


def sp(mu):
    """Symplectic character of Sp(2r) at (z_1,...,z_r), l(mu) <= r."""
    mu = list(mu) + [0] * (r - len(mu))
    num = alt([mu[i] + r - i for i in range(r)])
    # exact division, not quo_rem: in a multivariate LAURENT ring quo_rem is not the division we
    # want and reports a nonzero remainder even for sp of the empty partition, where numerator and
    # denominator are literally the same alternant.
    q = R(num / DEN)
    assert q * DEN == num, "sp(%s) did not divide" % (mu,)
    return q


def terms(f):
    """{exponent tuple: coefficient}.  Rebuilt with plain tuples as keys: the ring hands back
    ETuple objects, which do not hash equal to the tuples I build here, and a lookup by tuple
    raises KeyError on a monomial that is plainly present."""
    return dict((tuple(e), c) for e, c in f.dict().items())


def lead(f):
    """The dominance-highest monomial exponent of a Laurent polynomial in z0, z1."""
    best = None
    for e in terms(f):
        if best is None or (sum(e), e) > (sum(best), best):
            best = e
    return best


def peel(f):
    """Write f = sum c_mu sp_mu; returns the list of (mu, c) or None if it does not close."""
    out = []
    g = R(f)
    for _ in range(200):
        if g == 0:
            return out
        e = lead(g)
        if e[0] < e[1] or e[1] < 0:
            return None                      # leading term is not a dominant weight
        mu = (e[0], e[1])
        c = terms(g)[e]
        if c not in ZZ:
            return None
        out.append((mu, c))
        g = g - c * sp(mu)
    return None


def quotient2(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    out = []
    for i in (0, 1):
        e = sorted([b for b in beta if b % 2 == i], reverse=True)
        n = len(e)
        out.append(tuple((e[j] - i) // 2 - (n - 1 - j) for j in range(n)))
    return tuple(out)


print("=" * 78)
print("Psi_2 in the symplectic basis: how big is the support?")
print("=" * 78)

hist = {}
failed = 0
rows = []
for size in range(0, 13):
    for l in Partitions(size, max_length=N):
        v = psi(list(l))
        if v == 0:
            continue
        ex = peel(v)
        if ex is None:
            failed += 1
            continue
        hist[len(ex)] = hist.get(len(ex), 0) + 1
        rows.append((list(l), ex, quotient2(list(l))))

print("")
print("  CONTROL: expansions that failed to close over Z in this span: %d of %d"
      % (failed, failed + len(rows)))
print("")
print("  terms in the expansion   how many shapes")
for k in sorted(hist):
    print("    %-22d %d" % (k, hist[k]))

print("")
print("  the shortest expansions, with the 2-quotient beside them:")
rows.sort(key=lambda t: (len(t[1]), sum(t[0])))
for lam, ex, q in rows[:10]:
    print("    lam=%-16s  quot=%-26s  %s"
          % (str(lam), str([list(x) for x in q]),
             " + ".join("%d*sp%s" % (c, list(m)) for m, c in ex)))

print("")
print("  and the longest:")
for lam, ex, q in rows[-4:]:
    print("    lam=%-16s  %d terms: %s"
          % (str(lam), len(ex),
             " + ".join("%d*sp%s" % (c, list(m)) for m, c in ex[:5]) + (" + ..." if len(ex) > 5 else "")))

print("")
print("DONE")
