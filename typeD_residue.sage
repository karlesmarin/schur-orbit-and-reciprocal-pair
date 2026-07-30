# The residue class reduced by exact linear algebra, not by hunting King 1971.
#
# FJKW (math-ph/0508034, source read): "Classically one has only two cases dealt with by
# case-by-case studies [King 1971]. The general theory, having infinitely many cases, needs a not
# yet available formalism." So there is no general statement to find -- but their METHOD is to
# express non-standard characters in terms of standard ones using the determinant character, case
# by case. That is exactly what can be done here, exactly, by linear algebra.
#
# From typeD_rule.sage, on our coset alphabet A, the exact table is:
#   standard, self-associate (nu'_1 = N/2)  -> o_nu(A) = 0     (41/41)
#   standard, otherwise                     -> o_nu(A) != 0    (59/59)
#   standard with nu'_1 > N/2               -> o_nu = +- o_{nu*}  (23/23)
#   NON-standard (nu'_1 + nu'_2 > N)        -> 85 of 106 vanish; 21 do NOT   <-- THE RESIDUE
#
# The residue is what broke (C-univ). Here each residue value is expressed in the basis
#   B = { o_mu(A) : mu standard, mu'_1 < N/2 }
# by solving a linear system over Q. If the coefficients come out as small integers, that IS the
# missing rule, in the only form we need it.
#
# Carles Marin + Claude (AI assistant).

Sym = SymmetricFunctions(QQ); s = Sym.s(); o = Sym.o()

def alphabet(r):
    R = LaurentPolynomialRing(QQ, ['z%d' % i for i in range(1, r + 1)])
    zs = R.gens()
    A = [R(1), R(-1)]
    for z in zs:
        A += [z, z**-1]
    return R, A

def schur_at(lam, R, A):
    N = len(A)
    L = list(lam) + [0] * max(0, N - len(lam))
    if len([x for x in lam if x > 0]) > N:
        return R(0)
    L = L[:N]
    D = max(L[0], 1) + N + 1
    PK = PowerSeriesRing(R, 'u', default_prec=D + 2); u = PK.gen()
    G = prod(1 / (1 - a * u) for a in A)
    h = [G[k] for k in range(D + 1)]
    H = lambda k: R(0) if k < 0 else h[k]
    return matrix(R, N, N, lambda i, j: H(L[i] - i + j)).det()

def o_at(nu, R, A):
    tot = R(0)
    for mu, co in s(o(list(nu))):
        tot += R(QQ(co)) * schur_at(list(mu), R, A)
    return tot

def conj12(nu):
    c = list(Partition(list(nu)).conjugate())
    return (c[0] if c else 0), (c[1] if len(c) > 1 else 0)

def vec(p, R, mons):
    d = {} if p == 0 else {tuple(k) if hasattr(k, '__iter__') else (k,): QQ(v)
                           for k, v in p.dict().items()}
    return [d.get(m, QQ(0)) for m in mons]

for r in [1, 2]:
    N = 2 * r + 2
    half = N // 2
    R, A = alphabet(r)
    DMAX = 10 if r == 1 else 8
    print("=" * 76)
    print("r=%d  N=%d  half=%d" % (r, N, half))

    basis_nu, residue_nu, vals = [], [], {}
    for d in range(0, DMAX + 1):
        for nu in Partitions(d):
            k = tuple(nu)
            v = o_at(k, R, A)
            vals[k] = v
            c1, c2 = conj12(k)
            if c1 + c2 <= N and c1 < half:
                basis_nu.append(k)
            elif c1 + c2 > N and v != 0:
                residue_nu.append(k)
    print("  basis (standard, nu'_1 < half): %d ;  residue (non-standard, nonzero): %d"
          % (len(basis_nu), len(residue_nu)))

    mons = sorted(set(m for k in vals for m in
                      ([] if vals[k] == 0 else
                       [tuple(x) if hasattr(x, '__iter__') else (x,) for x in vals[k].dict()])))
    Bm = matrix(QQ, [vec(vals[k], R, mons) for k in basis_nu]).transpose()
    print("  basis rank: %d of %d columns" % (Bm.rank(), len(basis_nu)))

    print("  reductions of the residue, solved exactly:")
    intcoef = allsolved = 0
    for k in residue_nu:
        try:
            sol = Bm.solve_right(vector(QQ, vec(vals[k], R, mons)))
        except ValueError:
            print("     nu=%-22s NOT in the span of the standard values" % str(k))
            continue
        allsolved += 1
        terms = [(basis_nu[i], sol[i]) for i in range(len(basis_nu)) if sol[i] != 0]
        if all(c.denominator() == 1 for _, c in terms):
            intcoef += 1
        c1, c2 = conj12(k)
        pretty = " + ".join("%s*o%s" % (c, str(mu)) for mu, c in terms) if terms else "0"
        print("     nu=%-20s (c1,c2)=(%d,%d)  ->  %s" % (str(k), c1, c2, pretty[:96]))
    print("  residue values inside the standard span: %d/%d ; with INTEGER coefficients: %d"
          % (allsolved, len(residue_nu), intcoef))
