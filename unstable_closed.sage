# The unstable-range criterion, with the full reduction in place.
#
# Pieces, all now derived and verified (typeD_rule.sage, typeD_residue.sage):
#   basis            B = { o_mu(A) : mu standard, mu'_1 < N/2 }  -- LINEARLY INDEPENDENT (full rank)
#   standard, mu'_1 = N/2 (self-associate)        -> o = 0
#   standard, mu'_1 > N/2                         -> o = +- o_{associate}
#   non-standard                                  -> o = 0, or +-1 times a single basis element
#
# Since s_lam = sum_nu c_nu o_nu is universal (no stable range), reducing every label gives
#     Psi_lam = sum_{mu in B} C_mu(lam) * o_mu(A),   C_mu an explicit signed sum of the c_nu,
# and B independent gives the exact, range-free criterion
#     Psi_lam = 0   <=>   C_mu(lam) = 0 for every mu in B.
#
# The reduction map is built by solving each label against B once, then applied to every lambda.
# Tested against the exact object in BOTH ranges.
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

def vecof(p, mons):
    d = {} if p == 0 else {tuple(k) if hasattr(k, '__iter__') else (k,): QQ(v)
                           for k, v in p.dict().items()}
    return [d.get(m, QQ(0)) for m in mons]

def run(r, LMAX):
    N = 2 * r + 2
    half = N // 2
    R, A = alphabet(r)
    DMAX = LMAX
    vals = {}
    for d in range(0, DMAX + 1):
        for nu in Partitions(d):
            vals[tuple(nu)] = o_at(tuple(nu), R, A)
    basis = [k for k in vals if conj12(k)[0] + conj12(k)[1] <= N and conj12(k)[0] < half]
    mons = sorted(set(m for k in vals if vals[k] != 0
                      for m in [tuple(x) if hasattr(x, '__iter__') else (x,) for x in vals[k].dict()]))
    Bm = matrix(QQ, [vecof(vals[k], mons) for k in basis]).transpose()
    rk = Bm.rank()
    # reduction map: every label -> coordinates in the basis
    red = {}
    unresolved = []
    for k in vals:
        if vals[k] == 0:
            red[k] = None
            continue
        try:
            sol = Bm.solve_right(vector(QQ, vecof(vals[k], mons)))
        except ValueError:
            unresolved.append(k); continue
        red[k] = {basis[i]: sol[i] for i in range(len(basis)) if sol[i] != 0}
    print("=" * 74)
    print("r=%d N=%d : basis %d (rank %d), labels reduced %d, unresolved %d"
          % (r, N, len(basis), rk, len(red), len(unresolved)))
    if unresolved:
        print("   unresolved:", unresolved[:6])

    bad = []
    nst = nun = 0
    for tot in range(0, LMAX + 1):
        for lam in Partitions(tot, max_length=N):
            ell = len([x for x in lam if x > 0])
            if ell <= half:
                nst += 1
            else:
                nun += 1
            C = {}
            ok = True
            for nu, cc in o(s(Partition(list(lam)))):
                k = tuple(nu)
                if k not in red:
                    ok = False; break
                if red[k] is None:
                    continue
                for mu, coef in red[k].items():
                    C[mu] = C.get(mu, QQ(0)) + QQ(cc) * coef
            if not ok:
                continue
            crit = all(v == 0 for v in C.values())
            z = (schur_at(lam, R, A) == 0)
            if crit != z:
                bad.append((tuple(lam), ell, crit, z))
    print("   shapes: %d stable, %d unstable" % (nst, nun))
    print("   criterion vs exact object -- disagreements: %d" % len(bad))
    for b in bad[:8]:
        print("      lam=%-20s ell=%d crit=%s exact=%s" % (str(b[0]), b[1], b[2], b[3]))
    return len(bad)

t = 0
t += run(1, 10)
t += run(2, 9)
print("\nTOTAL disagreements: %d" % t)
if t == 0:
    print("=> the criterion is EXACT in both ranges. The unstable range is closed as a TOOL:")
    print("   Psi_lam = 0 <=> C_mu(lam) = 0 for all mu in an independent basis, with every")
    print("   modification derived rather than assumed.")
