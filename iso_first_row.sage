# Attacking the converse of Section 8 at the FIRST unstable height, ell(lambda) = r+2.
#
# Conjecture (W) asks for an isolating mu for every lambda outside branches (a) and (b).  In the
# stable range ell(lambda) <= N/2 the paper builds one explicitly.  The first height beyond it is
# ell(lambda) = r+2 = N/2 + 1, and there the competitors are pinched from both sides:
#
#   * the associate needs ell(mu*) = N - ell(mu) > ell(lambda) = r+2, i.e. ell(mu) <= r-1;
#   * a NON-STANDARD nu needs ell(nu) >= r+2 by fact (ii), so to fit inside lambda it must have
#     ell(nu) = r+2 EXACTLY and nu subset lambda.
#
# So at this one height the question is finite and tight.  This script asks, for every such lambda
# with Psi_r != 0: which mu isolate, and is the truncation mu = (lambda_1, ..., lambda_{r-1})
# -- the tallest mu the associate bound allows -- always among them?
#
# Carles Marin + Claude (AI assistant).

Sym = SymmetricFunctions(QQ); s = Sym.s(); o = Sym.o()

def alphabet(r):
    R = LaurentPolynomialRing(QQ, ['z%d' % i for i in range(1, r + 1)])
    A = [R(1), R(-1)]
    for z in R.gens():
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

_CU = {}
def c_univ(lam, nu):
    """c_nu(lam) = sum over even beta' of the Littlewood-Richardson coefficient."""
    key = (tuple(lam), tuple(nu))
    if key in _CU:
        return _CU[key]
    lam = Partition(list(lam)); nu = Partition(list(nu))
    d = lam.size() - nu.size()
    if d < 0 or d % 2 != 0:
        return 0
    tot = 0
    for hf in Partitions(d // 2):
        tot += (s(nu) * s(Partition([2 * x for x in hf]))).coefficient(lam)
    _CU[key] = tot
    return tot

def contains(lam, nu):
    L = list(lam) + [0] * 20
    M = list(nu) + [0] * 20
    return all(M[i] <= L[i] for i in range(20))

def run(r, LMAX):
    N = 2 * r + 2
    half = N // 2
    R, A = alphabet(r)
    vals = {}
    for d in range(0, LMAX + 1):
        for nu in Partitions(d):
            vals[tuple(nu)] = o_at(tuple(nu), R, A)
    basis = [k for k in vals if conj12(k)[0] + conj12(k)[1] <= N and conj12(k)[0] < half]
    mons = sorted(set(m for k in vals if vals[k] != 0
                      for m in [tuple(x) if hasattr(x, '__iter__') else (x,) for x in vals[k].dict()]))
    Bm = matrix(QQ, [vecof(vals[k], mons) for k in basis]).transpose()
    red = {}
    for k in vals:
        if vals[k] == 0:
            red[k] = None; continue
        try:
            sol = Bm.solve_right(vector(QQ, vecof(vals[k], mons)))
            red[k] = {basis[i]: sol[i] for i in range(len(basis)) if sol[i] != 0}
        except ValueError:
            red[k] = None
    # which labels feed each mu
    feeders = {}
    for k in vals:
        for mu in (red[k] or {}):
            feeders.setdefault(mu, set()).add(k)

    tested = trunc_ok = any_ok = none_ok = 0
    fails = []
    for tot in range(0, LMAX + 1):
        for lam in Partitions(tot, max_length=N):
            ell = len([x for x in lam if x > 0])
            if ell != half + 1:                       # the FIRST unstable height
                continue
            if schur_at(lam, R, A) == 0:
                continue
            tested += 1
            iso = []
            for mu in basis:
                # NOT "c_mu != 0": the surviving feeder may be any label reducing to mu, and for
                # tall thin lambda it always is -- c_mu itself vanishes there.
                fits = [k for k in feeders.get(mu, set()) if contains(lam, k) and c_univ(lam, k) != 0]
                if len(fits) == 1:
                    iso.append(mu)
            if iso:
                any_ok += 1
            else:
                none_ok += 1
                fails.append(tuple(lam))
            cand = tuple([x for x in list(lam)[:max(0, r - 1)] if x > 0])
            if cand in [tuple(m) for m in iso]:
                trunc_ok += 1
    print("=" * 74)
    print("r = %d,  N = %d,  ell(lambda) = %d,  |lambda| <= %d" % (r, N, half + 1, LMAX))
    print("  non-vanishing shapes at that height : %d" % tested)
    print("  some isolating mu exists            : %d" % any_ok)
    print("  none                                : %d   %s" % (none_ok, fails[:6]))
    print("  the truncation (lam_1..lam_{r-1}) works : %d" % trunc_ok)

run(2, 13)
