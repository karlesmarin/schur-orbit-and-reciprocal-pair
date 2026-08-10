# The one shape at height r+2 with no isolating mu -- is it real, or an edge of the range?
#
# iso_first_row.sage found exactly one lambda at ell = r+2 with Psi_r != 0 and no isolating mu:
# lambda = (5,5,2,1) at r = 2.  Its size 13 sat at the top of that run's range, where associates
# and feeders of size > 13 are simply absent from the tables -- an absence that makes isolation
# EASIER to find, not harder, so it cannot manufacture a failure.  Still, the shape is re-examined
# here with room above it, and the competitors are printed rather than counted.
#
# Carles Marin + Claude (AI assistant).

Sym = SymmetricFunctions(QQ); s = Sym.s(); o = Sym.o()
LMAX = 17
r = 2
N = 2 * r + 2
half = N // 2
TARGETS = [(5, 5, 2, 1)]

def alphabet(r):
    R = LaurentPolynomialRing(QQ, ['z%d' % i for i in range(1, r + 1)])
    A = [R(1), R(-1)]
    for z in R.gens():
        A += [z, z**-1]
    return R, A

def schur_at(lam, R, A):
    n = len(A)
    L = list(lam) + [0] * max(0, n - len(lam))
    if len([x for x in lam if x > 0]) > n:
        return R(0)
    L = L[:n]
    D = max(L[0], 1) + n + 1
    PK = PowerSeriesRing(R, 'u', default_prec=D + 2); u = PK.gen()
    G = prod(1 / (1 - a * u) for a in A)
    h = [G[k] for k in range(D + 1)]
    H = lambda k: R(0) if k < 0 else h[k]
    return matrix(R, n, n, lambda i, j: H(L[i] - i + j)).det()

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
    key = (tuple(lam), tuple(nu))
    if key in _CU:
        return _CU[key]
    la = Partition(list(lam)); nu2 = Partition(list(nu))
    d = la.size() - nu2.size()
    if d < 0 or d % 2 != 0:
        _CU[key] = 0; return 0
    tot = 0
    for hf in Partitions(d // 2):
        tot += (s(nu2) * s(Partition([2 * x for x in hf]))).coefficient(la)
    _CU[key] = tot
    return tot

def contains(lam, nu):
    L = list(lam) + [0] * 24
    M = list(nu) + [0] * 24
    return all(M[i] <= L[i] for i in range(24))

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
feeders = {}
for k in vals:
    for mu in (red[k] or {}):
        feeders.setdefault(mu, set()).add(k)

print("r = %d, N = %d, tables built to |nu| <= %d" % (r, N, LMAX))
for lam in TARGETS:
    print("")
    print("=" * 72)
    print("lambda = %s   |lambda| = %d   ell = %d   Psi_r != 0 : %s"
          % (str(lam), sum(lam), len([x for x in lam if x > 0]), schur_at(lam, R, A) != 0))
    print("=" * 72)
    best = None
    for mu in sorted(basis, key=lambda m: (sum(m), m)):
        fits = sorted(k for k in feeders.get(mu, set())
                      if contains(lam, k) and c_univ(lam, k) != 0)
        if not fits:
            continue
        mark = "  <-- ISOLATING" if len(fits) == 1 else ""
        print("  mu=%-12s competitors inside lambda: %d  %s%s"
              % (str(mu), len(fits), [tuple(f) for f in fits][:4], mark))
        if best is None or len(fits) < best[1]:
            best = (mu, len(fits))
    print("")
    print("  fewest competitors over all mu: %s" % (str(best)))
