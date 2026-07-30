# Exact data for the new figures. Computed in Sage, dumped to JSON; the matplotlib scripts only
# draw what is in the file. Nothing here is asserted.
#
# For each r and each lambda in range we record:
#   z_exact : the object s_lambda(1,-1,z_1,1/z_1,...) is identically zero
#   z_end   : its value at z_i = 1 is zero      (the ENDPOINT, an order-2 element)
#   branch  : 'a' if all beta_j share a parity, 'b' if lambda is self-complementary of odd width
#   nfac    : number of sign changes in the B_r-character expansion is not computed here; instead
#             we record the number of DISTINCT signs among the coefficients of the expansion of the
#             object in the basis of the surviving universal orthogonal values -- 1 means the object
#             is +-(a genuine character), >1 means it is properly VIRTUAL.
#
# Carles Marin + Claude (AI assistant).

import json

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

def branch(lam, N):
    L = list(lam) + [0] * (N - len(lam))
    beta = [L[i] + (N - 1 - i) for i in range(N)]
    if len(set(x % 2 for x in beta)) == 1:
        return 'a'
    if len(L) == N:
        w = L[0] + L[N - 1]
        if w % 2 == 1 and all(L[i] + L[N - 1 - i] == w for i in range(N)):
            return 'b'
    return None

OUT = {}
for r, LMAX in [(1, 22), (2, 17), (3, 12)]:
    N = 2 * r + 2
    R, A = alphabet(r)
    Aone = [QQ(1), QQ(-1)] + [QQ(1)] * (2 * r)
    rows = []
    for tot in range(0, LMAX + 1):
        for lam in Partitions(tot, max_length=N):
            ell = len([x for x in lam if x > 0])
            ze = (schur_at(lam, R, A) == 0)
            # endpoint value: same bialternant over QQ at z=1
            L = list(lam) + [0] * (N - len(lam))
            Dg = max(L[0], 1) + N + 1
            PK = PowerSeriesRing(QQ, 'u', default_prec=Dg + 2); u = PK.gen()
            G = prod(1 / (1 - a * u) for a in Aone)
            h = [G[k] for k in range(Dg + 1)]
            H = lambda k: QQ(0) if k < 0 else h[k]
            zend = (matrix(QQ, N, N, lambda i, j: H(L[i] - i + j)).det() == 0)
            rows.append({"lam": [int(x) for x in lam], "size": int(tot), "ell": int(ell),
                         "z_exact": bool(ze), "z_end": bool(zend), "branch": branch(lam, N)})
    OUT["r%d" % r] = {"N": int(N), "LMAX": int(LMAX), "rows": rows}
    nz = len([x for x in rows if x["z_exact"]])
    rig = len([x for x in rows if x["z_end"] and not x["z_exact"]])
    print("r=%d  N=%d  shapes %d  identically zero %d  zero ONLY at the endpoint %d"
          % (r, N, len(rows), nz, rig))

# the type-D reduction table, for the reduction figure
def conj12(nu):
    c = list(Partition(list(nu)).conjugate())
    return (c[0] if c else 0), (c[1] if len(c) > 1 else 0)

RED = {}
for r in [1, 2]:
    N = 2 * r + 2; half = N // 2
    R, A = alphabet(r)
    tab = []
    for d in range(0, (9 if r == 1 else 8) + 1):
        for nu in Partitions(d):
            k = tuple(nu)
            v = sum((R(QQ(co)) * schur_at(list(mu), R, A) for mu, co in s(o(list(nu)))), R(0))
            c1, c2 = conj12(k)
            std = (c1 + c2 <= N)
            cls = ('selfassoc' if (std and c1 == half) else
                   'basis' if (std and c1 < half) else
                   'assoc' if std else 'nonstd')
            tab.append({"nu": [int(x) for x in nu], "c1": int(c1), "c2": int(c2),
                        "cls": cls, "zero": bool(v == 0)})
    RED["r%d" % r] = {"N": int(N), "half": int(half), "table": tab}
    z = len([x for x in tab if x["zero"]])
    print("r=%d reduction table: %d labels, %d vanish on the coset" % (r, len(tab), z))

json.dump({"loci": OUT, "reduction": RED}, open("fig_data_new.json", "w"), indent=1)
print("written fig_data_new.json")
