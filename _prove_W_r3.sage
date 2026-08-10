# Proving (W): the isolating witness, and how far it reaches.
#
# C_mu(lam) = c_mu(lam) +- c_{mu*}(lam) + sum over NON-STANDARD nu reducing to mu of +- c_nu(lam).
#
# Two exact facts, proved:
#   (i)  |mu*| = |mu| + N - 2 ell(mu), and ell(mu) <= r with N = 2r+2 gives |mu*| >= |mu| + 2,
#        and ell(mu*) = N - ell(mu).
#   (ii) every non-standard nu has nu'_1 + nu'_2 > N with nu'_2 <= nu'_1, hence 2 ell(nu) > N,
#        i.e. ell(nu) >= r + 2.  So in the stable range NO non-standard nu fits inside lam.
#
# ISOLATION LEMMA (proved): if ell(mu) <= N - 1 - ell(lam) then ell(mu*) = N - ell(mu) > ell(lam),
# so mu* is not contained in lam and c_{mu*} = 0.
#
# What is left to bound is the non-standard part. This script measures:
#   (A) the true minimum of ell(nu) over NON-STANDARD nu with o_nu(A) != 0  -- if it exceeds
#       r+2 the isolation lemma reaches further than (ii) alone allows;
#   (B) for every unstable lam, whether an ISOLATING mu exists (all competing nu fail to fit in
#       lam) with c_mu != 0 -- i.e. a witness that is PROVED, not merely observed;
#   (C) the residue: unstable non-vanishing lam with no isolating mu.
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

def c_univ(lam, nu):
    lam = Partition(list(lam)); nu = Partition(list(nu))
    d = lam.size() - nu.size()
    if d < 0 or d % 2 != 0:
        return 0
    tot = 0
    for hf in Partitions(d // 2):
        tot += (s(nu) * s(Partition([2 * x for x in hf]))).coefficient(lam)
    return tot

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
    red, target = {}, {}
    for k in vals:
        if vals[k] == 0:
            red[k] = None; continue
        try:
            sol = Bm.solve_right(vector(QQ, vecof(vals[k], mons)))
            red[k] = {basis[i]: sol[i] for i in range(len(basis)) if sol[i] != 0}
        except ValueError:
            red[k] = None
    # (A) minimum ell over non-standard nu with o != 0, and which mu they feed
    ns = [k for k in vals if conj12(k)[0] + conj12(k)[1] > N and vals[k] != 0]
    minell = min((len([x for x in k if x > 0]) for k in ns), default=None)
    feeds = {}
    for k in ns:
        for mu in (red[k] or {}):
            feeds.setdefault(mu, []).append(len([x for x in k if x > 0]))
    print("=" * 76)
    print("r=%d N=%d : non-standard nu with o!=0 : %d ; min ell(nu) = %s  (fact (ii) gives >= %d)"
          % (r, N, len(ns), str(minell), r + 2))

    # (B),(C)
    proved = residue = nvan = 0
    reslist = []
    for tot in range(0, LMAX + 1):
        for lam in Partitions(tot, max_length=N):
            ell = len([x for x in lam if x > 0])
            if ell <= half:
                continue
            if schur_at(lam, R, A) == 0:
                nvan += 1
                continue
            # an ISOLATING mu: every competing nu must fail to fit inside lam
            found = None
            for mu in sorted(basis, key=lambda m: (-sum(m), m)):
                # ALL nu that feed mu, mu itself included: exactly one may be nonzero
                contribs = [tuple(mu)]
                cs = list(Partition(list(mu)).conjugate())
                first = N - (cs[0] if cs else 0)
                rest = cs[1:]
                if first >= 0 and first >= (rest[0] if rest else 0):
                    contribs.append(tuple(Partition([first] + rest).conjugate()))
                for k in ns:
                    if mu in (red[k] or {}):
                        contribs.append(k)
                nzc = [c for c in contribs if c_univ(lam, c) != 0]
                if len(nzc) == 1:
                    found = (mu, nzc[0])
                    break
            if found:
                proved += 1
            else:
                residue += 1
                if len(reslist) < 10:
                    reslist.append((tuple(lam), ell))
    print("   UNSTABLE shapes: %d vanish ; %d have a PROVED isolating witness ; %d residue"
          % (nvan, proved, residue))
    if reslist:
        print("   residue examples:", reslist)
    return residue

t = 0
t += run(3, 13)
print("\nTOTAL unstable residue without a proved isolating witness: %d" % t)
