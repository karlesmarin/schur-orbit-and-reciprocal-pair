# The type-D specialization rule, DERIVED from our own alphabet instead of hunted down.
#
# Koike-Terada Prop 2.4.1 (primary source, confirmed): pi_{O(2n)}(chi_O(lam)) = 0 or
# +- chi_{O(2n)}(mu), mu obtained from lam "by the procedure for type D". Our scan lost the
# explicit procedure to OCR, King 1971 is not open access, and neither 2208.05526 (Jing-Li-Wang,
# skew symplectic/orthogonal Schur functions) nor 2209.00767 restates it.
#
# But we do not need the general rule -- we need it for OUR alphabet, and there it is computable.
# o_nu(A) is obtained honestly from the Schur expansion of the universal character, since the
# specialization homomorphism IS evaluation at the alphabet and s_mu(A) = 0 whenever ell(mu) > N.
#
# My earlier reading ("non-standard => 0, otherwise o_nu = -o_{nu*}") was fitted to nu with
# |nu| <= 8 and it is FALSE -- unstable_criterion.sage disagreed with the exact object in both
# directions. So this time: tabulate EVERY nu in a complete range, classify by the invariants the
# rule is supposed to use (nu'_1, nu'_1 + nu'_2, self-associate or not), and read the rule off the
# table rather than guess it.
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

def conj12(nu, N):
    c = list(Partition(list(nu)).conjugate())
    c1 = c[0] if len(c) > 0 else 0
    c2 = c[1] if len(c) > 1 else 0
    return c1, c2

def associate(nu, N):
    c = list(Partition(list(nu)).conjugate())
    first = N - (c[0] if c else 0)
    rest = c[1:]
    if first < 0 or first < (rest[0] if rest else 0):
        return None
    return tuple(Partition([first] + rest).conjugate())

for r in [1, 2]:
    N = 2 * r + 2
    half = N // 2
    R, A = alphabet(r)
    print("=" * 76)
    print("r=%d  N=%d  half=%d   the type-D table on our coset alphabet" % (r, N, half))
    vals = {}
    DMAX = 10 if r == 1 else 8
    for d in range(0, DMAX + 1):
        for nu in Partitions(d):
            vals[tuple(nu)] = o_at(tuple(nu), R, A)

    # classify
    buckets = {}
    for nu, v in vals.items():
        c1, c2 = conj12(nu, N)
        std = (c1 + c2 <= N)
        selfassoc = (c1 == half)
        key = ('standard' if std else 'NON-standard',
               'selfassoc' if selfassoc else ('c1<half' if c1 < half else 'c1>half'))
        buckets.setdefault(key, [0, 0])
        buckets[key][0] += 1
        if v == 0:
            buckets[key][1] += 1
    print("  bucket                                  count   of which o_nu(A)=0")
    for k in sorted(buckets, key=str):
        n, z = buckets[k]
        print("    %-38s %5d   %5d %s" % (str(k), n, z, "(ALL)" if z == n else ""))

    # for the nonzero ones with c1 > half, is o_nu = +- o_{associate}?
    print("  labels with c1 > half and o_nu(A) != 0 : does the associate relation hold?")
    ok = fail = noass = 0
    examples = []
    for nu, v in vals.items():
        c1, c2 = conj12(nu, N)
        if v == 0 or c1 <= half:
            continue
        st = associate(nu, N)
        if st is None or st not in vals:
            noass += 1
            if len(examples) < 6:
                examples.append((nu, 'no associate in range', None))
            continue
        w = vals[st]
        if v == -w:
            ok += 1
        elif v == w:
            ok += 1
        else:
            fail += 1
            if len(examples) < 6:
                examples.append((nu, st, 'MISMATCH'))
    print("    associate relation holds: %d ; fails: %d ; associate unavailable: %d"
          % (ok, fail, noass))
    for e in examples:
        print("      ", e)

    # the honest statement: which invariants predict o_nu(A) = 0 ?
    print("  zero-prediction test: is  (c1+c2 > N) or (c1 == half)  exactly the zero set?")
    wrong = [(nu, conj12(nu, N), vals[nu] == 0) for nu in vals
             if ((conj12(nu, N)[0] + conj12(nu, N)[1] > N) or (conj12(nu, N)[0] == half))
             != (vals[nu] == 0)]
    print("    counterexamples: %d" % len(wrong))
    for w in wrong[:10]:
        print("      nu=%-20s (c1,c2)=%s  o_nu(A)==0: %s" % (str(w[0]), str(w[1]), w[2]))
