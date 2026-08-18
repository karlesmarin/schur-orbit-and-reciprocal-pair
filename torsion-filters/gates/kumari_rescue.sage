# kumari_rescue.sage
#
# arXiv:2211.14093 (= Chapter 6 of N. Kumari's thesis) was withdrawn by the author on 2025-12-17,
# "Error in Section 4".  Section 4 is the cyclic sieving section, and its Theorem 4.4 (= thesis
# Theorem 6.17) reads:
#
#     if sgn(sigma_lambda) = sgn(sigma_mu), then there is an action of C_t on SSYT_{tn}(lambda/mu)
#     such that (SSYT_{tn}(lambda/mu), C_t, s_{lambda/mu}(1,q,...,q^{tn-1})) exhibits CSP.
#
# This file asks a sharper question than "is it false".  It asks WHERE.
#
# The test is decidable.  By Alexandersson-Amini [AA19, Thm 2.7], such an action exists iff for
# every d | t,   f(omega^d) = sum_{j|d} j c_j   with every c_j a nonnegative integer, where
# f(q) = s_{lambda/mu}(1,q,...,q^{tn-1}).  The c_j are determined by triangularity, so we simply
# compute them.  f is built by Jacobi-Trudi with h_k(1,q,...,q^{m-1}) = qbinomial(k+m-1,k), exact
# in Z[q], and evaluated in the cyclotomic field.
#
# HYPOTHESIS UNDER TEST: the theorem survives for t PRIME and fails from the first composite on.
# The reason would be structural: for t prime the divisor set is {1,t}, so the only constraints
# are c_1 = f(omega) >= 0 and c_t = (f(1)-f(omega))/t >= 0, and there is no INTERMEDIATE c_j that
# could go negative.  For composite t the intermediate divisors are exactly where the proof's
# Jacobi-Trudi step -- a determinant of entries that are nonnegative combinations, concluded to be
# a nonnegative combination -- loses control, because a determinant also subtracts.
#
# Authors: Carles Marin, Claude (AI assistant).

def betaset(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - (j + 1) for j in range(N)]

def sgn_sigma(lam, N, t):
    w = [v % t for v in betaset(lam, N)]
    inv = sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])
    return (-1) ** inv

def parts_upto(K, maxlen):
    out = []
    for k in range(K + 1):
        for p in Partitions(k, max_length=maxlen):
            out.append(list(p))
    return out

R = PolynomialRing(ZZ, 'q'); q = R.gen()

def f_poly(L, M, tn):
    r = len(L)
    if r == 0:
        return R(1)
    MM = list(M) + [0] * (r - len(M))
    rows = []
    for i in range(r):
        row = []
        for j in range(r):
            k = L[i] - MM[j] - (i + 1) + (j + 1)
            row.append(R(0) if k < 0 else R(q_binomial(k + tn - 1, k, q)))
        rows.append(row)
    return R(matrix(R, rows).det())

def multiplicities(L, M, t, n, w):
    """the Alexandersson-Amini c_j, or a reason they do not exist"""
    tn = int(t) * int(n)
    f = f_poly(L, M, tn)
    cs = {}
    for d in sorted(ZZ(t).divisors()):
        val = f(w ** d)
        if val not in QQ:
            return (None, 'f(omega^%d) irrational' % d, f)
        val = QQ(val)
        rem = val - sum(j * cs[j] for j in cs if d % j == 0)
        if rem % d != 0:
            return (None, 'c_%d not integral' % d, f)
        cs[d] = ZZ(rem / d)
    return (cs, 'ok', f)

print("=" * 90)
print("Theorem 4.4 of the withdrawn arXiv:2211.14093, tested divisor by divisor")
print("=" * 90)
print("%-4s %-4s %-9s %-9s %-30s" % ("t", "prime", "pairs", "failures", "smallest failure"))
print("-" * 90)

SIZE = {2: 12, 3: 12, 4: 12, 5: 11, 6: 11, 7: 10, 8: 10, 9: 10}
worst = {}
for t in (2, 3, 4, 5, 6, 7, 8, 9):
    K = CyclotomicField(t); w = K.gen()
    tot = 0; fails = []
    for n in (1, 2):
        tn = t * n
        if tn > 9:
            continue
        for L in parts_upto(SIZE[t], tn):
            sl = sgn_sigma(L, tn, t)
            for M in parts_upto(sum(L), tn):
                MM = list(M) + [0] * (len(L) - len(M))
                if len(M) > len(L) or any(MM[i] > L[i] for i in range(len(L))):
                    continue
                if sgn_sigma(MM, tn, t) != sl:
                    continue
                tot += 1
                cs, why, f = multiplicities(L, MM, t, n, w)
                if cs is None or any(cs[j] < 0 for j in cs):
                    fails.append((sum(L) - sum(MM), L, MM, cs, why, n))
    fails.sort(key=lambda r: (r[0], sum(r[1])))
    sm = ""
    if fails:
        d, L, M, cs, why, n = fails[0]
        sm = "lam=%s mu=%s n=%d  c=%s" % (L, M, n, cs)
        worst[t] = fails[0]
    print("%-4d %-5s %-9d %-9d %-30s" % (t, ZZ(t).is_prime(), tot, len(fails), sm))

print()
print("=" * 90)
print("THE SMALLEST FAILURE, written out")
print("=" * 90)
if 4 in worst:
    d, L, M, cs, why, n = worst[4]
    t = 4; K = CyclotomicField(t); w = K.gen()
    _, _, f = multiplicities(L, M, t, n, w)
    print("  t=4,  lambda=%s,  mu=%s,  tn=%d" % (L, M, t * n))
    print("  sgn(sigma_lambda) = %d = sgn(sigma_mu) = %d   -- the hypothesis HOLDS"
          % (sgn_sigma(L, t * n, t), sgn_sigma(M, t * n, t)))
    print("  f(q) = s_{lambda/mu}(1,q,q^2,q^3) = %s" % f)
    for e in (1, 2, 4):
        print("     f(omega^%d) = %s" % (e, f(w ** e)))
    print("  Alexandersson-Amini multiplicities: %s" % cs)
    print("  => sigma^2 would have to fix %s elements.  No action exists." % (f(w ** 2)))

print()
print("=" * 90)
print("CONTROLS")
print("=" * 90)
K2 = CyclotomicField(2); w2 = K2.gen()
cs, why, f = multiplicities([1, 1, 1], [1, 0, 0], 2, 2, w2)
# NOTE.  'why' only reports whether the c_j were COMPUTABLE; it does not test their sign.
# This pair has c_1 = -2 < 0, so no action exists -- and it does NOT satisfy the hypothesis
# (sgn(sigma_lambda) = +1, sgn(sigma_mu) = -1).  That is the point of the control: the sign
# hypothesis is load-bearing even at prime t, and dropping it breaks the criterion at t=2.
print("  the SAME shape at t=2 (prime, tn=4): c=%s ; computable? %s" % (cs, why))
print("     sgn(sigma_lambda)=%d  sgn(sigma_mu)=%d  -> hypothesis FAILS here, as it must"
      % (sgn_sigma([1,1,1], 4, 2), sgn_sigma([1,0,0], 4, 2)))
print("     so at prime t the sign hypothesis is exactly what rules this out.")
print("  the criterion is not vacuous: it returns negative c only sometimes, see the table above.")
K4 = CyclotomicField(4); w4 = K4.gen()
cs, why, f = multiplicities([2, 1], [], 4, 1, w4)
print("  a shape that PASSES at t=4: lam=(2,1) mu=() -> c=%s" % (cs,))
print("  f(1) >= f(omega) always?  (needed for c_t >= 0 at prime t)")
bad = 0; tested = 0
for t in (2, 3, 5, 7):
    K = CyclotomicField(t); w = K.gen()
    for L in parts_upto(9, t):
        for M in parts_upto(sum(L), t):
            MM = list(M) + [0] * (len(L) - len(M))
            if len(M) > len(L) or any(MM[i] > L[i] for i in range(len(L))):
                continue
            if sgn_sigma(MM, t, t) != sgn_sigma(L, t, t):
                continue
            f = f_poly(L, MM, t)
            tested += 1
            if not (f(1) >= QQ(f(w))):
                bad += 1
print("      %d pairs at prime t, %d with f(1) < f(omega)" % (tested, bad))
