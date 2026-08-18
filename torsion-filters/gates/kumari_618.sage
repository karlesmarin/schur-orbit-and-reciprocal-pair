# kumari_618.sage -- does the SUPER version, Theorem 6.18 of the thesis (= Theorem 4.5 of the
# withdrawn arXiv:2211.14093), fail the same way?
#
# Theorem 6.18.  Suppose t is ODD.  Let lambda, mu have length at most tn with
# sgn(sigma_lambda) = sgn(sigma_mu).  Then there is an action of C_t on SSYT_{tn/tm}(lambda/mu)
# such that
#     ( SSYT_{tn/tm}(lambda/mu), C_t, hs_{lambda/mu}(1,q,...,q^{tn-1} / 1,q,...,q^{tm-1}) )
# exhibits the cyclic sieving phenomenon.
#
# PREDICTION.  "t odd" is not "t prime".  The divisor lattice of 9 is {1,3,9}, so 9 HAS an
# intermediate divisor, and the mechanism that kills Theorem 4.4 at composite t should kill this
# one at odd COMPOSITE t while leaving odd primes alone.  Her Theorem 6.5 (the hook analogue of
# Corollary 6.6) gives, at q = omega^d with e = t/d,
#     hs = sgn(sigma^(e)_lambda) sgn(sigma^(e)_mu) prod_i hs_{lambda^(i,e)/mu^(i,e)}(...),
# and for t odd the sign (-1)^{t-1} in her Theorem 6.1 is +1, so nothing else intervenes.
#
# INGREDIENTS (her (2.5.2) and (2.5.6)):
#     H_k(X/Y) = sum_{l=0..k} h_l(X) e_{k-l}(Y)
#     hs_{lambda/mu}(X/Y) = det( H_{lambda_i - mu_j - i + j}(X/Y) )
#     h_l(1,q,...,q^{a-1}) = qbinomial(l+a-1, l)
#     e_j(1,q,...,q^{b-1}) = q^{binom(j,2)} qbinomial(b, j)
#
# CONTROL.  At m = 0 the primed alphabet is empty and hs must reduce to s.  If the m=0 column does
# not reproduce the Schur results (0 failures at prime t, failures at composite t), the
# implementation is wrong and nothing below means anything.
#
# Authors: Carles Marin, Claude (AI assistant).

R = PolynomialRing(ZZ, 'q'); q = R.gen()

def betaset(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - (j + 1) for j in range(N)]

def sgn_sigma(lam, N, e):
    w = [v % e for v in betaset(lam, N)]
    return (-1) ** sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])

def parts_upto(K, maxlen):
    out = []
    for k in range(K + 1):
        for p in Partitions(k, max_length=maxlen):
            out.append(list(p))
    return out

def h_spec(l, a):
    """h_l(1,q,...,q^{a-1})"""
    if l < 0 or a <= 0:
        return R(0) if l != 0 else R(1 if l == 0 else 0)
    return R(q_binomial(l + a - 1, l, q))

def e_spec(j, b):
    """e_j(1,q,...,q^{b-1})"""
    if j < 0 or j > b:
        return R(0)
    return R(q ** binomial(j, 2) * q_binomial(b, j, q))

def H_spec(k, a, b):
    """H_k(X/Y) with |X|=a, |Y|=b, both principally specialized"""
    if k < 0:
        return R(0)
    tot = R(0)
    for l in range(k + 1):
        hl = R(1) if l == 0 else h_spec(l, a)
        if a == 0 and l > 0:
            hl = R(0)
        tot += hl * e_spec(k - l, b)
    return tot

def hs_poly(lam, mu, a, b):
    L = list(lam); r = len(L)
    if r == 0:
        return R(1)
    M = list(mu) + [0] * (r - len(mu))
    return R(matrix(R, [[H_spec(L[i] - M[j] - (i + 1) + (j + 1), a, b)
                         for j in range(r)] for i in range(r)]).det())

def mults(f, t):
    K = CyclotomicField(t); w = K.gen()
    cs = {}
    for d in sorted(ZZ(t).divisors()):
        val = f(w ** d)
        if val not in QQ:
            return None
        rem = QQ(val) - sum(j * cs[j] for j in cs if d % j == 0)
        if rem % d != 0:
            return None
        cs[d] = ZZ(rem / d)
    return cs

print("=" * 100)
print("CONTROL FIRST: m = 0 must reproduce the Schur case (Theorem 4.4)")
print("=" * 100)
print("%-4s %-7s %-10s %-10s" % ("t", "prime", "pairs", "failures"))
for t in (3, 4, 5, 9):
    n = 1; tn = t * n
    tot = bad = 0
    for L in parts_upto(9, tn):
        for mu in parts_upto(sum(L), tn):
            M = list(mu) + [0] * (len(L) - len(mu))
            if len(mu) > len(L) or any(M[i] > L[i] for i in range(len(L))):
                continue
            if sgn_sigma(M, tn, t) != sgn_sigma(L, tn, t):
                continue
            tot += 1
            cs = mults(hs_poly(L, M, tn, 0), t)
            if cs is None or any(cs[d] < 0 for d in cs):
                bad += 1
    print("%-4d %-7s %-10d %-10d" % (t, ZZ(t).is_prime(), tot, bad))
print("  (must be: 0 failures at t=3,5 ; failures at t=4,9 -- otherwise the code is wrong)")

print()
print("=" * 100)
print("THEOREM 6.18, t ODD, with a genuine primed alphabet")
print("=" * 100)
print("%-4s %-7s %-4s %-4s %-10s %-10s %-14s %-10s"
      % ("t", "prime", "n", "m", "pairs", "failures", "repaired H*", "fail H*"))
print("-" * 100)
for t in (3, 5, 7, 9, 15):
    for (n, m) in ((1, 1), (1, 2)):
        tn, tm = t * n, t * m
        if tn > 9 or tm > 15:
            continue
        tot = bad = rtot = rbad = 0
        smallest = None
        for L in parts_upto(8 if t < 9 else 7, tn):
            for mu in parts_upto(sum(L), tn):
                M = list(mu) + [0] * (len(L) - len(mu))
                if len(mu) > len(L) or any(M[i] > L[i] for i in range(len(L))):
                    continue
                if sgn_sigma(M, tn, t) != sgn_sigma(L, tn, t):
                    continue
                tot += 1
                cs = mults(hs_poly(L, M, tn, tm), t)
                isbad = (cs is None) or any(cs[d] < 0 for d in cs)
                if isbad:
                    bad += 1
                    if smallest is None:
                        smallest = (L, M, cs)
                if all(sgn_sigma(L, tn, e) == sgn_sigma(M, tn, e)
                       for e in ZZ(t).divisors() if e > 1):
                    rtot += 1
                    if isbad:
                        rbad += 1
        print("%-4d %-7s %-4d %-4d %-10d %-10d %-14d %-10d"
              % (t, ZZ(t).is_prime(), n, m, tot, bad, rtot, rbad))
        if smallest:
            print("       smallest failure: lambda=%s mu=%s  c=%s"
                  % (smallest[0], smallest[1], smallest[2]))

print()
print("=" * 100)
print("READING")
print("=" * 100)
print("  If t=3,5,7 show 0 failures and t=9,15 show failures, then 'odd' in Theorem 6.18 is")
print("  doing the job 'prime' does in Theorem 4.4, and the same repair -- the sign condition at")
print("  EVERY divisor e | t -- should send the last column to 0.")
