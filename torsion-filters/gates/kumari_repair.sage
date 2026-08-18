# kumari_repair.sage -- does the REPAIRED theorem still do the job the original was written for?
#
# Repaired hypothesis (H*):  sgn(sigma^(e)_lambda) = sgn(sigma^(e)_mu) for EVERY e | t, e > 1.
#
# Q1  Does Lee-Oh's hypothesis imply H*?   Lee-Oh [LO22] assume t | lambda_i - mu_i for all i.
#     Then beta_i(lambda) - beta_i(mu) = lambda_i - mu_i is divisible by t, hence by every e | t,
#     so lambda and mu have the SAME residue word mod e and sigma^(e)_lambda = sigma^(e)_mu.
#     So H* should hold trivially.  Checked, not assumed.
#
# Q2  Is H* STRICTLY weaker than Lee-Oh?  If not, the repair rescues nothing beyond Lee-Oh and the
#     paper's contribution is gone.  We count the pairs satisfying H* but NOT Lee-Oh.
#
# Q3  Is H* equivalent to imposing the sign condition only at the PRIME divisors of t?  That would
#     be a cleaner statement of the same theorem.
#
# Authors: Carles Marin, Claude (AI assistant).

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

R = PolynomialRing(ZZ, 'q'); q = R.gen()

def f_poly(L, M, tn):
    r = len(L)
    if r == 0:
        return R(1)
    MM = list(M) + [0] * (r - len(M))
    return R(matrix(R, [[R(0) if (L[i] - MM[j] - i + j) < 0
                         else R(q_binomial(L[i] - MM[j] - i + j + tn - 1, L[i] - MM[j] - i + j, q))
                         for j in range(r)] for i in range(r)]).det())

def mults(f, t, w):
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

SIZE = {4: 12, 6: 11, 8: 10, 9: 10, 10: 10, 12: 10}
print("=" * 100)
print("%-4s %-9s %-11s %-11s %-16s %-14s %-14s"
      % ("t", "pairs", "LeeOh", "H*", "H* and not LO", "LO not H*", "H* vs primes-only"))
print("=" * 100)
for t in (4, 6, 8, 9, 10, 12):
    K = CyclotomicField(t); w = K.gen()
    primes = [p for p in ZZ(t).prime_divisors()]
    tot = lo = hs = hs_not_lo = lo_not_hs = prime_mismatch = 0
    hs_fail = 0
    for n in (1, 2):
        tn = t * n
        if tn > 12:
            continue
        for L in parts_upto(SIZE[t], tn):
            for M in parts_upto(sum(L), tn):
                MM = list(M) + [0] * (len(L) - len(M))
                if len(M) > len(L) or any(MM[i] > L[i] for i in range(len(L))):
                    continue
                tot += 1
                LO = all((L[i] - MM[i]) % t == 0 for i in range(len(L)))
                HS = all(sgn_sigma(L, tn, e) == sgn_sigma(MM, tn, e)
                         for e in ZZ(t).divisors() if e > 1)
                PR = all(sgn_sigma(L, tn, p) == sgn_sigma(MM, tn, p) for p in primes) and \
                     sgn_sigma(L, tn, t) == sgn_sigma(MM, tn, t)
                if LO: lo += 1
                if HS: hs += 1
                if HS and not LO: hs_not_lo += 1
                if LO and not HS: lo_not_hs += 1
                if HS != PR: prime_mismatch += 1
                if HS:
                    cs = mults(f_poly(L, MM, tn), t, w)
                    if cs is None or any(cs[d] < 0 for d in cs):
                        hs_fail += 1
    print("%-4d %-9d %-11d %-11d %-16d %-14d %-14d"
          % (t, tot, lo, hs, hs_not_lo, lo_not_hs, prime_mismatch))
    print("     -> H* pairs failing the CSP criterion: %d   (must be 0)" % hs_fail)

print()
print("=" * 100)
print("Q1  Lee-Oh => H* :   the 'LO not H*' column must be 0 everywhere.")
print("Q2  H* strictly weaker : the 'H* and not LO' column must be large, or the repair is empty.")
print("Q3  primes suffice? :  'H* vs primes-only' counts disagreements; 0 means the condition at")
print("     the prime divisors (plus t itself) already implies it at every divisor.")
