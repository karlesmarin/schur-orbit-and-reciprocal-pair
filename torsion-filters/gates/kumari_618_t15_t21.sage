# kumari_618_t15_t21.sage -- the supertableaux theorem at the next two odd composites.
#
# kumari_618.sage found that Theorem 6.18 of the thesis (= Theorem 4.5 of the withdrawn
# arXiv:2211.14093) is clean at the odd PRIMES 3, 5, 7 and fails at t = 9, with the negative
# multiplicity sitting at c_3 -- the intermediate divisor of 9.  That is the same mechanism that
# breaks Theorem 4.4: the hypothesis fixes the sorting sign modulo t, and the criterion tests one
# condition per divisor d | t, at which the relevant sign is taken modulo e = t/d.
#
# t = 9 has one intermediate divisor.  t = 15 has two (3 and 5) and t = 21 has two (3 and 7).  So:
#
#   PREDICTION 1  both fail, and every negative multiplicity sits at an intermediate divisor,
#                 never at c_1 and never at c_t.
#   PREDICTION 2  the repaired hypothesis -- sgn(sigma^(e)_lambda) = sgn(sigma^(e)_mu) for every
#                 e | t, e > 1 -- sends the failures to 0.
#
# CONTROL.  m = 0 must reproduce the Schur case (Theorem 4.4) at the same t.  If the m = 0 column
# is clean where the m = 1 column is not, the implementation is treating the two cases
# inconsistently and nothing here is evidence.
#
# H_k(X/Y) = sum_l h_l(X) e_{k-l}(Y), Jacobi-Trudi of size l(lambda); h_l(1,q,...,q^{a-1}) =
# qbinomial(l+a-1, l) and e_j(1,q,...,q^{b-1}) = q^{binom(j,2)} qbinomial(b, j).  Exact in Z[q].
#
# Authors: Carles Marin, Claude (AI assistant).

import sys
sys.stdout.reconfigure(line_buffering=True)

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
    if l == 0:
        return R(1)
    if l < 0 or a <= 0:
        return R(0)
    return R(q_binomial(l + a - 1, l, q))


def e_spec(j, b):
    if j == 0:
        return R(1)
    if j < 0 or j > b:
        return R(0)
    return R(q ** binomial(j, 2) * q_binomial(b, j, q))


def H_spec(k, a, b):
    if k < 0:
        return R(0)
    if b == 0:
        return h_spec(k, a)
    return sum((h_spec(l, a) * e_spec(k - l, b) for l in range(k + 1)), R(0))


def f_poly(lam, mu, a, b):
    L = list(lam); r = len(L)
    if r == 0:
        return R(1)
    M = list(mu) + [0] * (r - len(mu))
    return R(matrix(R, [[H_spec(L[i] - M[j] - (i + 1) + (j + 1), a, b)
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


print("=" * 100)
print("Theorem 6.18 at the next two odd composites.  n = m = 1, |lambda| <= 8.")
print("=" * 100)
print("%-4s %-14s %-4s %-8s %-9s %-9s %-11s %-9s"
      % ("t", "divisors>1", "m", "pairs", "failures", "neg at c_1", "neg at 1<d<t", "H* fails"))
print("-" * 100)

for t in (9, 15, 21):
    K = CyclotomicField(t); w = K.gen()
    divs = [d for d in ZZ(t).divisors() if d > 1]
    tn = t
    shapes = parts_upto(8, tn)
    for m in (0, 1):
        tm = t * m
        tot = fails = at1 = atmid = hstar_tot = hstar_bad = 0
        smallest = None
        for L in shapes:
            sl = sgn_sigma(L, tn, t)
            for mu in parts_upto(sum(L), tn):
                M = list(mu) + [0] * (len(L) - len(mu))
                if len(mu) > len(L) or any(M[i] > L[i] for i in range(len(L))):
                    continue
                if sgn_sigma(M, tn, t) != sl:
                    continue
                tot += 1
                cs = mults(f_poly(L, M, tn, tm), t, w)
                isbad = (cs is None) or any(cs[d] < 0 for d in cs)
                if isbad:
                    fails += 1
                    if cs is not None:
                        neg = [d for d in cs if cs[d] < 0]
                        if 1 in neg:
                            at1 += 1
                        if any(1 < d < t for d in neg):
                            atmid += 1
                        if smallest is None:
                            smallest = (L, M, cs)
                if all(sgn_sigma(L, tn, e) == sgn_sigma(M, tn, e) for e in divs):
                    hstar_tot += 1
                    if isbad:
                        hstar_bad += 1
        print("%-4d %-14s %-4d %-8d %-9d %-9d %-11d %-9d"
              % (t, str(divs), m, tot, fails, at1, atmid, hstar_bad))
        if smallest:
            print("       smallest: lambda=%-16s mu=%-16s c=%s"
                  % (smallest[0], smallest[1], smallest[2]))
        if m == 1:
            print("       (H* admitted %d of the %d pairs)" % (hstar_tot, tot))

print()
print("=" * 100)
print("READING")
print("=" * 100)
print("  PREDICTION 1 holds if 'neg at c_1' is 0 and 'neg at 1<d<t' equals 'failures'.")
print("  PREDICTION 2 holds if 'H* fails' is 0 in every row.")
print("  CONTROL: the m=0 rows must fail too -- they are Theorem 4.4 at the same t.")
