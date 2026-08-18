# kumari_root.sage -- WHY does Theorem 4.4 of the withdrawn arXiv:2211.14093 fail exactly at the
# composites, and what is the repair?
#
# THE ARGUMENT.  f(q) = s_{lambda/mu}(1,q,...,q^{tn-1}).  For d | t put e = t/d.  The multiset
# {omega^{dk} : 0 <= k < tn} is the set of e-th roots of unity, each repeated dn times.  That is
# the alphabet of her Corollary 6.6 with t replaced by e and n by dn, so
#
#     f(omega^d) = sgn(sigma^(e)_lambda) sgn(sigma^(e)_mu) prod_{i<e} s_{lambda^(i,e)/mu^(i,e)}(1^{dn})
#
# where sigma^(e) sorts beta by residue MODULO e.  The product is a count, hence >= 0; the sign is
# the only thing that can make f(omega^d) negative.
#
#   d = 1  -> e = t : the sign IS her hypothesis.        controlled
#   d = t  -> e = 1 : f(1) = |SSYT|.                     trivially fine
#   1<d<t  -> 1<e<t : the sign is taken modulo e, and her hypothesis says nothing about it.
#
# A prime t has no intermediate divisor.  That is the whole of it: the hypothesis controls ONE
# channel of the tau(t) channels the Alexandersson-Amini criterion tests.
#
# THREE TESTS.
#   R1  the identity above, checked directly at every divisor -- is the mod-e sign really what
#       governs the sign of f(omega^d)?
#   R2  is every failure explained by a NEGATIVE mod-e sign at some intermediate divisor?
#   R3  THE REPAIR.  Impose sgn(sigma^(e)_lambda) = sgn(sigma^(e)_mu) for EVERY e | t, e > 1.
#       Prediction: zero failures at every t, composite included.
#       Control: the repaired hypothesis must be strictly stronger, i.e. it must actually
#       throw away pairs that the original admitted, or it is not a repair but a relabelling.
#
# Authors: Carles Marin, Claude (AI assistant).

def betaset(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - (j + 1) for j in range(N)]

def sgn_sigma(lam, N, e):
    """sorting sign of the residue word MODULO e"""
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

TS = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
SIZE = {2: 12, 3: 12, 4: 12, 5: 11, 6: 11, 7: 10, 8: 10, 9: 10, 10: 10, 12: 10}

def pairs(t):
    for n in (1, 2):
        tn = t * n
        if tn > 12:
            continue
        for L in parts_upto(SIZE[t], tn):
            for M in parts_upto(sum(L), tn):
                MM = list(M) + [0] * (len(L) - len(M))
                if len(M) > len(L) or any(MM[i] > L[i] for i in range(len(L))):
                    continue
                yield (L, MM, tn, n)

print("=" * 96)
print("R1  f(omega^d) and the sorting sign MODULO e = t/d")
print("=" * 96)
print("%-4s %-10s %-30s %-30s" % ("t", "checked", "sign(f(om^d)) = mod-e sign?", "f(om^d)=0 when e-cores differ?"))
print("-" * 96)
for t in TS:
    K = CyclotomicField(t); w = K.gen()
    agree = tot = zero_ok = zero_tot = 0
    for (L, M, tn, n) in pairs(t):
        if sgn_sigma(M, tn, t) != sgn_sigma(L, tn, t):
            continue
        f = f_poly(L, M, tn)
        for d in sorted(ZZ(t).divisors()):
            e = t // d
            if e == 1:
                continue
            val = f(w ** d)
            if val not in QQ:
                continue
            val = QQ(val)
            s = sgn_sigma(L, tn, e) * sgn_sigma(M, tn, e)
            pe = [sum(1 for v in betaset(L, tn) if v % e == i) for i in range(e)]
            qe = [sum(1 for v in betaset(M, tn) if v % e == i) for i in range(e)]
            tot += 1
            if pe != qe:
                zero_tot += 1
                if val == 0:
                    zero_ok += 1
                agree += 1                      # vacuously: value is 0
            elif val == 0 or (val > 0) == (s > 0):
                agree += 1
    print("%-4d %-10d %-30s %-30s" % (t, tot, "%d/%d" % (agree, tot), "%d/%d" % (zero_ok, zero_tot)))

print()
print("=" * 96)
print("R2 + R3   the original hypothesis vs the repaired one")
print("=" * 96)
print("%-4s %-6s %-12s %-10s %-14s %-10s %-14s"
      % ("t", "prime", "pairs(orig)", "fail", "pairs(repair)", "fail", "thrown away"))
print("-" * 96)
for t in TS:
    K = CyclotomicField(t); w = K.gen()
    o_tot = o_bad = r_tot = r_bad = 0
    for (L, M, tn, n) in pairs(t):
        if sgn_sigma(M, tn, t) != sgn_sigma(L, tn, t):
            continue
        o_tot += 1
        cs = mults(f_poly(L, M, tn), t, w)
        bad = (cs is None) or any(cs[d] < 0 for d in cs)
        if bad:
            o_bad += 1
        repaired = all(sgn_sigma(L, tn, e) == sgn_sigma(M, tn, e)
                       for e in ZZ(t).divisors() if e > 1)
        if repaired:
            r_tot += 1
            if bad:
                r_bad += 1
    print("%-4d %-6s %-12d %-10d %-14d %-10d %-14d"
          % (t, ZZ(t).is_prime(), o_tot, o_bad, r_tot, r_bad, o_tot - r_tot))

print()
print("PREDICTION: 'fail' under the repaired hypothesis is 0 at EVERY t.")
print("CONTROL   : 'thrown away' must be 0 exactly at prime t and positive at composite t,")
print("            or the repair is a relabelling rather than a strengthening.")
