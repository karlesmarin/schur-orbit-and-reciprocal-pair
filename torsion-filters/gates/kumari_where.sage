# kumari_where.sage -- WHERE does the negative multiplicity sit?
#
# The diagnosis offered for the failure of Theorem 4.4 of the withdrawn arXiv:2211.14093 is that
# the proof's Jacobi-Trudi step loses sign control at the INTERMEDIATE divisors of t: a determinant
# of entries that are nonnegative combinations sum_j j a_j need not itself be one, because a
# determinant subtracts.  A prime t has no intermediate divisor, which is why it survives.
#
# That diagnosis makes a sharp, falsifiable prediction: in EVERY failure, the negative
# Alexandersson-Amini multiplicity c_d sits at a PROPER divisor 1 < d < t -- never at c_1
# (which is f(omega), nonnegative by her own Corollary 6.6 under the sign hypothesis) and never
# at c_t (the top one).  This file checks that over every failure it can find.
#
# Authors: Carles Marin, Claude (AI assistant).

def betaset(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - (j + 1) for j in range(N)]

def sgn_sigma(lam, N, t):
    w = [v % t for v in betaset(lam, N)]
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

print("=" * 88)
print("In every failure, WHERE is the negative multiplicity?")
print("=" * 88)
print("%-4s %-9s %-9s %-14s %-14s %-14s" % ("t", "pairs", "failures", "neg at c_1", "neg at PROPER d", "neg at c_t"))
print("-" * 88)

SIZE = {4: 12, 6: 11, 8: 10, 9: 10}
grand = [0, 0, 0]
for t in (4, 6, 8, 9):
    K = CyclotomicField(t); w = K.gen()
    divs = sorted(ZZ(t).divisors())
    at1 = atmid = attop = 0
    tot = fails = 0
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
                f = f_poly(L, MM, tn)
                cs = {}
                broke = False
                for d in divs:
                    val = f(w ** d)
                    if val not in QQ:
                        broke = True
                        break
                    rem = QQ(val) - sum(j * cs[j] for j in cs if d % j == 0)
                    if rem % d != 0:
                        broke = True
                        break
                    cs[d] = ZZ(rem / d)
                if broke:
                    fails += 1
                    continue
                neg = [d for d in cs if cs[d] < 0]
                if neg:
                    fails += 1
                    if 1 in neg:
                        at1 += 1
                    if any(1 < d < t for d in neg):
                        atmid += 1
                    if t in neg:
                        attop += 1
    grand[0] += at1; grand[1] += atmid; grand[2] += attop
    print("%-4d %-9d %-9d %-14d %-14d %-14d" % (t, tot, fails, at1, atmid, attop))

print("-" * 88)
print("%-4s %-9s %-9s %-14d %-14d %-14d" % ("all", "", "", grand[0], grand[1], grand[2]))
print()
print("PREDICTION: the middle column carries every failure; the outer two are 0.")
print("RESULT    : c_1 %s, proper divisor %s, c_t %s"
      % (grand[0], grand[1], grand[2]))
print()
print("CONTROL -- c_1 = f(omega) must be >= 0 whenever the sign hypothesis holds,")
print("           by her own Corollary 6.6 (Section 3, NOT the withdrawn section).")
bad = tot = 0
for t in (4, 6, 8, 9):
    K = CyclotomicField(t); w = K.gen()
    for L in parts_upto(9, t):
        sl = sgn_sigma(L, t, t)
        for M in parts_upto(sum(L), t):
            MM = list(M) + [0] * (len(L) - len(M))
            if len(M) > len(L) or any(MM[i] > L[i] for i in range(len(L))):
                continue
            if sgn_sigma(MM, t, t) != sl:
                continue
            v = f_poly(L, MM, t)(w)
            tot += 1
            if v in QQ and QQ(v) < 0:
                bad += 1
print("           %d pairs, %d with f(omega) < 0" % (tot, bad))
