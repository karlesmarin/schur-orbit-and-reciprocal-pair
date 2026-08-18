# -*- coding: utf-8 -*-
# Simplify the day's formulas to canonical form and cross them.
#
# Everything today has been a statement about ONE object: an involution sigma on Z/t and its orbit
# space.  Written in that language the formulas collapse, and two of them then say opposite things,
# which is the crossing worth checking.
#
# ---- the canonical forms -----------------------------------------------------------------------
#
# F1  NPP's type-C count, simplified.  With n_k the number of m_j in the folded class k,
#         N(m) = #{i<j: t|m_i-m_j} + #{i<j: t|m_i+m_j} + #{i: t|2m_i}
#     For a NON-FIXED class k the n_k elements split a of them = k and b = -k mod t, and
#     C(a,2)+C(b,2)+ab = C(a+b,2), independent of the split.  For a FIXED class every element is its
#     own negative, so the first two terms are each C(n_k,2) and the third adds n_k, giving
#     2C(n_k,2)+n_k = n_k^2.  Hence
#                  N(m)  =  sum_{k non-fixed} C(n_k,2)  +  n_0^2  +  n_{t/2}^2.
#     The fixed classes cost a SQUARE and the others a binomial: that asymmetry is the whole reason
#     the count parts company with the covering condition.
#
# F2  the folded lemma, as a rank:  rank(frozen block) = #{distinct non-fixed classes present},
#     so  sp_nu != 0  <=>  rank is full  <=>  no non-fixed class is MISSED.
#
# F3  the GL criterion.  Let sigma_C(k) = C-k on Z/t and let E be the set of EXCESS classes.
#     Condition (i) says sigma_C preserves the excess values, hence E is sigma_C-stable.
#     Condition (ii) says both fixed points of sigma_C lie in E.  So (i)+(ii) is
#                  E is sigma_C-stable and contains BOTH fixed points of sigma_C,
#     and since a sigma_C-stable set is the two fixed points plus whole 2-element orbits,
#                  |E| is EVEN.                                    <-- a prediction, tested below
#     Also s := #fixed points of sigma_C inside E, and 2k = C mod t has gcd(2,t) solutions, so
#     s <= 2 always and s <= 1 for odd t.  Writing q = (-1)^{binom(s,2)} was over-parametrised:
#     with s in {0,1,2} that function is just the indicator of s = 2.
#
# ---- the crossing ------------------------------------------------------------------------------
#
# X   Both live on the orbit space of an involution on Z/t, and they ask OPPOSITE things:
#         GL side  (Phi_t = 0):   the FIXED classes must be OCCUPIED (by excess)
#         Sp side  (sp_nu = 0):   some NON-FIXED class must be EMPTY
#     If that is more than a pun there should be a measurable trace: on the GL side the vanishing
#     shapes should ALSO cover their non-fixed sigma_C-classes, i.e. the two conditions would be
#     the two halves of one statement.  If instead the non-fixed classes of a zero are arbitrary,
#     the resemblance is a coincidence of vocabulary and gets said so.
#
# CONTROLS able to fail:
#   C1  F1 checked as an identity against the literal triple sum, on every m in range.
#   C2  |E| even, on every zero.  A single odd one kills F3's reading.
#   C3  X measured, not asserted: the distribution of "non-fixed sigma_C classes missed" over the
#       zeros.  A spread distribution refutes the crossing.
#   C4  a control that can fail: the same distribution over the NON-vanishing concentric shapes.
#       If the two look alike, the crossing says nothing.
#
# Authors: Carles Marin, Claude (AI assistant).

TS = [4, 6, 8, 10, 12]
L = lcm(TS)
p = next_prime(10 ** 9)
while (p - 1) % L != 0:
    p = next_prime(p)
F = GF(p)
G = F.multiplicative_generator()
for t in TS:
    assert (p - 1) % t == 0
print("field GF(%d); guard -> PASS" % p)


def zeta(t):
    z = G ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t))
    return z


def cls(m, t):
    return min(m % t, (-m) % t)


def npp_literal(m, t):
    n = len(m)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (m[i] - m[j]) % t == 0:
                c += 1
            if (m[i] + m[j]) % t == 0:
                c += 1
        if (2 * m[i]) % t == 0:
            c += 1
    return c


def npp_closed(m, t):
    n = {}
    for x in m:
        n[cls(x, t)] = n.get(cls(x, t), 0) + 1
    tot = 0
    for k, v in n.items():
        if k in (0, t // 2):
            tot += v * v
        else:
            tot += v * (v - 1) // 2
    return tot


print("")
print("=" * 92)
print("C1  F1 as an identity:  N(m) = sum_{k non-fixed} C(n_k,2) + n_0^2 + n_{t/2}^2")
print("=" * 92)
bad = 0
tested = 0
for t in TS:
    for R in (t // 2, t // 2 + 2):
        for size in range(0, 13):
            for l in Partitions(size, max_length=R):
                nu = list(l) + [0] * (R - len(l))
                m = [nu[j] + R - j for j in range(R)]
                tested += 1
                if npp_literal(m, t) != npp_closed(m, t):
                    bad += 1
print("  %d weight vectors tested, mismatches: %d  ->  %s"
      % (tested, bad, "IDENTITY" if bad == 0 else "FALSE"))
if bad:
    raise SystemExit(1)


# ------------------------------------------------------------------- the GL side: zeros of Phi_t
def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


print("")
print("=" * 92)
print("C2/C3/C4  the GL side: |E| parity, and how the zeros treat their NON-fixed sigma_C classes")
print("=" * 92)
print("")
print("   t   r   N |lam|<=  ZEROS  |E| odd   zeros missing a non-fixed class   controls missing one")
print("  " + "-" * 92)

for t, r, MAX in ((4, 1, 26), (4, 2, 22), (4, 3, 20), (6, 2, 18), (6, 3, 16), (8, 2, 16)):
    N = t + 2 * r
    zt = zeta(t)
    RT = [zt ** k for k in range(t)]
    set_random_seed(5500 + 10 * t + r)
    PTS = []
    tries = 0
    while len(PTS) < 3 and tries < 600:
        tries += 1
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        al = RT + [y for x in zz for y in (x, 1 / x)]
        if len(set(al)) == N:
            PTS.append(al)
    POW = [[[a ** e for e in range(MAX + N + 1)] for a in al] for al in PTS]
    nz = odd = miss = 0
    ctrl = ctrl_miss = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            clsmap = {}
            for k, b in enumerate(beta):
                clsmap.setdefault(b % t, []).append(k)
            if len(clsmap) < t:
                continue
            S = set(beta[i] for v in clsmap.values() if len(v) >= 2 for i in v)
            if not S:
                continue
            C = min(S) + max(S)
            if set(C - b for b in S) != S:
                continue                                  # not concentric: no sigma_C to speak of
            E = set(k for k in clsmap if len(clsmap[k]) >= 2)
            fixed = set(k for k in range(t) if (2 * k - C) % t == 0)
            nonfixed_orbits = set()
            for k in range(t):
                if k not in fixed:
                    nonfixed_orbits.add(frozenset([k % t, (C - k) % t]))
            missed = sum(1 for o in nonfixed_orbits if not (o & E))
            v = True
            for q in range(3):
                if matrix(F, [[POW[q][b][beta[a]] for b in range(N)]
                              for a in range(N)]).det() != 0:
                    v = False
                    break
            if v:
                nz += 1
                if len(E) % 2:
                    odd += 1
                if missed:
                    miss += 1
            else:
                ctrl += 1
                if missed:
                    ctrl_miss += 1
    print("  %2d %3d %3d %6d %6d %8d %33d %10d/%d"
          % (t, r, N, MAX, nz, odd, miss, ctrl_miss, ctrl))

print("")
print("  C2: '|E| odd' must be 0 -- a sigma_C-stable set containing both fixed points has even size.")
print("  C3: 'zeros missing a non-fixed class' is the crossing.  If it is 0 the GL condition really")
print("      is 'occupy the fixed classes AND cover the non-fixed ones', the same shape as the Sp")
print("      lemma read the other way round.")
print("  C4: the control is the concentric NON-vanishing shapes.  If they miss classes just as")
print("      rarely, the crossing measures the concentric stratum and not the zeros.")
print("")
print("DONE")
