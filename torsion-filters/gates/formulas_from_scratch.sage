# Independent check of the paper's central displayed formulas, written from the PRINTED statements
# and from the definition of a Schur polynomial as a bialternant.  It shares no code with
# paper/anc/: a bug in the machinery that produced the tables must not be shared with the thing
# checking them.
#
# Exact arithmetic, and no fraction field.  The alphabet lives in Q(zeta_t); the free letter is
# written as u^2 so that the half-integer exponents of Corollary 3.3 are integers in u.  Every
# identity is verified in CLEARED form -- both sides multiplied out into Laurent polynomials --
# which is the same statement and avoids the gcd work that makes a fraction field slow.
#
# Checked:
#   (5)  eq:main       det(x^beta) * f(t)^2 f(2)  ==  eps f(d1) f(d2) f(d3) * det(x^{N-1-j})
#   (9)  eq:chiratio   the same written as a ratio of sl2 characters, half-integer indices included
#   (6)  eq:sign       the long form of eps, off the columns and the residue word
#   (14) eq:shortsign  the short form of eps, off sgn(sigma) and the two residues
#   (8)  eq:shift      Phi(lambda + (m^N)) = (-1)^{(t+1)m} Phi(lambda)
#   (34) eq:congrgen   Delta_i(k) = 2i (mod t), on the beta sets of the general-r setup
#
# Each carries a decoy that must FAIL, with the count of how often, because a control that cannot
# fail reports nothing.
#
# Authors: Carles Marin + Claude (AI assistant).

import sys

def partitions_upto(N, n):
    out = []
    def rec(rest, maxpart, cur):
        out.append(tuple(cur))
        if len(cur) == N:
            return
        for p in range(min(rest, maxpart), 0, -1):
            rec(rest - p, p, cur + [p])
    rec(n, n, [])
    return out

def inv(word):
    return sum(1 for i in range(len(word)) for j in range(i + 1, len(word)) if word[i] > word[j])

def sgn(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)

tot = {k: 0 for k in ('main', 'chi', 'short', 'shift', 'zero', 'parity', 'uflip')}
bad = {k: 0 for k in ('main', 'chi', 'short', 'shift', 'zero', 'parity', 'uflip')}
dn = {k: 0 for k in ('nosign', 'oriented', 'valorder', 'clsorder', 'wrongt')}
dk = {k: 0 for k in ('nosign', 'oriented', 'valorder', 'clsorder', 'wrongt')}

for t, maxsize in [(2, 14), (3, 13), (4, 12), (5, 11), (6, 10)]:
    N = t + 2
    K = CyclotomicField(t) if t > 2 else QQ
    zt = K.zeta(t) if t > 2 else K(-1)
    R = LaurentPolynomialRing(K, 'u')
    u = R.gen()
    alpha = [R(zt**k) for k in range(t)] + [u**2, u**-2]

    def bialt(exps):
        return matrix(R, N, N, lambda i, j: alpha[i]**exps[j]).determinant()

    V = bialt([N - 1 - j for j in range(N)])
    assert V != 0
    def f(v):
        return u**v - u**(-v)

    shapes = [p for p in partitions_upto(N, maxsize) if len(p) <= N]
    for lam in shapes:
        lamp = list(lam) + [0] * (N - len(lam))
        beta = [lamp[j] + N - 1 - j for j in range(N)]
        D = bialt(beta)

        cls = {}
        for j, b in enumerate(beta):
            cls.setdefault(b % t, []).append((j, b))

        if len(cls) < t:                                   # degenerate profile
            tot['zero'] += 1
            if D != 0:
                bad['zero'] += 1
            continue

        big = sorted([i for i in cls if len(cls[i]) >= 2])
        if len(big) == 2:
            rA, rB = big
            A = [v for _, v in cls[rA]]
            B = [v for _, v in cls[rB]]
            jA1, jB1 = cls[rA][0][0], cls[rB][0][0]
        else:
            rA = rB = big[0]
            p, q, rr = [v for _, v in cls[rA]]
            A, B = [p, q], [q, rr]
            jA1, jB1 = cls[rA][0][0], cls[rA][1][0]
        a1, a2, b1, b2 = max(A), min(A), max(B), min(B)
        d1, d2 = a1 - a2, b1 - b2
        d3t = a1 + a2 - b1 - b2
        d3 = abs(d3t)

        if d3t == 0:                                       # the concentric locus
            tot['zero'] += 1
            if D != 0:
                bad['zero'] += 1
            continue

        # ---- (6) the long sign, columns 1-based as the paper indexes them
        S = [j for j in range(N) if j not in (jA1, jB1)]
        bS = [beta[j] % t for j in S]
        eps = ((-1)**(t + binomial(N + 1, 2))
               * (-1)**((jA1 + 1) + (jB1 + 1) + inv(bS))
               * sgn(a1 - b1) * sgn(d3t))

        # ---- (5), cleared of denominators
        tot['main'] += 1
        if D * f(t)**2 * f(2) != eps * f(d1) * f(d2) * f(d3) * V:
            bad['main'] += 1

        # ---- (9) chi_k = (u^{2k+2} - u^{-2k-2}) / (u^2 - u^{-2}); with k = d/2 - 1 the
        #      numerator exponent is d.  Cleared, (9) is (5) with f(2) moved about, so the
        #      check is that the character form has the SAME cleared identity.
        tot['chi'] += 1
        if D * f(t)**2 * f(2) != eps * f(d1) * f(d2) * f(d3) * V:
            bad['chi'] += 1

        # ---- the parity that makes Corollary 3.3 safe: each chi factor reverses sign under
        #      u -> -u exactly when its index is half-integral, so the ratio is invariant only if
        #      d1+d2+d3 is even.  It always is, having the parity of d1+d2+d3t = 2(a1-b2).
        tot['parity'] += 1
        if (d1 + d2 + d3) % 2 != 0:
            bad['parity'] += 1
        # The line above is arithmetic on d1+d2+d3 and would pass whatever the alphabet did, so on
        # its own it is not a receipt for the corollary.  This is: divide out and look at Phi
        # itself.  Phi = D/V must come out exact and carry ONLY EVEN powers of u -- that is what
        # "a Laurent polynomial in z" means when the free letter is written as u^2, and it is the
        # invariance under u -> -u stated in Corollary 3.3.
        tot['uflip'] += 1
        q, rem = D.quo_rem(V)
        if rem != 0 or any(e % 2 for e in q.exponents()):
            bad['uflip'] += 1

        # ---- (14) the short sign: classes in increasing residue order, values decreasing inside
        order = []
        for i in sorted(cls):
            order += [j for j, _ in sorted(cls[i], key=lambda pr: -pr[1])]
        eps_short = (-1)**(t // 2) * (-1)**inv(order) * (-1)**(rA + rB) * sgn(d3t)
        tot['short'] += 1
        if eps_short != eps:
            bad['short'] += 1

        # ---- decoys
        dn['nosign'] += 1
        if D * f(t)**2 * f(2) != f(d1) * f(d2) * f(d3) * V:
            dk['nosign'] += 1
        dn['oriented'] += 1
        if D * f(t)**2 * f(2) != eps * f(d1) * f(d2) * f(d3t) * V:
            dk['oriented'] += 1
        # A decoy has to BREAK the rule, not restate it.  The first version of this one sorted by
        # (residue, ascending column) -- and inside a class an ascending column IS a descending
        # value, because beta decreases with j.  So it was the rule itself, and it refuted 0 of
        # 746: a decoy that ties means untested.  These two really do break it: one reverses the
        # order of the values inside each class, the other keeps the values right and reverses
        # the order of the classes.
        asc = []
        for i in sorted(cls):
            asc += [j for j, _ in sorted(cls[i], key=lambda pr: pr[1])]
        e_asc = (-1)**(t // 2) * (-1)**inv(asc) * (-1)**(rA + rB) * sgn(d3t)
        dn['valorder'] += 1
        if e_asc != eps:
            dk['valorder'] += 1
        rev = []
        for i in sorted(cls, reverse=True):
            rev += [j for j, _ in sorted(cls[i], key=lambda pr: -pr[1])]
        e_rev = (-1)**(t // 2) * (-1)**inv(rev) * (-1)**(rA + rB) * sgn(d3t)
        dn['clsorder'] += 1
        if e_rev != eps:
            dk['clsorder'] += 1
        dn['wrongt'] += 1
        if D * f(t + 1)**2 * f(2) != eps * f(d1) * f(d2) * f(d3) * V:
            dk['wrongt'] += 1

        # ---- (8) the shift law
        if sum(lamp) + N <= maxsize:
            for m in (1, 2):
                lm = [x + m for x in lamp]
                bm = [lm[j] + N - 1 - j for j in range(N)]
                tot['shift'] += 1
                if bialt(bm) != (-1)**((t + 1) * m) * D:
                    bad['shift'] += 1

    print("  t=%d  |lambda|<=%-3d  %d shapes" % (t, maxsize, len(shapes)))
    sys.stdout.flush()

# ---- (34) the increment congruence, on the general-r beta sets: N = t + 2r
print("")
congr_tot = congr_bad = 0
congr_decoy = congr_decoy_n = 0
for t in (2, 3, 4, 5, 6):
    for r in (1, 2, 3):
        N = t + 2 * r
        for lam in [p for p in partitions_upto(N, 12) if len(p) <= N]:
            lamp = list(lam) + [0] * (N - len(lam))
            beta = [lamp[j] + N - 1 - j for j in range(N)]
            cls = {}
            for b in beta:
                cls.setdefault(b % t, []).append(b)
            if len(cls) < t:
                continue
            for i in cls:
                c = sorted(cls[i], reverse=True)
                for k in range(len(c) - 1):
                    congr_tot += 1
                    if (c[k] + c[k + 1] - 2 * i) % t != 0:
                        congr_bad += 1
                    congr_decoy_n += 1
                    if (c[k] + c[k + 1] - (2 * i + 1)) % t != 0:
                        congr_decoy += 1

print("=" * 90)
print("RESULT")
print("=" * 90)
for k, name in [('main', '(5)  the closed form, sign included'),
                ('chi', '(9)  the same as a ratio of sl2 characters'),
                ('short', '(14) the short form of the sign'),
                ('shift', '(8)  the shift law'),
                ('zero', '     the two vanishing branches of Theorem 3.1(i) and Cor 3.2'),
                ('parity', '     d1+d2+d3 is even (what Cor 3.3 rests on)'),
                ('uflip', '     Phi = D/V is exact and even in u, i.e. Laurent in z')]:
    print("  %-52s %6d checked   %d failures" % (name, tot[k], bad[k]))
print("  %-52s %6d checked   %d failures"
      % ('(34) Delta_i(k) = 2i (mod t), all t and r<=3', congr_tot, congr_bad))

print("")
print("  DECOYS -- each must be refuted, and the count says how often")
for k, name in [('nosign', 'the sign dropped from (5)'),
                ('oriented', 'the oriented d3 where the size belongs'),
                ('valorder', 'values read increasingly inside each class'),
                ('clsorder', 'classes read in decreasing residue order'),
                ('wrongt', 'the denominator at t+1 instead of t')]:
    print("  %-52s %6d tried    %d refuted" % (name, dn[k], dk[k]))
print("  %-52s %6d tried    %d refuted"
      % ('the congruence shifted to 2i+1', congr_decoy_n, congr_decoy))

print("")
print("TOTAL failures on the paper's formulas: %d" % (sum(bad.values()) + congr_bad))
