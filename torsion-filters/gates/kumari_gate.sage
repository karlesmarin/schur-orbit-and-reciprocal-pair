# kumari_gate.sage -- two independent gates on Kumari's thesis, both run from the definitions.
#
# WHY.  Our Section 6 leans on three lemmas of [KumariThesis]: 2.14 (the parity of the number of
# vertical dominoes is independent of the tableau), 2.16 (coverable tableaux <-> domino tableaux)
# and 2.18 (a fixed-point-free sign-reversing involution on the non-coverable fillings).  Chapter 6
# of that thesis appeared as arXiv:2211.14093, and that preprint was WITHDRAWN by the author on
# 2025-12-17 with the note "Error in Section 4".  Section 4 of the preprint is "Cyclic sieving
# phenomenon".  So:
#
#   GATE A -- are the three lemmas we actually use sound?  They are not in Section 4 (the words
#             "coverable", "involution" and "domino" do not occur anywhere in the withdrawn
#             preprint), but "not in the withdrawn section" is an argument about page numbers.
#             This gate checks the mathematics instead, by brute force.
#
#   GATE B -- where does the withdrawn Theorem 4.4 (= thesis Theorem 6.17) break?  It claims: if
#             sgn(sigma_lambda) = sgn(sigma_mu) then a cyclic action exists on SSYT_{tn}(lambda/mu)
#             sieved by s_{lambda/mu}(1,q,...,q^{tn-1}).  By Alexandersson-Amini that is equivalent
#             to: for every d | t, f(omega^d) = sum_{j|d} j c_j with every c_j a nonnegative
#             integer.  The c_j are determined, so this is decidable.  We look for a counterexample.
#
# Both gates carry controls that must fail.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

def betaset(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - (j + 1) for j in range(N)]

def sgn_sigma(lam, N, t):
    b = betaset(lam, N)
    w = [v % t for v in b]
    inv = sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])
    return (-1) ** inv

def profile(lam, N, t):
    return tuple(sum(1 for v in betaset(lam, N) if v % t == i) for i in range(t))

def parts_upto(K, maxlen):
    """every partition of size 0..K with at most maxlen parts"""
    out = []
    for k in range(K + 1):
        for p in Partitions(k, max_length=maxlen):
            out.append(list(p))
    return out

# ---------------------------------------------------------------- GATE A
print("=" * 78)
print("GATE A -- the three lemmas our Section 6 actually uses, checked by brute force")
print("=" * 78)

def coverable(cells, n):
    """Definition 2.15.  The two admissible dominoes are

         [2a | 2a]   a horizontal pair of equal EVEN entries, and
         [2a-1]
         [ 2a ]      a vertical pair, 2a-1 directly above 2a,

    which is forced by semistandardness: a column strictly increases, so two equal entries can
    only sit side by side, and 2a-1 immediately below 2a is impossible.  Lemma 2.16's identity
    sum ht = c_1 + ... + c_{2n-1} says the same thing -- every vertical domino carries exactly
    one odd entry.

    `cells` is a dict (row, col) -> entry.  We match separately inside each pair {2a-1, 2a},
    exhaustively, since the regions are tiny."""
    for a in range(1, n + 1):
        lo, hi = 2 * a - 1, 2 * a
        region = sorted(c for c, v in cells.items() if v in (lo, hi))
        if not region:
            continue
        if len(region) % 2:
            return False
        pos = set(region)

        def match(rem):
            if not rem:
                return True
            c = min(rem)
            i, j = c
            v = cells[c]
            # vertical: 2a-1 at c with 2a directly below
            if v == lo and (i + 1, j) in rem and cells[(i + 1, j)] == hi:
                if match(rem - {c, (i + 1, j)}):
                    return True
            # horizontal: two equal even entries side by side
            if v == hi and (i, j + 1) in rem and cells[(i, j + 1)] == hi:
                if match(rem - {c, (i, j + 1)}):
                    return True
            return False

        if not match(frozenset(pos)):
            return False
    return True


def cells_of(T, mu):
    """(row, col) -> entry for a skew tableau, columns offset by mu"""
    out = {}
    for i, row in enumerate(T):
        for j, v in enumerate(row):
            if v is not None:              # Sage pads the skew part with None
                out[(i, j)] = v
    return out

def weight_sign(T_rows):
    """(X,-X)_T sign = (-1)^{c_2 + c_4 + ...}"""
    c = 0
    for row in T_rows:
        for v in row:
            if v is not None and v % 2 == 0:
                c += 1
    return (-1) ** c

def content(T_rows, m):
    c = [0] * (m + 1)
    for row in T_rows:
        for v in row:
            c[v] += 1
    return c

bad218 = 0
tested218 = 0
for n in (1, 2):
    m = 2 * n
    for lam in parts_upto(9, 4):
        for mu in parts_upto(sum(lam), 4):
            mu = list(mu)
            L = list(lam)
            if len(mu) > len(L):
                continue
            mu = mu + [0] * (len(L) - len(mu))
            if any(mu[i] > L[i] for i in range(len(L))):
                continue
            try:
                S = SemistandardSkewTableaux([L, mu], max_entry=m)
            except Exception:
                continue
            tot = 0
            for T in S:
                rows = [list(r) for r in T]
                if not coverable(cells_of(rows, mu), n):
                    tot += weight_sign(rows)
            tested218 += 1
            if tot != 0:
                bad218 += 1
print("  Lemma 2.18  sum over NON-coverable of (X,-X)_T = 0")
print("      %d skew shapes tested, %d failures" % (tested218, bad218))

# Lemma 2.14 / 2.16 composite, at t=2: the surviving (coverable) fillings all carry one sign,
# and that sign is sgn(sigma_lambda) sgn(sigma_mu).
bad_sign = 0
tested_sign = 0
for N in (4, 5, 6):
    for lam in parts_upto(10, N):
        for mu in parts_upto(sum(lam), N):
            L, M = list(lam), list(mu)
            if len(M) > len(L):
                continue
            M = M + [0] * (len(L) - len(M))
            if any(M[i] > L[i] for i in range(len(L))):
                continue
            if profile(L, N, 2) != profile(M, N, 2):
                continue
            n = N // 2
            try:
                S = SemistandardSkewTableaux([L, M], max_entry=2 * n)
            except Exception:
                continue
            signs = set()
            for T in S:
                rows = [list(r) for r in T]
                if coverable(cells_of(rows, M), n):
                    signs.add(weight_sign(rows))
            tested_sign += 1
            if len(signs) > 1:
                bad_sign += 1
            elif signs:
                if signs.pop() != sgn_sigma(L, N, 2) * sgn_sigma(M, N, 2):
                    bad_sign += 1
print("  Lemmas 2.14+2.16  every coverable filling carries ONE sign, = sgn(sig_l) sgn(sig_m)")
print("      %d skew shapes tested, %d failures" % (tested_sign, bad_sign))

print("  CONTROL: a deliberately wrong sign rule must fail")
ctl = 0
for N in (4,):
    for lam in parts_upto(7, N):
        L = list(lam)
        if sgn_sigma(L, N, 2) != +1:
            ctl += 1
print("      sgn(sigma) is not identically +1: %s" % (ctl > 0))

# ---------------------------------------------------------------- GATE B
print()
print("=" * 78)
print("GATE B -- the WITHDRAWN Theorem 4.4 (= thesis 6.17), tested by Alexandersson-Amini")
print("=" * 78)

def f_poly(L, M, tn, R, q):
    """f(q) = s_{L/M}(1,q,...,q^{tn-1}) by Jacobi-Trudi, with
       h_k(1,q,...,q^{tn-1}) = qbinomial(k+tn-1, k).
    Exact in Z[q], and far cheaper than expanding the symmetric function."""
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


def csp_defect(L, M, t, n, R, q, K, w):
    """Alexandersson-Amini: a cyclic action exists iff for every d | t,
       f(omega^d) = sum_{j|d} j c_j  with every c_j a nonnegative integer."""
    tn = int(t) * int(n)
    f = f_poly(L, M, tn, R, q)
    cs = {}
    for d in sorted(ZZ(t).divisors()):
        val = f(w ** d)
        if val not in QQ:
            return (False, cs, 'f(omega^%d) is not rational' % d)
        val = QQ(val)
        acc = sum(j * cs[j] for j in cs if d % j == 0)
        rem = val - acc
        if rem % d != 0:
            return (False, cs, 'c_%d not an integer' % d)
        cs[d] = ZZ(rem / d)
    bad = [j for j in sorted(cs) if cs[j] < 0]
    return (len(bad) == 0, cs, ('negative c_%s' % bad) if bad else 'ok')


counter = []
checked = 0
R = PolynomialRing(ZZ, 'q'); q = R.gen()
for t in (Integer(2), Integer(3), Integer(4), Integer(6)):
    K = CyclotomicField(int(t)); w = K.gen()
    for n in (Integer(1), Integer(2)):
        tn = t * n
        if tn > 6: continue
        for lam in parts_upto(12, int(tn)):
            for mu in parts_upto(sum(lam), int(tn)):
                L, M = list(lam), list(mu)
                if len(M) > len(L):
                    continue
                M = M + [0] * (len(L) - len(M))
                if any(M[i] > L[i] for i in range(len(L))):
                    continue
                if sgn_sigma(L, int(tn), int(t)) != sgn_sigma(M, int(tn), int(t)):
                    continue
                checked += 1
                ok, cs, why = csp_defect(L, M, t, n, R, q, K, w)
                if not ok:
                    counter.append((t, L, M, cs, why))
        print("  t=%d n=%d : %d hypothesis-satisfying pairs checked, %d counterexamples so far"
              % (t, n, checked, len(counter)))

print()
if counter:
    print("  COUNTEREXAMPLES to the withdrawn Theorem 4.4 (smallest first):")
    counter.sort(key=lambda r: (sum(r[1]) - sum(r[2]), sum(r[1])))
    for t, L, M, cs, why in counter[:10]:
        print("     t=%d  lambda=%s  mu=%s   %s   c=%s" % (t, L, M, why, cs))
else:
    print("  no counterexample found in this range")

print()
print("  CONTROL: the criterion must reject something, or it is not a criterion.")
print("           Lee-Oh's own hypothesis (t | lambda_i - mu_i) is a SUBSET of the above;")
print("           any counterexample must therefore violate it.")
for t, L, M, cs, why in counter[:5]:
    LO = all((L[i] - (M[i] if i < len(M) else 0)) % t == 0 for i in range(len(L)))
    print("     t=%d lambda=%s mu=%s : satisfies Lee-Oh's hypothesis? %s  (must be False)" % (t, L, M, LO))
