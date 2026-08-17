# -*- coding: utf-8 -*-
"""Two statements about what the evaluation invariant sees, and what it does not.

Proposition 3.10 says the invariant is MINIMAL: at a fixed t, two partitions take the same nonzero
value if and only if they agree on the multiset {d_1,d_2,d_3} and on the sign.  The paper proves it
by unique factorization in Z[u]; this script tests the conclusion directly, and over a wider
population than the partitions realising it, because the map to check for injectivity is

        (multiset, eps)  |-->  eps * prod_i (u^{d_i} - u^{-d_i}),      u = z^(1/2),

the denominator of (2) depending on t alone.  That is exact integer polynomial arithmetic, so no
tolerance enters anywhere.  The population is every triple with entries up to a bound, whether or
not some partition realises it: a collision there would refute the proposition a fortiori.

Identity (10) is the shift law.  A common shift of the beta set is lambda -> lambda + (m^N), and
s_{lambda+(m^N)}(A) = det(A)^m s_lambda(A), so the value picks up (-1)^{(t+1)m} while the interval
triple does not move at all.  For odd t that is invisible; for even t an odd shift reverses the
value, and since the triple is fixed the whole reversal has to be carried by eps_lambda.  That
makes (10) a free check on the sign formula (7), and it is the reason the three-factor slogan of
section 3 is stated about the TRIPLE and not about the value: an earlier draft said a common shift
was invisible to (2), which is false for even t.

CONTROL.  A test of injectivity can pass by measuring nothing, so the run also feeds the same
machinery a deliberately lossy invariant -- the SUM d_1+d_2+d_3 with the sign -- which must
collide, and prints how often.  And the shift law is re-run with the exponent (-1)^{t+1} replaced
by +1, which must fail at every even t.

Authors: Carles Marin, Claude (AI assistant)."""
import sys
from itertools import combinations_with_replacement

from mpmath import mp, mpf, mpc, exp, pi

mp.dps = 40
THETA = mpf("0.41")
DBOUNDS = (12, 30, 60)


# ------------------------------------------------------------------ the invariant, from (2)-(7)

def partitions(n, k):
    if n == 0:
        yield ()
        return
    if k == 0:
        return
    for first in range(n, 0, -1):
        for rest in partitions(n - first, k - 1):
            if not rest or rest[0] <= first:
                yield (first,) + rest


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - 1 - j for j in range(N)]


def inv_count(word):
    return sum(1 for i in range(len(word)) for j in range(i + 1, len(word)) if word[i] > word[j])


def invariant(lam, t):
    """(d_1, d_2, d_3, eps) by (2) and (7); None when Phi_t vanishes identically."""
    N = t + 2
    beta = beta_of(lam, N)
    cls = {}
    for j, b in enumerate(beta):
        cls.setdefault(b % t, []).append(j)
    if len(cls) < t:
        return None
    big = [v for v in cls.values() if len(v) >= 2]
    if len(big) == 1:
        j1, j2, j3 = big[0]
        (jA1, jA2), (jB1, jB2) = (j1, j2), (j2, j3)
    else:
        big.sort()
        (jA1, jA2), (jB1, jB2) = tuple(big[0]), tuple(big[1])
    a1, a2, b1, b2 = beta[jA1], beta[jA2], beta[jB1], beta[jB2]
    c = a1 + a2 - b1 - b2                       # the ORIENTED coupling of Theorem 3.1
    if c == 0:
        return None
    S = [j for j in range(N) if j not in (jA1, jB1)]
    k11 = (-1) ** ((jA1 + 1) + (jB1 + 1) + inv_count([beta[j] % t for j in S]))
    k11 *= 1 if a1 > b1 else -1
    eps = (-1) ** (t + (N + 1) * N // 2) * k11 * (1 if c > 0 else -1)
    return (a1 - a2, b1 - b2, abs(c), eps)


def mydet(M, n):
    M = [row[:] for row in M]
    d = mpc(1)
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        if abs(M[p][k]) == 0:
            return mpc(0)
        if p != k:
            M[k], M[p] = M[p], M[k]
            d = -d
        d *= M[k][k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    return d


def phi_bialt(lam, t, theta=THETA):
    """the object itself, from the bialternant; used only for the shift law"""
    N = t + 2
    z = exp(theta)
    w = exp(2j * pi / t)
    xs = [w ** k for k in range(t)] + [z, 1 / z]
    beta = beta_of(lam, N)
    num = [[xs[i] ** beta[j] for j in range(N)] for i in range(N)]
    den = [[xs[i] ** (N - 1 - j) for j in range(N)] for i in range(N)]
    return mydet(num, N) / mydet(den, N)


# ------------------------------------------------------- exact Laurent arithmetic for the value

def lau_mul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            out[a + b] = out.get(a + b, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def numerator(d, eps):
    """eps * prod_i (u^{d_i} - u^{-d_i}), exactly, as a Laurent polynomial in u = z^(1/2)"""
    p = {0: eps}
    for x in d:
        p = lau_mul(p, {x: 1, -x: -1})
    return tuple(sorted(p.items()))


def lossy(d, eps):
    """the CONTROL invariant: the SUM instead of the multiset.  Must collide."""
    return (sum(d), eps)


# ------------------------------------------------------------------------------------ the runs

def run_minimality():
    print("=" * 78)
    print("MINIMALITY: distinct (multiset, sign) must give distinct values")
    print("=" * 78)
    bad = 0
    for D in DBOUNDS:
        seen, coll = {}, 0
        seen2, coll2 = {}, 0
        n = 0
        for d in combinations_with_replacement(range(1, D + 1), 3):
            for eps in (1, -1):
                n += 1
                k = numerator(d, eps)
                if k in seen and seen[k] != (d, eps):
                    coll += 1
                seen[k] = (d, eps)
                k2 = lossy(d, eps)
                if k2 in seen2 and seen2[k2] != (d, eps):
                    coll2 += 1
                seen2[k2] = (d, eps)
        bad += coll
        print("  d_i <= %2d : %6d invariants, %6d distinct values, collisions %d"
              "     [control, graded by the SUM: %d collisions]" % (D, n, len(seen), coll, coll2))
    print()
    print("  and over the partitions that realise them:")
    for t in range(2, 6):
        N, vals, n = t + 2, {}, 0
        for size in range(0, 17):
            for lam in partitions(size, N):
                iv = invariant(lam, t)
                if iv is None:
                    continue
                n += 1
                vals.setdefault(numerator(sorted(iv[:3]), iv[3]), set()).add(
                    (tuple(sorted(iv[:3])), iv[3]))
        wide = sum(1 for v in vals.values() if len(v) > 1)
        bad += wide
        print("  t=%d, |lambda| <= 16 : %4d shapes, %3d distinct values, "
              "values carrying two invariants %d" % (t, n, len(vals), wide))
    return bad


def run_shift():
    print()
    print("=" * 78)
    print("THE SHIFT LAW (10):  Phi_t(lambda + (1^N)) = (-1)^(t+1) Phi_t(lambda)")
    print("=" * 78)
    print("   t   triple fixed   eps carries (-1)^(t+1)   value law   [control: sign law without")
    print("                                                            the exponent]")
    bad = 0
    tot = tot_ok = 0
    for t in range(2, 8):
        N = t + 2
        n = ok_tr = ok_eps = ok_val = ctl = 0
        for size in range(0, 13):
            for lam in partitions(size, N):
                iv = invariant(lam, t)
                lam2 = tuple(x + 1 for x in (list(lam) + [0] * (N - len(lam))))
                iv2 = invariant(lam2, t)
                if (iv is None) != (iv2 is None):
                    bad += 1
                    continue
                if iv is None:
                    continue
                n += 1
                ok_tr += (iv[:3] == iv2[:3])
                ok_eps += (iv2[3] == iv[3] * (-1) ** (t + 1))
                ok_val += bool(abs(phi_bialt(lam2, t) - (-1) ** (t + 1) * phi_bialt(lam, t))
                               < mpf(10) ** -20)
                ctl += (iv2[3] == iv[3])              # the control: no exponent at all
        bad += (n - ok_tr) + (n - ok_eps) + (n - ok_val)
        tot += n
        tot_ok += min(ok_tr, ok_eps, ok_val)
        print("   %d      %4d/%4d       %4d/%4d              %4d/%4d     %4d/%4d %s"
              % (t, ok_tr, n, ok_eps, n, ok_val, n, ctl, n,
                 "(must be n at odd t, 0 at even)" if t == 2 else ""))
    print()
    print("  the shift law over t <= 7 and |lambda| <= 12: %d shapes, %d fail"
          % (tot, tot - tot_ok))
    print()
    print("  The triple is translation-invariant; the value is not, for even t.  The control column")
    print("  is the same test with the exponent dropped: it agrees at odd t, where the law is")
    print("  trivial, and disagrees everywhere at even t, which is where the law has content.")
    return bad


if __name__ == "__main__":
    rc = run_minimality() + run_shift()
    print()
    print("TOTAL FAILURES: %d" % rc)
    sys.exit(1 if rc else 0)
