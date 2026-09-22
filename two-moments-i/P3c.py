#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""P3c.py --- every count quoted as measured in the note, and a self-test against the
published values.

    Subsets of F_p with two prescribed moments, Part I: the pair {1,3}.
    Carles Marin, 2026.   Written with Claude (Anthropic) as an assistant.

WHAT THIS PROGRAM IS FOR.  The note states four closed forms, two tables of strata, three
propositions about individual surfaces, and one classification.  Every one of those is a number
that can be recomputed from scratch, and this file recomputes all of them and checks them
against what the note prints.  Running it is the whole verification:

    python P3c.py

It needs only the standard library and numpy.  It prints one line per claim and exits non-zero
if any of them fails, so it can be wired into a build without being read.

WHAT IT COMPUTES, in the order the note uses it.

  C(n, p)              the census itself, by dynamic programming over (size, sum, sum of cubes).
                       This is the definition, and it is the thing everything else must match.
  N_lambda(lam, p)     the affine count of the stratum lambda, by brute force over F_p^(l-1).
  sieve(n, p)          the Li--Wan partition sieve, eq. (2) of the note, assembled from those
                       N_lambda.  Three independent routes to the same number: census, sieve,
                       closed form.
  sieve_dp(n, p)       the same sieve with each N_lambda computed by dynamic programming over
                       the two sums (N_lambda_dp), a route that never looks at V_lambda.
  singular_locus(lam)  Lemma 2 read as Corollary 4: the sign vectors that split the weights into
                       two blocks of equal sum.  Computed by subset-sum DP, not by 2^l search.
  classification(nmax) Theorem 6: C(n,p) is elementary exactly when every lambda of length four
                       is biquanimous.

A NOTE ON WHAT IS *NOT* HERE.  The symbolic steps --- the pencil identities, the change of
coordinates to Cayley's standard form --- are not redone symbolically here: the counts they lead
to are checked over F_p, and the identity F(Mu) = 24(...) at every point of F_p^4.  Their
symbolic verification is in the scripts under gates/ that the note names.
"""
from __future__ import print_function

import sys
from fractions import Fraction
from itertools import combinations, product
from math import factorial

try:
    import numpy as np
except ImportError:                                            # pragma: no cover
    sys.exit("P3c.py needs numpy:  pip install numpy")

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:                                              # pragma: no cover
    pass


# =====================================================================================
#  0.  Small arithmetic
# =====================================================================================
def chi(a, p):
    u"""The quadratic character of a mod p, with chi(0) = 0."""
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def partitions(n, largest=None):
    u"""The partitions of n, as weakly decreasing tuples."""
    largest = largest or n
    if n == 0:
        yield ()
        return
    for w in range(min(n, largest), 0, -1):
        for rest in partitions(n - w, w):
            yield (w,) + rest


def z_lambda(lam):
    u"""z_lambda = prod_j j^{m_j} m_j!, the centraliser order of the cycle type."""
    z = 1
    for w in set(lam):
        m = lam.count(w)
        z *= (w ** m) * factorial(m)
    return z


# =====================================================================================
#  1.  The census C(n, p), by dynamic programming
# =====================================================================================
def C(n, p):
    u"""#{T subset F_p : |T| = n, sum_{t in T} t = sum_{t in T} t^3 = 0}.

    Dynamic programming over (size, sum, sum of cubes): each element of F_p is offered once,
    so subsets are counted without order and without repetition.  O(n p^3) states.
    """
    dp = [dict() for _ in range(n + 1)]
    dp[0][(0, 0)] = 1
    for t in range(p):
        c3 = (t * t * t) % p
        for k in range(min(n - 1, t), -1, -1):
            for (s, c), v in list(dp[k].items()):
                key = ((s + t) % p, (c + c3) % p)
                dp[k + 1][key] = dp[k + 1].get(key, 0) + v
    return dp[n].get((0, 0), 0)


def C_bruteforce(n, p):
    u"""The same count by direct enumeration.  Only for small p: it is the control on C()."""
    return sum(1 for T in combinations(range(p), n)
               if sum(T) % p == 0 and sum(t ** 3 for t in T) % p == 0)


# =====================================================================================
#  2.  The strata, and the sieve
# =====================================================================================
def V_points(lam, p):
    u"""#V_lambda(F_p), the projective count, by brute force.

    The last variable is eliminated with the linear equation; the loop runs over the first
    coordinate so the live array is p^(l-2) and not p^(l-1).
    """
    l = len(lam)
    w = [x % p for x in lam]
    # Only the weight being eliminated has to be invertible.  A *different* weight vanishing
    # is legitimate and interesting: the variety degenerates, the count is still well defined,
    # and Corollary 4 stops applying --- which is exactly the (5,1,1,1) at p = 5 that the note
    # uses to show the hypothesis p > n cannot be dropped on its own.
    if w[-1] == 0:
        raise ValueError("p divides the eliminated weight of %s" % (lam,))
    if l == 1:
        return 0
    if l == 2:
        n = sum(1 for a in range(p) for b in range(p)
                if (a or b) and (w[0] * a + w[1] * b) % p == 0
                and (w[0] * pow(a, 3, p) + w[1] * pow(b, 3, p)) % p == 0)
        return n // (p - 1)
    inv = pow(w[-1], p - 2, p)
    cub = (np.arange(p, dtype=np.int64) ** 3) % p
    g = np.arange(p, dtype=np.int64)
    d = l - 2
    base = np.zeros((p,) * d, dtype=np.int64)
    lin = np.zeros((p,) * d, dtype=np.int64)
    for i in range(d):
        f = g.reshape((1,) * i + (p,) + (1,) * (d - 1 - i))
        base = (base + w[i + 1] * cub[f]) % p
        lin = (lin + w[i + 1] * f) % p
    zeros = 0
    for a in range(p):
        xl = (-(w[0] * a + lin) * inv) % p
        F = (w[0] * cub[a] + base + w[-1] * cub[xl]) % p
        zeros += int(np.count_nonzero(F == 0))
    return (zeros - 1) // (p - 1)


def N_lambda(lam, p):
    u"""#A_lambda(F_p), the affine cone: 1 + (p-1) #V_lambda."""
    return 1 + (p - 1) * V_points(lam, p)


def sieve(n, p):
    u"""C(n,p) by the Li--Wan partition sieve, eq. (2) of the note."""
    total = Fraction(0)
    for lam in partitions(n):
        total += Fraction((-1) ** (n - len(lam)), z_lambda(lam)) * N_lambda(lam, p)
    assert total.denominator == 1, "the sieve did not return an integer at n=%d, p=%d" % (n, p)
    return int(total)


def N_lambda_dp(lam, p):
    u"""N_lambda by a second, independent route: dynamic programming over the pair
    (sum w_i x_i, sum w_i x_i^3) in F_p x F_p, one coordinate at a time.  O(l p^3), and it
    never looks at V_lambda, so it does not share a line of code with V_points()."""
    state = np.zeros((p, p), dtype=np.int64)
    state[0, 0] = 1
    cubes = [pow(x, 3, p) for x in range(p)]
    for w in lam:
        new = np.zeros((p, p), dtype=np.int64)
        for x in range(p):
            new += np.roll(np.roll(state, (w * x) % p, axis=0), (w * cubes[x]) % p, axis=1)
        state = new
    return int(state[0, 0])


def sieve_dp(n, p):
    u"""The sieve, eq. (2), with every N_lambda taken from N_lambda_dp()."""
    return sum(Fraction((-1) ** (n - len(lam)), z_lambda(lam)) * N_lambda_dp(lam, p)
               for lam in partitions(n))


# =====================================================================================
#  3.  The singular locus: Lemma 2, read as Corollary 4
# =====================================================================================
def is_biquanimous(lam):
    u"""Can the multiset of weights be split into two blocks of equal total?"""
    n = sum(lam)
    if n % 2:
        return False
    reach = 1                                    # bitmask of reachable partial sums
    for w in lam:
        reach |= reach << w
    half = n // 2
    if not (reach >> half) & 1:
        return False
    # a split must be proper: the empty block and the whole multiset do not count
    return any(sum(lam[i] for i in A) == half
               for k in range(1, len(lam)) for A in combinations(range(len(lam)), k))


def singular_locus(lam):
    u"""#Sing(V_lambda) for l >= 4 and p > n: half the number of equal-sum splittings.

    Corollary 4 of the note.  Counted by subset-sum dynamic programming over the partial sums,
    which is O(l n) and not O(2^l).
    """
    n = sum(lam)
    if n % 2:
        return 0
    count = [0] * (n + 1)
    count[0] = 1
    for w in lam:
        for s in range(n, w - 1, -1):
            count[s] += count[s - w]
    return count[n // 2] // 2


def singular_bruteforce(lam):
    u"""The same, by direct search over sign vectors.  The control on singular_locus()."""
    n, l = sum(lam), len(lam)
    return sum(1 for eps in product((1, -1), repeat=l)
               if eps[0] == 1 and 2 * sum(w for w, e in zip(lam, eps) if e == 1) == n)


def classification(nmax):
    u"""Theorem 6: the n for which every lambda |- n of length four is biquanimous."""
    out = []
    for n in range(3, nmax + 1):
        four = [lam for lam in partitions(n) if len(lam) == 4]
        if all(is_biquanimous(lam) for lam in four):
            out.append(n)
    return out


# =====================================================================================
#  4.  The published closed forms
# =====================================================================================
def a_p(p):
    u"""The Frobenius trace of E: y^2 = x^3 + x^2 - x, of conductor 20."""
    count = 1
    for x in range(p):
        count += 1 + chi((x * x * x + x * x - x) % p, p)
    return p + 1 - count


def C3_closed(p):
    return Fraction(p - 1, 2)


def C4_closed(p):
    return Fraction((p - 1) * (p - 3), 8)


def C5_closed(p):
    return Fraction(p - 1, 120) * (p * p - 4 * p + 61 + chi(5, p) * (p + 15)
                                   + 20 * chi(-15, p) + 10 * a_p(p))


def C6_closed(p):
    return Fraction(p - 1, 720) * (p ** 3 - 9 * p * p + 81 * p - 399
                                   - 90 * chi(-1, p) - 40 * chi(-3, p))


CLOSED = {3: C3_closed, 4: C4_closed, 5: C5_closed, 6: C6_closed}


# =====================================================================================
#  5.  The three propositions about individual surfaces
# =====================================================================================
P22111 = (1, 8, 16, 10)                          # r^3 + 8r^2 + 16r + 10


def P_roots(p):
    u"""The distinct roots in F_p of P(r) = r^3 + 8r^2 + 16r + 10."""
    return [r for r in range(p)
            if (r * r * r + 8 * r * r + 16 * r + 10) % p == 0]


def V22111_closed(p):
    u"""Proposition 9: the count of V_(2,2,1,1,1), valid for p > 3."""
    return p * p + 1 + p * (3 + chi(5, p) + sum(chi(-2 * r, p) for r in P_roots(p)))


def V22111_collapsed(p):
    u"""Corollary 10: the same, for p > 5, with the character sum collapsed."""
    return p * p + 1 + p * (3 + chi(5, p) * (1 + len(P_roots(p))))


def Va1111_closed(a, p):
    u"""Proposition 11: the abelian family, for p not dividing 6a(a^2-4)(a^2-16)."""
    D = -3 * (a * a - 4)
    E = (a * a - 4) * (a * a - 16)
    return p * p + (3 + 3 * chi(D, p) + chi(E, p)) * p + 1


CAYLEY_NODES = [(1, -1, -1, -1), (1, -1, -1, 1), (1, -1, 1, -1), (1, 1, -1, -1)]


def cayley_identity_mod_p(p):
    u"""Proposition 16(2) over F_p: F(Mu) = 24 (u0u1u2 + u0u1u3 + u0u2u3 + u1u2u3).

    Checked at every point of F_p^4, which for p > 3 is more than a polynomial identity of
    degree three needs.
    """
    M = [[CAYLEY_NODES[j][i] for j in range(4)] for i in range(4)]
    g = np.arange(p, dtype=np.int64)
    U = [g.reshape((1,) * i + (p,) + (1,) * (3 - i)) for i in range(4)]
    X = [sum(M[i][j] * U[j] for j in range(4)) % p for i in range(4)]
    x5 = (-(2 * X[0] + X[1] + X[2] + X[3])) % p
    F = (2 * pow(X[0], 3) + pow(X[1], 3) + pow(X[2], 3) + pow(X[3], 3)
         + pow(x5, 3)) % p
    cay = (U[0] * U[1] * U[2] + U[0] * U[1] * U[3]
           + U[0] * U[2] * U[3] + U[1] * U[2] * U[3]) % p
    return bool(np.all((F - 24 * cay) % p == 0))


# =====================================================================================
#  6.  The self-test
# =====================================================================================
_ok = _bad = 0


def check(label, cond, detail=""):
    global _ok, _bad
    if cond:
        _ok += 1
        print(u"  ok   %s%s" % (label, ("   " + detail) if detail else ""))
    else:
        _bad += 1
        print(u"  FAIL %s%s" % (label, ("   " + detail) if detail else ""))


def self_test(primes=(7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)):
    print(u"P3c.py --- self-test against the published values")
    print(u"=" * 74)

    # -- the census, by two routes -------------------------------------------------
    bad = [(n, p) for n in (3, 4, 5) for p in (7, 11, 13)
           if C(n, p) != C_bruteforce(n, p)]
    check(u"the census C(n,p): dynamic programming agrees with direct enumeration",
          not bad, u"checked n = 3,4,5 and p = 7,11,13")

    # -- the four closed forms -----------------------------------------------------
    for n in (3, 4, 5, 6):
        bad = [(p, C(n, p), CLOSED[n](p)) for p in primes
               if p > n and Fraction(C(n, p)) != CLOSED[n](p)]
        check(u"closed form for C(%d,p) reproduces the census" % n, not bad,
              u"%d primes" % len([p for p in primes if p > n]))

    # -- the sieve -----------------------------------------------------------------
    for n in (5, 6):
        bad = [p for p in primes if p > n and sieve(n, p) != C(n, p)]
        check(u"the sieve, assembled from the strata, reproduces C(%d,p)" % n, not bad)
    odd = [p for p in range(3, 54) if all(p % q for q in range(2, p))]
    cases = [(n, p) for n in range(3, 8) for p in odd]
    bad = [(n, p) for (n, p) in cases if sieve_dp(n, p) != C(n, p)]
    check(u"the sieve with N_lambda by DP reproduces the census, 3 <= n <= 7, odd p <= 53",
          not bad, u"%d cases" % len(cases))

    # -- the singular locus --------------------------------------------------------
    lams = [lam for n in range(4, 11) for lam in partitions(n) if len(lam) >= 4]
    bad = [lam for lam in lams if singular_locus(lam) != singular_bruteforce(lam)]
    check(u"the singular locus by DP agrees with the sign-vector search",
          not bad, u"%d partitions with l >= 4, n <= 10" % len(lams))
    bad = [lam for lam in lams if (singular_locus(lam) > 0) != is_biquanimous(lam)]
    check(u"   ... and V_lambda is singular exactly when lambda is biquanimous", not bad)

    # -- the two tables ------------------------------------------------------------
    TABLE5 = {(3, 1, 1): lambda p: 2 + chi(-15, p),
              (2, 2, 1): lambda p: 2 + chi(5, p),
              (2, 1, 1, 1): lambda p: p + 1 - a_p(p),
              (1, 1, 1, 1, 1): lambda p: p * p + (6 + chi(5, p)) * p + 1}
    TABLE6 = {(1,) * 6: lambda p: p ** 3 + 6 * p * p - 4 * p + 1,
              (2,) + (1,) * 4: lambda p: p * p + 3 * p + 1,
              (3, 1, 1, 1): lambda p: p + 1 - chi(-3, p),
              (2, 2, 1, 1): lambda p: 2 * p,
              (4, 1, 1): lambda p: 2 + chi(-1, p),
              (3, 2, 1): lambda p: 2,
              (2, 2, 2): lambda p: 3}
    for name, table in ((u"Table 1", TABLE5), (u"Table 2", TABLE6)):
        bad = [(lam, p) for lam, f in table.items() for p in primes
               if V_points(lam, p) != f(p)]
        check(u"every row of %s with l >= 3 matches a brute-force count of its stratum" % name,
              not bad, u"%d rows x %d primes" % (len(table), len(primes)))
    bad = [p for p in primes if N_lambda((3, 3), p) != p]
    check(u"   ... including the short stratum (3,3), which the corollary does not cover",
          not bad)
    SHORT = [(5,), (4, 1), (3, 2), (6,), (5, 1), (4, 2)]
    bad = [(lam, p) for lam in SHORT for p in primes if N_lambda(lam, p) != 1]
    check(u"   ... and the rows with l <= 2 and N_lambda = 1 in both tables", not bad,
          u"%d rows x %d primes" % (len(SHORT), len(primes)))

    # -- the classification --------------------------------------------------------
    check(u"Theorem 6: C(n,p) is elementary exactly for n = 3, 4, 6",
          classification(40) == [3, 4, 6], u"searched n <= 40")
    witness = {n: (n - 3, 1, 1, 1) for n in range(5, 61) if n != 6}
    check(u"   ... and the witness (n-3,1,1,1) is never biquanimous, 5 <= n <= 60, n != 6",
          not any(is_biquanimous(w) for w in witness.values()))

    # -- the three propositions ----------------------------------------------------
    bad = [p for p in primes if V_points((2, 2, 1, 1, 1), p) != V22111_closed(p)]
    check(u"Proposition 9: the count of V_(2,2,1,1,1)", not bad)
    bad = [p for p in primes if p > 5
           and V_points((2, 2, 1, 1, 1), p) != V22111_collapsed(p)]
    check(u"Corollary 10: the collapsed form, for p > 5", not bad)
    check(u"   ... and the collapse genuinely fails at p = 5",
          V_points((2, 2, 1, 1, 1), 5) == V22111_closed(5) != V22111_collapsed(5),
          u"census %d, uncollapsed %d, collapsed %d"
          % (V_points((2, 2, 1, 1, 1), 5), V22111_closed(5), V22111_collapsed(5)))
    bad = [(a, p) for a in (1, 3, 5, 6, 7, 8, 9, 10, 11) for p in primes
           if (6 * a * (a * a - 4) * (a * a - 16)) % p != 0
           and V_points((a, 1, 1, 1, 1), p) != Va1111_closed(a, p)]
    check(u"Proposition 11: the family V_(a,1,1,1,1), under its hypothesis on p", not bad)
    check(u"   ... and it genuinely fails at a = p = 5, which is why the hypothesis is there",
          V_points((5, 1, 1, 1, 1), 5) != Va1111_closed(5, 5),
          u"census %d, formula %d" % (V_points((5, 1, 1, 1, 1), 5), Va1111_closed(5, 5)))
    bad = [p for p in primes if not cayley_identity_mod_p(p)]
    check(u"Proposition 16: F(Mu) = 24 (u0u1u2 + ...) over F_p", not bad)
    # the four listed nodes must actually be singular points of V_(2,1^4), not just plausible
    bad = []
    for v in CAYLEY_NODES:
        x5 = -(2 * v[0] + v[1] + v[2] + v[3])
        if 2 * v[0] ** 3 + sum(x ** 3 for x in v[1:]) + x5 ** 3 != 0:
            bad.append((v, "not on the surface"))
        # partial derivatives of F = 2x1^3 + x2^3 + x3^3 + x4^3 - (2x1+x2+x3+x4)^3
        s = 2 * v[0] + v[1] + v[2] + v[3]
        grad = [6 * v[0] ** 2 - 6 * s ** 2] + [3 * v[i] ** 2 - 3 * s ** 2 for i in (1, 2, 3)]
        if any(g != 0 for g in grad):
            bad.append((v, "not singular"))
    check(u"   ... the four nodes really are singular points of V_(2,1^4), and are rational",
          not bad, u"%s" % (bad[:2],))
    check(u"   ... and |det M| = 8, so M is invertible over any field of characteristic != 2",
          abs(_det4([[CAYLEY_NODES[j][i] for j in range(4)]
                     for i in range(4)])) == 8)

    print(u"=" * 74)
    print(u"%d ok, %d FAIL" % (_ok, _bad))
    return 1 if _bad else 0


def _det4(M):
    u"""Determinant of a 4x4 integer matrix, by cofactor expansion.  No numpy rounding."""
    def det(m):
        if len(m) == 1:
            return m[0][0]
        return sum((-1) ** j * m[0][j]
                   * det([row[:j] + row[j + 1:] for row in m[1:]])
                   for j in range(len(m)))
    return det(M)


if __name__ == "__main__":
    sys.exit(self_test())
