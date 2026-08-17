# -*- coding: utf-8 -*-
"""The sign of [thm:extra], on both of its families.

The theorem used to read  Phi_t(lambda; z) = +- s_lambda(z, 1/z),  and the sign was a computation
recorded in the remark beside it.  It is now an equality, and the two families get it for different
reasons:

  the t-cores    equality is (4), the criterion of Ayyer-Kumari, which is an equality already;
  the extras     eps_lambda = +1, because the residue word of beta has an EVEN number of
                 inversions.  With t = 2m and lambda_2 = m + j,

                     beta = (6m+j, 3m+j, 2m-1, ..., 1, 0),
                     w    = (j, m+j, 2m-1, ..., 1, 0),
                     inv(w) = C(2m,2) + j + (m+j) = 2m^2 + 2j,

                 which is (11).  The other three factors of the short sign (8) are
                 (-1)^floor(t/2), (-1)^(r_A+r_B) = (-1)^(m+2j) and sgn(oriented coupling) = +1,
                 and the first two cancel.

This script recomputes all of that from beta -- the word, its inversion count, the closed form
2m^2+2j, eps by the LONG formula (7) rather than the short one, and finally the equality itself
against a bialternant evaluation.  Nothing here reads the proof.

CONTROLS.  Three, because "eps = +1 on this family" is the kind of claim a broken sign routine
satisfies for free.

  (i)   eps must NOT be identically +1.  The run reports how many two-row shapes in the same range
        carry eps = -1; if that count were zero the theorem's sign claim would be vacuous and this
        script would be measuring nothing.
  (ii)  the closed form 2m^2 + 2j is compared against a deliberate near-miss, 2m^2 + j, which must
        disagree on every member with j > 0.
  (iii) the equality must FAIL off the classification: the run sweeps the two-row shapes that are
        neither a core nor an extra and confirms that not one of them satisfies it.

A fourth check is not a control but a route check: eps is computed by the LONG formula (7), while
the proof in the paper argues through the SHORT one (8).  The two are compared on the family, so a
divergence between the proof's route and the theorem's statement would show up here.

Authors: Carles Marin, Claude (AI assistant)."""
import sys

from mpmath import mp, mpf, mpc, exp, pi

mp.dps = 40
THETA = mpf("0.41")


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - 1 - j for j in range(N)]


def inv_count(word):
    return sum(1 for i in range(len(word)) for j in range(i + 1, len(word)) if word[i] > word[j])


def blocks(lam, t):
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
        return beta, (j1, j2), (j2, j3)
    big.sort()
    return beta, tuple(big[0]), tuple(big[1])


def eps_long(lam, t):
    """eps_lambda by (7), the long form: the theorem's own statement of the sign."""
    st = blocks(lam, t)
    if st is None:
        return None
    beta, (jA1, jA2), (jB1, jB2) = st
    N = t + 2
    a1, a2, b1, b2 = beta[jA1], beta[jA2], beta[jB1], beta[jB2]
    c = a1 + a2 - b1 - b2
    if c == 0:
        return None
    S = [j for j in range(N) if j not in (jA1, jB1)]
    k11 = (-1) ** ((jA1 + 1) + (jB1 + 1) + inv_count([beta[j] % t for j in S]))
    k11 *= 1 if a1 > b1 else -1
    return (-1) ** (t + (N + 1) * N // 2) * k11 * (1 if c > 0 else -1)


def eps_short(lam, t):
    """eps_lambda by (8), the short form: the route the proof of [thm:extra] actually takes."""
    st = blocks(lam, t)
    if st is None:
        return None
    beta, (jA1, jA2), (jB1, jB2) = st
    a1, a2, b1, b2 = beta[jA1], beta[jA2], beta[jB1], beta[jB2]
    c = a1 + a2 - b1 - b2
    if c == 0:
        return None
    rA, rB = sorted((a1 % t, b1 % t))
    sgn_sigma = (-1) ** inv_count([b % t for b in beta])
    return (-1) ** (t // 2) * sgn_sigma * (-1) ** (rA + rB) * (1 if c > 0 else -1)


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
    N = t + 2
    z = exp(theta)
    w = exp(2j * pi / t)
    xs = [w ** k for k in range(t)] + [z, 1 / z]
    beta = beta_of(lam, N)
    num = [[xs[i] ** beta[j] for j in range(N)] for i in range(N)]
    den = [[xs[i] ** (N - 1 - j) for j in range(N)] for i in range(N)]
    return mydet(num, N) / mydet(den, N)


def chi(k, theta=THETA):
    return (exp((k + 1) * theta) - exp(-(k + 1) * theta)) / (exp(theta) - exp(-theta))


def close(x, y):
    return abs(x - y) < mpf(10) ** -20 * max(1, abs(x))


def is_core(lam, t):
    L = list(lam) + [0, 0]
    A, B = L[0] + 1, L[1]
    return B < t and (A < t or A == B + t)


def run():
    print("=" * 78)
    print("THE EXTRA FAMILY: word, inversion count, eps, and the equality")
    print("=" * 78)
    print("   t    lambda          w = residue word          inv(w)  2m^2+2j  eps  Phi = s(z,1/z)")
    bad = nearmiss = shortdiff = 0
    n = njpos = 0
    for t in range(2, 15, 2):
        m = t // 2
        for j in range(m):
            lam = (m + j + 3 * t // 2 - 1, m + j)
            beta = beta_of(lam, t + 2)
            w = [b % t for b in beta]
            wpred = [j, m + j] + list(range(2 * m - 1, -1, -1))
            iw, pred = inv_count(w), 2 * m * m + 2 * j
            e = eps_long(lam, t)
            eq = close(phi_bialt(lam, t), chi(lam[0] - lam[1]))
            n += 1
            if not (w == wpred and iw == pred and e == 1 and eq):
                bad += 1
            if eps_short(lam, t) != e:                       # route check, (8) against (7)
                shortdiff += 1
            if j > 0:                                        # control (ii): the near-miss formula
                njpos += 1
                nearmiss += (iw != 2 * m * m + j)
            if t <= 8:
                print("  %3d   %-14s %-24s %6d  %7d  %+3d  %s"
                      % (t, str(lam), str(w) if t <= 6 else "(as predicted)", iw, pred, e, eq))
    print("  ... t = 2,4,...,14: %d members, %d failures" % (n, bad))
    print("  ROUTE:      the short form (8) disagrees with the long form (7) on %d of %d"
          % (shortdiff, n))
    print("  CONTROL ii: the near-miss 2m^2+j is refuted on %d of the %d members with j>0"
          % (nearmiss, njpos))

    print()
    print("=" * 78)
    print("THE t-CORE FAMILY, and the shapes that are neither")
    print("=" * 78)
    ncore = badcore = nother = badother = nneg = ntwo = 0
    for t in range(2, 11):
        for l2 in range(0, 2 * t + 2):
            for l1 in range(l2, 4 * t + 4):
                lam = (l1, l2) if l2 else ((l1,) if l1 else ())
                m = t // 2
                extra = (t % 2 == 0 and l2 and l1 == l2 + 3 * t // 2 - 1 and m <= l2 <= t - 1)
                eq = close(phi_bialt(lam, t), chi(l1 - l2))
                e = eps_long(lam, t)
                if e is not None:
                    ntwo += 1
                    nneg += (e == -1)
                if is_core(lam, t):
                    ncore += 1
                    badcore += (not eq)
                elif not extra:
                    nother += 1
                    badother += eq
    print("  t <= 10:  %d cores, %d failures of the EQUALITY" % (ncore, badcore))
    print("  CONTROL i:   eps = -1 on %d of the %d two-row shapes in range" % (nneg, ntwo))
    print("               <-- a nonzero count; eps is not identically +1, so 'eps = +1 on the")
    print("                   extra family' is a statement and not a property of the routine")
    print("  CONTROL iii: %d two-row shapes are neither core nor extra, and %d of them satisfy"
          % (nother, badother))
    print("               the equality anyway  <-- 0 means the classification is not vacuous")
    dead = (1 if nearmiss == 0 else 0) + (1 if nneg == 0 else 0)
    return bad + badcore + badother + shortdiff + dead


if __name__ == "__main__":
    rc = run()
    print()
    print("TOTAL FAILURES: %d" % rc)
    sys.exit(1 if rc else 0)
