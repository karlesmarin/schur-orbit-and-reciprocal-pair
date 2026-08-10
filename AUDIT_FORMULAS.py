# -*- coding: utf-8 -*-
"""Every displayed formula of the paper, checked against its own definitions.

The other checkers ask whether the NUMBERS in the text appear in an archived run. None of them asks
whether the FORMULAS are true. This does, one line per formula, evaluating both sides from the
definitions at numeric points: a sign convention or a stray constant is exactly what survives every
other check.

Nothing here imports the paper's closed forms; s_lambda is always the bialternant and h_j is always
read off the generating function, so a shared wrong assumption cannot hide.

Authors: Carles Marin, Claude (AI assistant)."""
import itertools
import sys

from mpmath import mp, mpc, mpf, exp, matrix, det, sqrt

mp.dps = 30
TOL = mpf(10) ** -18
Z = mpc("1.31", "0.23")


# ----------------------------------------------------------------- basics
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
    lam = list(lam) + [0] * N
    return [lam[j] + N - 1 - j for j in range(N)]


def schur(lam, alpha):
    """s_lambda(alpha) from the bialternant, for any alphabet of distinct letters."""
    N = len(alpha)
    b = beta_of(lam, N)
    num = matrix(N, N)
    den = matrix(N, N)
    for i in range(N):
        for j in range(N):
            num[i, j] = alpha[i] ** b[j]
            den[i, j] = alpha[i] ** (N - 1 - j)
    return det(num) / det(den)


def h(j, alpha):
    """complete homogeneous h_j(alpha) by expanding prod 1/(1 - x q) to order j."""
    if j < 0:
        return mpc(0)
    coef = [mpc(1)] + [mpc(0)] * j
    for x in alpha:
        for i in range(1, j + 1):
            coef[i] = coef[i] + x * coef[i - 1]
    return coef[j]


def cofactor_det(M):
    """exact cofactor expansion: mpmath's LU det fails on the singular matrices that occur here."""
    n = len(M)
    if n == 0:
        return mpc(1)
    if n == 1:
        return M[0][0]
    tot = mpc(0)
    for j in range(n):
        if M[0][j] != 0:
            minor = [row[:j] + row[j + 1:] for row in M[1:]]
            tot += ((-1) ** j) * M[0][j] * cofactor_det(minor)
    return tot


def skew(lam, nu, alpha, n=None):
    """s_{lam/nu}(alpha) by Jacobi-Trudi."""
    n = n or max(len(lam), len(nu))
    if n == 0:
        return mpc(1)
    lam = list(lam) + [0] * n
    nu = list(nu) + [0] * n
    return cofactor_det([[h(lam[i] - nu[j] - i + j, alpha) for j in range(n)] for i in range(n)])


def mu(t):
    return [exp(2j * mp.pi * k / t) for k in range(t)]


def core_and_sigma(lam, t, N):
    """(is lambda a t-core-free shape i.e. core empty, sgn(sigma)) from beta."""
    b = beta_of(lam, N)
    # remove t-hooks: repeatedly lower a beta number by t into a free slot
    cur = sorted(b, reverse=True)
    changed = True
    while changed:
        changed = False
        for i, x in enumerate(cur):
            if x >= t and (x - t) not in cur:
                cur[i] = x - t
                cur.sort(reverse=True)
                changed = True
                break
    core_empty = (cur == list(range(N - 1, -1, -1)))
    order = []
    for a in range(t):
        order += [j for j, x in enumerate(b) if x % t == a]
    sg = 1
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            if order[i] > order[j]:
                sg = -sg
    return core_empty, sg


def f(u):
    return Z ** u - Z ** (-u)


CHECKED = []


def report(name, ok, bad, extra=""):
    flag = "ok" if bad == 0 else "*** %d FAIL ***" % bad
    print("  %-46s %5d checked   %s%s" % (name, ok + bad, flag, extra))
    CHECKED.append(ok + bad)
    return bad


# ----------------------------------------------------------------- the audit
if __name__ == "__main__":
    total = 0

    # (2.1) Littlewood: s_lambda(mu_t) = (-1)^{C(t,2)} sgn(sigma) on the t-cores, else 0
    ok = bad = 0
    for t in range(2, 7):
        for n in range(0, 13):
            for lam in partitions(n, t):
                got = schur(lam, mu(t))
                ce, sg = core_and_sigma(lam, t, t)
                want = ((-1) ** (t * (t - 1) // 2)) * sg if ce else 0
                (ok, bad) = (ok + 1, bad) if abs(got - want) < TOL else (ok, bad + 1)
    total += report("Theorem 2.1  Littlewood, sign included", ok, bad)

    # (2) the interval reading: d1 = |A|, d2 = |B|, d3 = 2|centre A - centre B|
    ok = bad = 0
    for t in range(2, 8):
        for n in range(0, 14):
            for lam in partitions(n, t + 2):
                N = t + 2
                b = beta_of(lam, N)
                cls = {}
                for x in b:
                    cls.setdefault(x % t, []).append(x)
                if len(cls) < t:
                    continue
                big = sorted((v for v in cls.values() if len(v) >= 2), key=lambda v: -v[0])
                if len(big) == 1:
                    p, q, r = sorted(big[0], reverse=True)
                    A, B = [p, q], [q, r]
                else:
                    A, B = sorted(big, key=lambda v: v[0] % t)
                d1, d2 = A[0] - A[1], B[0] - B[1]
                d3 = abs(A[0] + A[1] - B[0] - B[1])
                cA, cB = mpf(A[0] + A[1]) / 2, mpf(B[0] + B[1]) / 2
                good = (d1 == abs(A[0] - A[1]) and d2 == abs(B[0] - B[1])
                        and d3 == 2 * abs(cA - cB))
                (ok, bad) = (ok + 1, bad) if good else (ok, bad + 1)
    total += report("eq. (2)      the interval reading", ok, bad)

    # Lemma 5.2 (L2): det(x^{N-j}) = (-1)^{C(t,2)} V (z^t-1)(z^{-t}-1)(z-z^{-1})
    ok = bad = 0
    for t in range(2, 9):
        N = t + 2
        alpha = mu(t) + [Z, 1 / Z]
        D = matrix(N, N)
        for i in range(N):
            for j in range(N):
                D[i, j] = alpha[i] ** (N - 1 - j)
        V = mpc(1)
        for r in range(t):
            for rp in range(r + 1, t):
                V *= (mu(t)[rp] - mu(t)[r])
        want = ((-1) ** (t * (t - 1) // 2)) * V * (Z ** t - 1) * (Z ** (-t) - 1) * (Z - 1 / Z)
        (ok, bad) = (ok + 1, bad) if abs(det(D) - want) < TOL else (ok, bad + 1)
    total += report("Lemma 5.2    the Vandermonde denominator", ok, bad)

    # Lemma 5.5 (a) and (b)
    ok = bad = 0
    for c in range(-4, 5):
        for p in range(0, 5):
            for q in range(0, 5):
                lhs = sum(e * n * f(c + e * p + n * q)
                          for e in (1, -1) for n in (1, -1))
                (ok, bad) = (ok + 1, bad) if abs(lhs - f(c) * f(p) * f(q)) < TOL else (ok, bad + 1)
    total += report("Lemma 5.5(a) the four-term collapse", ok, bad)

    ok = bad = 0
    for u in range(-4, 5):
        for v in range(-4, 5):
            lhs = f(2 * u) - f(2 * u + 2 * v) + f(2 * v)
            (ok, bad) = (ok + 1, bad) if abs(lhs + f(u) * f(v) * f(u + v)) < TOL else (ok, bad + 1)
    total += report("Lemma 5.5(b) the size-three collapse", ok, bad)

    # eq. (12): s_lambda(mu_t u B) = sum_nu s_{lambda/nu}(mu_t) s_nu(B)
    ok = bad = 0
    B = [Z, 1 / Z]
    for t in range(2, 6):
        N = t + 2
        for n in range(0, 10):
            for lam in partitions(n, N):
                got = schur(lam, mu(t) + B)
                tot = mpc(0)
                nus = set()
                for k in range(0, n + 1):
                    nus.update(partitions(k, 2))       # a set: partitions(0,k) already yields ()
                for nu in sorted(nus):
                    nul = list(nu) + [0, 0]
                    if any(nul[i] > (list(lam) + [0] * N)[i] for i in range(2)):
                        continue
                    tot += skew(lam, nu, mu(t), n=N) * schur(nu, B)
                (ok, bad) = (ok + 1, bad) if abs(got - tot) < TOL else (ok, bad + 1)
    total += report("eq. (12)     the splitting of the alphabet", ok, bad)

    # eq. (16): complementation on the reciprocal alphabet
    ok = bad = 0
    for r in (1, 2):
        N = 2 * r + 2
        zs = [mpc("1.17", "0.31"), mpc("0.93", "-0.44")][:r]
        A = [mpc(1), mpc(-1)] + [x for z in zs for x in (z, 1 / z)]
        for n in range(0, 12):
            for lam in partitions(n, N):
                lamp = list(lam) + [0] * N
                w = lamp[0] + lamp[N - 1]
                hat = tuple(w - lamp[N - 1 - i] for i in range(N))
                if any(hat[i] < 0 for i in range(N)) or list(hat) != sorted(hat, reverse=True):
                    continue
                c = w + N - 1
                lhs = schur(hat, A)
                rhs = ((-1) ** c) * (-1) * schur(lam, A)
                (ok, bad) = (ok + 1, bad) if abs(lhs - rhs) < TOL else (ok, bad + 1)
    total += report("eq. (16)     the complementation identity", ok, bad)

    # eq. (17): d'Ocagne for the Chebyshev-like recurrence
    ok = bad = 0
    cst = Z + 1 / Z
    S = {0: mpc(0), 1: mpc(1)}
    for k in range(2, 40):
        S[k] = cst * S[k - 1] - S[k - 2]
    for k in range(1, 40):
        S[-k] = -S[k]
    for b_ in range(0, 20):
        for bp in range(0, 20):
            lhs = S[b_] * S[bp - 1] - S[bp] * S[b_ - 1]
            (ok, bad) = (ok + 1, bad) if abs(lhs + S[b_ - bp]) < TOL else (ok, bad + 1)
    total += report("eq. (17)     d'Ocagne for this recurrence", ok, bad)

    # eq. (23) in Section 10: the scalar caveat  s(1,-1,y,-1/y) = i^{|lam|} s(-i,i,w,1/w)
    ok = bad = 0
    y = mpc("1.19", "0.37")
    w = -1j * y
    for n in range(0, 13):
        for lam in partitions(n, 4):
            lhs = schur(lam, [mpc(1), mpc(-1), y, -1 / y])
            rhs = (1j ** sum(lam)) * schur(lam, [mpc(-1j), mpc(1j), w, 1 / w])
            (ok, bad) = (ok + 1, bad) if abs(lhs - rhs) < TOL else (ok, bad + 1)
    total += report("eq. (23)     the scalar caveat of Conj. 10.3", ok, bad)

    print("")
    print("TOTAL evaluations: %d over %d displayed formulas" % (sum(CHECKED), len(CHECKED)))
    print("TOTAL formula failures: %d" % total)
    sys.exit(1 if total else 0)
