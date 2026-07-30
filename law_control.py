"""
Independent control of the general-m law for  s_lambda(mu_m, t, 1/t).

Authors: Carles Marin, Claude (AI assistant).

Written from the STATEMENT only (TWINING.md 2d/2e), not from any existing script:
a fresh bialternant evaluation in plain python + mpmath, compared against

    Phi_m(lambda;t) = eps * prod_i sinh(d_i*theta/2) / [ sinh(m*theta/2)^2 * sinh(theta) ]

with d read off the residue classes mod m of the beta-set.
Purpose: confirm magnitude + that the ratio is exactly +-1, out of our own code path.
"""
import itertools
from mpmath import mp, mpf, mpc, exp, sinh, e, pi, matrix, det, sqrt

mp.dps = 40


def partitions(n, maxlen):
    """all partitions of n with at most maxlen parts, decreasing"""
    def rec(rem, cap, k):
        if rem == 0:
            yield []
            return
        if k == 0:
            return
        for p in range(min(rem, cap), 0, -1):
            for tail in rec(rem - p, p, k - 1):
                yield [p] + tail
    return list(rec(n, n, maxlen))


def mydet(M, n):
    """Gaussian elimination with partial pivoting; returns 0 on an exactly singular matrix
    (mpmath's own det raises there, and the singular case is precisely the one we must score)."""
    A = [row[:] for row in M]
    sign = 1
    d = mpc(1)
    for c in range(n):
        piv = max(range(c, n), key=lambda i: abs(A[i][c]))
        if abs(A[piv][c]) < mpf("1e-30"):
            return mpc(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            sign = -sign
        d *= A[c][c]
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            for j in range(c, n):
                A[i][j] -= f * A[c][j]
    return sign * d


def phi_numeric(lam, m, theta):
    """s_lambda evaluated at (1, w, ..., w^{m-1}, t, 1/t) by the bialternant."""
    N = m + 2
    t = exp(theta)
    w = exp(2j * pi / m)
    xs = [w ** k for k in range(m)] + [t, 1 / t]
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[j] + N - 1 - j for j in range(N)]          # decreasing
    num = [[xs[i] ** beta[j] for j in range(N)] for i in range(N)]
    den = [[xs[i] ** (N - 1 - j) for j in range(N)] for i in range(N)]
    return mydet(num, N) / mydet(den, N)


def law_prediction(lam, m, theta):
    """(d-multiset, |value|) from the residue profile; None if a class is empty (=> 0)."""
    N = m + 2
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[j] + N - 1 - j for j in range(N)]
    classes = {}
    for b in beta:
        classes.setdefault(b % m, []).append(b)
    if len(classes) < m:                       # a residue class is empty
        return None
    sizes = sorted(len(v) for v in classes.values())
    big = [sorted(v, reverse=True) for v in classes.values() if len(v) >= 2]
    if sizes[-1] == 3:                          # one class of size 3
        p, q, r = big[0]
        A, B = (p, q), (q, r)
    else:                                       # two classes of size 2
        A, B = tuple(big[0]), tuple(big[1])
    a1, a2 = A
    b1, b2 = B
    d = (a1 - a2, b1 - b2, abs(a1 + a2 - b1 - b2))
    val = (sinh(d[0] * theta / 2) * sinh(d[1] * theta / 2) * sinh(d[2] * theta / 2)
           / (sinh(m * theta / 2) ** 2 * sinh(theta)))
    return d, val


def run():
    theta = mpf("0.41")                          # generic
    tol = mpf("1e-20")
    tot_nz = tot_z = fail = badratio = 0
    for m in range(2, 8):
        N = m + 2
        maxsize = {2: 14, 3: 14, 4: 13, 5: 12, 6: 11, 7: 10}[m]
        for n in range(0, maxsize + 1):
            for lam in partitions(n, N):
                phi = phi_numeric(lam, m, theta)
                pred = law_prediction(lam, m, theta)
                if pred is None:
                    if abs(phi) > 1e-12:
                        fail += 1
                        print("  ZERO-PREDICTION FAILS", m, lam, phi)
                    else:
                        tot_z += 1
                    continue
                d, val = pred
                if abs(val) < 1e-18:             # law says 0 via d_3 = 0
                    if abs(phi) > 1e-12:
                        fail += 1
                        print("  CONCENTRIC-ZERO FAILS", m, lam, phi)
                    else:
                        tot_z += 1
                    continue
                ratio = phi / val
                if abs(ratio.imag) > 1e-12 or abs(abs(ratio.real) - 1) > 1e-12:
                    badratio += 1
                    if badratio < 6:
                        print("  RATIO NOT +-1", m, lam, "d=", d, "ratio=", ratio)
                else:
                    tot_nz += 1
        print(f"m={m}: running totals  nonzero-match={tot_nz}  zeros={tot_z} "
              f"  fail={fail}  bad-ratio={badratio}")
    print()
    print(f"TOTAL: {tot_nz} nonzero matched (ratio exactly +-1), {tot_z} zeros predicted, "
          f"{fail} vanishing failures, {badratio} magnitude failures")


if __name__ == "__main__":
    run()
