# -*- coding: utf-8 -*-
"""In what sense is (1) the smallest alphabet of its kind?

The introduction calls mu_t u {z, 1/z} the smallest alphabet containing one letter of each kind.
That is a slogan unless the alternatives are enumerated, and inversion-closure enumerates them: an
inversion-closed extension of mu_t can only add reciprocal pairs {x, 1/x} and the fixed points of
inversion, which are 1 and -1.  Now 1 lies in mu_t for every t, and -1 lies in mu_t exactly when t
is even.  So the only extension smaller than a pair is

        mu_t u {-1},   t odd,   N = t+1,   excess ONE,

and this script evaluates it.  Its values are 0 and +-1: a Littlewood-type statement with no free
variable and nothing to evaluate as a function of anything.  A free reciprocal pair is therefore the
smallest extension that leaves something to compute, which is the precise content of the slogan.

Authors: Carles Marin, Claude (AI assistant)."""
import sys

from mpmath import mp, mpc, exp, matrix, det

mp.dps = 30
TMAX, NMAX = 8, 15


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


def schur(lam, alpha):
    N = len(alpha)
    lam = list(lam) + [0] * N
    b = [lam[j] + N - 1 - j for j in range(N)]
    num = matrix(N, N)
    den = matrix(N, N)
    for i in range(N):
        for j in range(N):
            num[i, j] = alpha[i] ** b[j]
            den[i, j] = alpha[i] ** (N - 1 - j)
    return det(num) / det(den)


if __name__ == "__main__":
    print("the excess-one alphabet  mu_t u {-1}  for odd t")
    print("")
    total = 0
    offending = 0
    for t in range(3, TMAX, 2):
        N = t + 1
        alpha = [exp(2j * mp.pi * k / t) for k in range(t)] + [mpc(-1)]
        counts = {}
        for n in range(0, NMAX):
            for lam in partitions(n, N):
                v = schur(lam, alpha)
                total += 1
                if abs(v) < mp.mpf(10) ** -12:
                    counts["0"] = counts.get("0", 0) + 1
                elif abs(v - 1) < mp.mpf(10) ** -12:
                    counts["+1"] = counts.get("+1", 0) + 1
                elif abs(v + 1) < mp.mpf(10) ** -12:
                    counts["-1"] = counts.get("-1", 0) + 1
                else:
                    offending += 1
                    counts["other"] = counts.get("other", 0) + 1
        print("  t=%d : %s" % (t, ", ".join("%s x%d" % (k, counts[k]) for k in sorted(counts))))
    print("")
    print("  %d shapes over t = 3,5,7 and |lambda| <= %d" % (total, NMAX - 1))
    print("  values outside {0,+1,-1}: %d   <-- 0 means the excess-one alphabet is frozen"
          % offending)
    sys.exit(1 if offending else 0)
