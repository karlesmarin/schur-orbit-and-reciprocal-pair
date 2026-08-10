# -*- coding: utf-8 -*-
"""Is the fixed set of Kumari's involution a single tableau whenever sigma != 0?

Her Lemma 2.18 defines a fixed-point-free involution gamma on the NON-coverable tableaux of
SSYT_{2n}(lambda/mu), so the signed sum reduces to the COVERABLE ones; Lemma 2.16 puts those in
bijection with domino tableaux, and Lemma 2.14 says the parity of the number of vertical dominoes
is independent of the tableau -- so all survivors carry the same sign.

At n = 1 the alphabet is (x, -x) = (1, -1) and the weight is (-1)^{c_2(T)}, which is exactly
sigma = s_{lambda/nu}(1,-1).  The paper's Section 6 asserts that the surviving set is then a single
tableau when sigma != 0 and empty otherwise.  That is checked here by brute force: enumerate every
2-letter semistandard filling of the skew shape, keep the coverable ones, and compare the count with
sigma computed independently by Jacobi-Trudi.

A tableau with entries in {1,2} is coverable when it can be tiled by dominoes of the shapes
    [2|2]  (two 2's)   and   [1 over 2]  (a 1 above a 2),
which are Definition 2.15 at a = 1.
"""
import itertools
import sys


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


def cells(lam, nu):
    nu = list(nu) + [0] * len(lam)
    return [(i, j) for i in range(len(lam)) for j in range(nu[i], lam[i])]


def fillings(lam, nu):
    """every semistandard filling with entries in {1,2}: rows weakly increase, columns strictly."""
    cs = cells(lam, nu)
    idx = {c: k for k, c in enumerate(cs)}
    for vals in itertools.product((1, 2), repeat=len(cs)):
        f = {c: vals[idx[c]] for c in cs}
        ok = True
        for (i, j) in cs:
            if (i, j + 1) in f and f[(i, j)] > f[(i, j + 1)]:
                ok = False
                break
            if (i + 1, j) in f and f[(i, j)] >= f[(i + 1, j)]:
                ok = False
                break
        if ok:
            yield f


def coverable(f):
    """can the filling be tiled by [2|2] horizontal and [1 over 2] vertical dominoes?"""
    cs = sorted(f)
    n = len(cs)
    if n % 2:
        return False
    pos = set(cs)

    def rec(remaining):
        if not remaining:
            return True
        c = min(remaining)
        i, j = c
        # vertical domino: 1 on top of 2
        if f[c] == 1 and (i + 1, j) in remaining and f[(i + 1, j)] == 2:
            if rec(remaining - {c, (i + 1, j)}):
                return True
        # horizontal domino: two 2's side by side
        if f[c] == 2 and (i, j + 1) in remaining and f[(i, j + 1)] == 2:
            if rec(remaining - {c, (i, j + 1)}):
                return True
        return False

    return rec(frozenset(pos))


def det(M):
    n = len(M)
    if n <= 1:
        return M[0][0] if n else 1
    tot = 0
    for j in range(n):
        if M[0][j]:
            tot += (-1) ** j * M[0][j] * det([r[:j] + r[j + 1:] for r in M[1:]])
    return tot


def sigma(lam, nu, n):
    """s_{lam/nu}(1,-1) by Jacobi-Trudi with h_j(1,-1) = [j >= 0 and j even]."""
    lam = list(lam) + [0] * n
    nu = list(nu) + [0] * n
    return det([[1 if (lam[i] - nu[j] - i + j) >= 0 and (lam[i] - nu[j] - i + j) % 2 == 0 else 0
                 for j in range(n)] for i in range(n)])


if __name__ == "__main__":
    NROWS, MAXSIZE = 4, 11
    checked = mismatch = 0
    nonzero_many = 0
    worst = []
    for size in range(1, MAXSIZE):
        for lam in partitions(size, NROWS):
            for k in range(0, size + 1):
                for nu in list(partitions(k, 2)) if k else [()]:
                    nul = list(nu) + [0] * NROWS
                    laml = list(lam) + [0] * NROWS
                    if any(nul[i] > laml[i] for i in range(NROWS)):
                        continue
                    s = sigma(lam, nu, NROWS)
                    cov = sum(1 for f in fillings(lam, nu) if coverable(f))
                    checked += 1
                    if abs(s) != cov:
                        mismatch += 1
                        if len(worst) < 6:
                            worst.append((lam, nu, s, cov))
                    if s != 0 and cov != 1:
                        nonzero_many += 1
    print("skew shapes checked (at most %d rows, |lambda| <= %d): %d" % (NROWS, MAXSIZE - 1, checked))
    print("")
    print("  |sigma| != #coverable :          %d   <-- 0 means the survivors carry one sign" % mismatch)
    print("  sigma != 0 but #coverable != 1 : %d   <-- 0 CONFIRMS the paper's claim" % nonzero_many)
    for w in worst:
        print("     lam=%s nu=%s sigma=%s coverable=%s" % w)
    sys.exit(1 if (mismatch or nonzero_many) else 0)
