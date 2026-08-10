# -*- coding: utf-8 -*-
"""The second half of the involution question, made explicit (t = 2).

Splitting the alphabet as in (13),

    s_lambda(1,-1,z,1/z) = sum_nu s_{lambda/nu}(1,-1) * chi_{nu1-nu2}(z),   ell(nu) <= 2,

an involution carrying the free parameter must preserve k = nu1 - nu2, because that is the index of
the chi its term sits on.  For fixed k the admissible nu are exactly (m+k, m) for m >= 0 -- a single
one-parameter family -- so the cancellation across nu is a statement about the word

    sigma_m(k) = s_{lambda/(m+k,m)}(1,-1),      m = 0, 1, 2, ...

computed here by Jacobi-Trudi with h_j(1,-1) = [j >= 0 and j even], which is exact.

Three measurements:

  A  max |sigma|.  If it is 1, the terms are a SET and the required cancellation is a matching
     rather than something carrying multiplicity.

  B  the sign alternates inside every maximal run of consecutive nonzero sigma.  This is what makes
     "pair m with m+1 inside a run" sign-reversing, and it is the candidate involution.

  C  the survivors -- one per odd-length run, carrying that run's first sign -- rebuilt into the chi
     basis and compared with the paper's closed form eps * chi_{d1/2-1} chi_{d2/2-1} chi_{d3/2-1},
     coefficient by coefficient.  This is the control: it is what would fail if the matching threw
     away the wrong terms.

Authors: Carles Marin, Claude (AI assistant)."""
from collections import defaultdict

from theorem_full import setup

MAXSIZE, NROWS = 17, 4


def partitions(n, maxlen):
    if n == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for first in range(n, 0, -1):
        for rest in partitions(n - first, maxlen - 1):
            if not rest or rest[0] <= first:
                yield (first,) + rest


def det(M):
    n = len(M)
    if n <= 1:
        return M[0][0] if n else 1
    tot = 0
    for j in range(n):
        if M[0][j]:
            tot += (-1) ** j * M[0][j] * det([r[:j] + r[j + 1:] for r in M[1:]])
    return tot


def skew_at_pm1(lam, nu, n=NROWS):
    """s_{lam/nu}(1,-1), by Jacobi-Trudi with h_j(1,-1) = [j >= 0 and j even]."""
    lam = list(lam) + [0] * n
    nu = list(nu) + [0] * n
    return det([[1 if (lam[i] - nu[j] - i + j) >= 0 and (lam[i] - nu[j] - i + j) % 2 == 0 else 0
                 for j in range(n)] for i in range(n)])


def sequences(lam):
    """{k: [sigma_0, sigma_1, ...]} for every k carrying a nonzero term."""
    out = {}
    for k in range(sum(lam) + 1):
        seq, m = [], 0
        while m + k <= (list(lam) + [0])[0] and m <= (list(lam) + [0, 0])[1]:
            seq.append(skew_at_pm1(lam, (m + k, m)))
            m += 1
        if any(seq):
            out[k] = seq
    return out


def chi_mult(poly, a):
    out = defaultdict(int)
    for k, c in poly.items():
        for j in range(abs(k - a), k + a + 1, 2):
            out[j] += c
    return {k: v for k, v in out.items() if v}


def closed_form_chi(lam):
    st = setup(lam, 2)
    if st is None:
        return {}
    beta, Ac, Bc = st
    a1, a2, b1, b2 = beta[Ac[0]], beta[Ac[1]], beta[Bc[0]], beta[Bc[1]]
    d1, d2, d3 = a1 - a2, b1 - b2, abs(a1 + a2 - b1 - b2)
    if d3 == 0:
        return {}
    return chi_mult(chi_mult({d1 // 2 - 1: 1}, d2 // 2 - 1), d3 // 2 - 1)


def runs_of(seq):
    runs, cur = [], []
    for s in seq:
        if s:
            cur.append(s)
        else:
            if cur:
                runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    return runs


if __name__ == "__main__":
    maxabs = nruns = same_adj = odd = even = ok = bad = 0
    shapes = 0
    for size in range(1, MAXSIZE):
        for lam in partitions(size, NROWS):
            seqs = sequences(lam)
            if not seqs:
                continue
            shapes += 1
            survived = {}
            for k, seq in seqs.items():
                maxabs = max(maxabs, max(abs(s) for s in seq))
                tot = 0
                for r in runs_of(seq):
                    nruns += 1
                    if any(r[i] == r[i + 1] for i in range(len(r) - 1)):
                        same_adj += 1
                    if len(r) % 2:
                        odd += 1
                        tot += r[0]
                    else:
                        even += 1
                if tot:
                    survived[k] = tot
            want = closed_form_chi(lam)
            same = ({k: abs(v) for k, v in survived.items()}
                    == {k: abs(v) for k, v in want.items() if v})
            ok, bad = (ok + 1, bad) if same else (ok, bad + 1)

    print("range: |lambda| <= %d, at most %d rows, %d shapes with a nonzero term"
          % (MAXSIZE - 1, NROWS, shapes))
    print()
    print("A  max |sigma|:                                   %d   <-- 1 means a set, not a multiset"
          % maxabs)
    print("B  runs of consecutive nonzeros:                  %d" % nruns)
    print("   runs with two EQUAL adjacent signs:            %d   <-- 0 means the sign alternates"
          % same_adj)
    print("   odd-length runs %d, even-length runs %d" % (odd, even))
    print()
    print("C  survivors vs the closed form (up to the global sign): %d agree, %d disagree"
          % (ok, bad))
