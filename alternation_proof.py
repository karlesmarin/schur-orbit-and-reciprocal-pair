# -*- coding: utf-8 -*-
"""The alternation of Question 6.4, proved -- each step of the proof checked separately.

Write beta_i = lam_i - i and c_j = nu_j - j.  Jacobi-Trudi with h_j(1,-1) = [j >= 0 and j even] is

    M_ij = [beta_i >= c_j] . [beta_i = c_j mod 2],        sigma = det M .

STEP 1  that rewriting reproduces s_{lam/nu}(1,-1).

STEP 2  every column is supported on ONE parity class of rows, and there it is the indicator of an
        initial segment (the rows of a class are ordered by decreasing beta).  So M is block
        diagonal, and det M != 0 forces the number of columns of each parity to equal the size of
        that row class, and the segment lengths inside a class to be a permutation of 1..size.
        Hence |sigma| <= 1, which was measurement A.

STEP 3  with nu = (m+k, m) one has c_1 - c_2 = k+1.  For k ODD the two moving columns share a
        parity, so passing m -> m+1 moves both to the other class and changes that class's column
        count by 2; the count cannot match the row class twice.  So for k odd no two consecutive
        sigma are nonzero -- every run has length one and the alternation is vacuous there.

STEP 4  for k EVEN the two columns have opposite parities and swap classes.  If sigma_m and
        sigma_{m+1} are both nonzero then the segment lengths are exchanged, l'_1 = l_2 and
        l'_2 = l_1, so each class sees the SAME length sequence in the same column order and the two
        block determinants are unchanged.  The only thing that changes is the shuffle separating the
        two classes of columns, and columns 1 and 2 are adjacent, so its sign flips exactly once.

STEP 5  therefore sigma_{m+1} = -sigma_m whenever both are nonzero, which is the alternation.

Authors: Carles Marin, Claude (AI assistant)."""
from involution_runs import partitions, det, skew_at_pm1

MAXSIZE, NROWS = 15, 4


def cols(lam, nu, n=NROWS):
    lam = list(lam) + [0] * n
    nu = list(nu) + [0] * n
    beta = [lam[i] - i for i in range(n)]
    c = [nu[j] - j for j in range(n)]
    return beta, c


def M_beta(lam, nu, n=NROWS):
    beta, c = cols(lam, nu, n)
    return [[1 if (beta[i] >= c[j] and (beta[i] - c[j]) % 2 == 0) else 0
             for j in range(n)] for i in range(n)]


def classes(lam, nu, n=NROWS):
    """(row class sizes by parity, per-column (parity, segment length))."""
    beta, c = cols(lam, nu, n)
    size = {0: sum(1 for b in beta if b % 2 == 0), 1: sum(1 for b in beta if b % 2 == 1)}
    info = []
    for j in range(n):
        e = c[j] % 2
        info.append((e, sum(1 for b in beta if b % 2 == e and b >= c[j])))
    return size, info


def perfect(lam, nu, n=NROWS):
    """does every parity class carry exactly its own size in columns, with lengths 1..size?"""
    size, info = classes(lam, nu, n)
    for e in (0, 1):
        lens = sorted(l for (p, l) in info if p == e)
        if lens != list(range(1, size[e] + 1)):
            return False
    return True


if __name__ == "__main__":
    s1 = s2 = s3 = s4 = 0
    b1 = b2 = b3 = b4 = 0
    maxabs = 0
    for size in range(1, MAXSIZE):
        for lam in partitions(size, NROWS):
            for k in range(sum(lam) + 1):
                seq = []
                m = 0
                while m + k <= (list(lam) + [0])[0] and m <= (list(lam) + [0, 0])[1]:
                    nu = (m + k, m)
                    # STEP 1
                    a, b = det(M_beta(lam, nu)), skew_at_pm1(lam, nu)
                    s1, b1 = (s1 + 1, b1) if a == b else (s1, b1 + 1)
                    # STEP 2
                    maxabs = max(maxabs, abs(a))
                    if (a != 0) == perfect(lam, nu):
                        s2 += 1
                    else:
                        b2 += 1
                    seq.append(a)
                    m += 1
                # STEP 3
                if k % 2 == 1:
                    if any(seq[i] and seq[i + 1] for i in range(len(seq) - 1)):
                        b3 += 1
                    else:
                        s3 += 1
                # STEP 4 and 5
                for i in range(len(seq) - 1):
                    if seq[i] and seq[i + 1]:
                        _, i0 = classes(lam, (i + k, i))
                        _, i1 = classes(lam, (i + 1 + k, i + 1))
                        swapped = (i1[0][1] == i0[1][1] and i1[1][1] == i0[0][1])
                        s4, b4 = (s4 + 1, b4) if swapped else (s4, b4 + 1)
                        if seq[i + 1] != -seq[i]:
                            b3 += 1

    print("range: |lambda| <= %d, at most %d rows" % (MAXSIZE - 1, NROWS))
    print()
    print("STEP 1  beta-form of Jacobi-Trudi = s_{lam/nu}(1,-1):     %d ok, %d fail" % (s1, b1))
    print("STEP 2  sigma != 0  <=>  every class perfect:             %d ok, %d fail" % (s2, b2))
    print("        max |sigma| (STEP 2 predicts 1):                  %d" % maxabs)
    print("STEP 3  k odd => no two consecutive nonzeros:             %d ok, %d fail" % (s3, b3))
    print("STEP 4  consecutive nonzeros => lengths exchanged:        %d ok, %d fail" % (s4, b4))
    print()
    print("STEP 5  sigma_{m+1} = -sigma_m on every consecutive pair: %d pairs, %d fail"
          % (s4, b4))
