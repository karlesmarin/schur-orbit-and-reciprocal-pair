"""extra_locus_kernel.py -- Problem 10.6 asks whether the extra locus can be read off the KERNEL of
the specialization.  This answers part of it: not from the kernel's existence, but from its size.

Three free parts of size two are compared, each substituted into the criterion
s_lambda(mu_t, F) = +- s_lambda(F) for l(lambda) <= 2:

    F = (z, 1/z)   inversion-closed.  s_mu(F) = chi_{mu_1-mu_2}(z) keeps one grading.
    F = (z, -z)    a zeta_2-orbit.    s_mu(F) = z^{|mu|} s_mu(1,-1) keeps only |mu| and a sign.
    F = (z, w)     free.              no kernel at all.

The first is Theorem 5.2.  The third is Remark 5.3's own control and must return the t-cores and
nothing else -- if it does not, this file is measuring the wrong thing.

Everything is exact over Q: an alphabet enters only through k -> p_k, and s_lambda is the
Jacobi-Trudi determinant with h from the power sums by Newton.  No closed form of the paper is used.

Authors: Carles Marin, Claude (AI assistant).
"""
from fractions import Fraction as F


def two_row_shapes(maxpart):
    out = [()]
    for a in range(1, maxpart + 1):
        out.append((a,))
        for b in range(1, a + 1):
            out.append((a, b))
    return out


def det_frac(M):
    n = len(M)
    M = [row[:] for row in M]
    d = F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d *= M[c][c]
        pv = M[c][c]
        for r in range(c + 1, n):
            if M[r][c] == 0:
                continue
            f = M[r][c] / pv
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d


def schur(lam, P):
    L = [k for k in lam if k > 0]
    n = len(L)
    if n == 0:
        return F(1)
    M = max(L) + n
    h = [F(1)]
    for m in range(1, M + 1):
        h.append(sum(P(i) * h[m - i] for i in range(1, m + 1)) / m)
    H = lambda d: F(0) if d < 0 else h[d]
    return det_frac([[H(L[i] - (i + 1) + (j + 1)) for j in range(n)] for i in range(n)])


def is_core(lam, t):
    N = t + 2
    L = list(lam) + [0] * (N - len(lam))
    b = set(L[j] + N - (j + 1) for j in range(N))
    return not any(v >= t and (v - t) not in b for v in b)


def is_extra(lam, t):
    """the family of [thm:extra]"""
    L = list(lam) + [0, 0]
    return (t % 2 == 0 and len([k for k in lam if k > 0]) == 2
            and L[0] == L[1] + 3 * t // 2 - 1 and t // 2 <= L[1] <= t - 1)


PTS = [(F(3, 2), F(5, 3)), (F(7, 4), F(2)), (F(5, 2), F(3))]

FREE = [
    ("(z,1/z)", lambda z, w: (lambda k: z ** k + F(1) / z ** k)),
    ("(z,-z)", lambda z, w: (lambda k: z ** k + (-z) ** k)),
    ("(z,w)", lambda z, w: (lambda k: z ** k + w ** k)),
]


def main():
    bar = "=" * 76
    print(bar)
    print("Problem 10.6: the extra locus against the kernel of the specialization")
    print(bar)
    print("  %-9s %-4s %-10s %-8s %-9s %-8s" % ("free part", "t", "solutions", "cores", "[thm:extra]", "OTHER"))
    print("  " + "-" * 52)
    tally = {}
    for name, mk in FREE:
        for t in range(2, 9):
            cores = extras = other = 0
            for lam in two_row_shapes(3 * t + 4):
                agree = True
                for (z, w) in PTS:
                    Pf = mk(z, w)
                    Pa = (lambda k, Pf=Pf: F(t if k % t == 0 else 0) + Pf(k))
                    if abs(schur(lam, Pa)) != abs(schur(lam, Pf)):
                        agree = False
                        break
                if not agree:
                    continue
                if is_core(lam, t):
                    cores += 1
                elif is_extra(lam, t):
                    extras += 1
                else:
                    other += 1
            tally[(name, t)] = (cores, extras, other)
            print("  %-9s %-4d %-10d %-8d %-9d %-8d"
                  % (name, t, cores + extras + other, cores, extras, other))
        print("  " + "-" * 52)

    print()
    print(bar)
    print("WHAT THE TABLE SAYS")
    print(bar)
    zz = [tally[("(z,-z)", t)][2] for t in range(2, 9)]
    print("  (z,1/z)  extra beyond cores and [thm:extra]: %s"
          % [tally[("(z,1/z)", t)][2] for t in range(2, 9)])
    print("  (z,-z)   extra beyond cores and [thm:extra]: %s" % zz)
    print("  (z,-z) at t=2: %d further solutions" % tally[("(z,-z)", 2)][2])
    print("  (z,-z) at t=8: %d further solutions" % tally[("(z,-z)", 8)][2])
    print("  the reciprocal pair is the collapsing specialization of smallest kernel.")

    print()
    print(bar)
    print("CONTROLS")
    print(bar)
    free_clean = all(tally[("(z,w)", t)][1] == 0 and tally[("(z,w)", t)][2] == 0
                     for t in range(2, 9))
    print("  (z,w) free returns the t-cores and nothing else : %s   (Remark 5.3's control)"
          % ("yes" if free_clean else "*** NO -- this file is wrong ***"))
    recip_clean = all(tally[("(z,1/z)", t)][2] == 0 for t in range(2, 9))
    print("  (z,1/z) returns cores + [thm:extra] and nothing else : %s"
          % ("yes" if recip_clean else "*** NO ***"))
    counts = [tally[("(z,1/z)", t)][1] for t in range(2, 9)]
    print("  and its extras, t = 2..8 : %s   (t/2 at even t, none at odd)" % counts)
    zz_nonzero = any(v > 0 for v in zz)
    print("  (z,-z) does acquire further solutions : %s   (or the comparison is vacuous)"
          % ("yes" if zz_nonzero else "*** NO ***"))
    problems = (0 if free_clean else 1) + (0 if recip_clean else 1) + (0 if zz_nonzero else 1)
    print()
    print("TOTAL problems: %d" % problems)


if __name__ == "__main__":
    main()
