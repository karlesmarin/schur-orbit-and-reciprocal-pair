"""sign_one_object.py -- the sorting sign of Theorem 3.1 is the sign that makes Section 6 alternate.

Section 2 attaches to a partition the permutation sigma that sorts beta(lambda) by residue class,
the classes increasing and the values decreasing inside each class.  Proposition 3.11 puts sgn(sigma)
inside eps_lambda.  Section 6 separately proves that the terms sigma_m = s_{lambda/(m+k,m)}(1,-1)
alternate along each run.  This file checks that those are the same object:

  (1)  s_{lambda/mu}(1,-1) = sgn(sigma_lambda) sgn(sigma_mu) when core_2(lambda) = core_2(mu) and
       both components of the skew 2-quotient are horizontal strips, and 0 otherwise;

  (2)  consequently Proposition 6.4(ii) and (iii) are one statement about sigma alone: as
       nu = (m+k,m) moves to (m+k+1,m+1), sgn(sigma_nu) reverses for k even and is constant for
       k odd.

The identity behind (1) is Macdonald's ribbon sign [Mac95, Ex. I.1.8] written as a sorting sign;
it is Lemma 6.3 of [KumariThesis].  Nothing here is taken from either -- both statements are
recomputed from the definitions, in exact rational arithmetic, through Jacobi-Trudi with the
power sums of the alphabet (1,-1).

CONTROLS.  Four, and each must be reported as refuted.  Without them a bug that made every
comparison vacuous would leave this file printing all-clear.

Authors: Carles Marin, Claude (AI assistant).
"""
from fractions import Fraction as F


def parts(n, maxp=None, maxlen=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for k in range(min(n, maxp), 0, -1):
        for rest in parts(n - k, k, None if maxlen is None else maxlen - 1):
            yield (k,) + rest


def all_parts(nmax, maxlen=None):
    for n in range(nmax + 1):
        for l in parts(n, maxlen=maxlen):
            yield l


def betaset(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - (j + 1) for j in range(N)]


def sgn_sigma(lam, N, t=2):
    """the sign of the permutation of Section 2, at modulus t"""
    w = [v % t for v in betaset(lam, N)]
    return (-1) ** sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])


def profile(lam, N, t=2):
    return [sum(1 for v in betaset(lam, N) if v % t == i) for i in range(t)]


def quotient2(lam, N):
    b = betaset(lam, N)
    out = []
    for i in (0, 1):
        a = sorted([(v - i) // 2 for v in b if v % 2 == i], reverse=True)
        out.append(tuple(a[j] - (len(a) - 1 - j) for j in range(len(a))))
    return out


def is_hstrip(nu, kap):
    nu = list(nu)
    kap = list(kap) + [0] * (len(nu) - len(kap))
    if len(kap) > len(nu):
        if any(kap[i] for i in range(len(nu), len(kap))):
            return False
        kap = kap[:len(nu)]
    if any(kap[i] > nu[i] for i in range(len(nu))):
        return False
    return all(nu[i + 1] <= kap[i] for i in range(len(nu) - 1))


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


def skew_at_pm1(lam, mu):
    """s_{lam/mu}(1,-1) by Jacobi-Trudi; h_d(1,-1) = 1 for d even >= 0 and 0 otherwise"""
    L = list(lam)
    n = len(L)
    if n == 0:
        return F(1)
    M = list(mu) + [0] * (n - len(mu))
    H = lambda d: F(1) if (d >= 0 and d % 2 == 0) else F(0)
    return det_frac([[H(L[i] - M[j] - (i + 1) + (j + 1)) for j in range(n)] for i in range(n)])


def main():
    bar = "=" * 74
    print(bar)
    print("(1)  s_{lam/mu}(1,-1) = sgn(sigma_lam) sgn(sigma_mu) . [cores agree, both strips]")
    print(bar)
    grand = grandbad = 0
    for N in (4, 5, 6):
        tot = bad = nonzero = 0
        for lam in all_parts(12, maxlen=N):
            for mu in all_parts(sum(lam), maxlen=N):
                L = list(lam) + [0] * (N - len(lam))
                M = list(mu) + [0] * (N - len(mu))
                if len(mu) > len(lam) or any(M[i] > L[i] for i in range(N)):
                    continue
                tot += 1
                lhs = skew_at_pm1(L, M)
                ok = (profile(lam, N) == profile(mu, N))
                if ok:
                    ql, qm = quotient2(lam, N), quotient2(mu, N)
                    ok = is_hstrip(ql[0], qm[0]) and is_hstrip(ql[1], qm[1])
                rhs = sgn_sigma(lam, N) * sgn_sigma(mu, N) if ok else 0
                if lhs != rhs:
                    bad += 1
                if lhs != 0:
                    nonzero += 1
        grand += tot
        grandbad += bad
        print("  N=%d : %5d skew pairs, %4d of them nonzero, %d disagreements"
              % (N, tot, nonzero, bad))
    print("  TOTAL: %d skew pairs, %d disagreements" % (grand, grandbad))

    print()
    print(bar)
    print("(2)  Proposition 6.4(ii)+(iii) as a statement about sigma alone")
    print(bar)
    bad2 = tot2 = 0
    for N in (4, 5, 6):
        for k in range(0, 8):
            for m in range(0, 10):
                a = sgn_sigma((m + k, m), N)
                b = sgn_sigma((m + k + 1, m + 1), N)
                tot2 += 1
                if k % 2 == 0 and b != -a:
                    bad2 += 1
                if k % 2 == 1 and b != a:
                    bad2 += 1
    print("  N=4,5,6, k<=7, m<=9 : %d steps, %d disagreements" % (tot2, bad2))
    print("     k even -> sgn(sigma_nu) reverses ;  k odd -> sgn(sigma_nu) is constant")

    print()
    print(bar)
    print("(3)  and those two reproduce Proposition 6.4 on real lambda")
    print(bar)
    N = 4
    bi = bii = biii = shapes = 0
    for lam in all_parts(14, maxlen=N):
        for k in range(0, 8):
            sig = []
            for m in range(0, 10):
                nu = (m + k, m) if k + m > 0 else ()
                L = list(lam) + [0] * (N - len(lam))
                M = list(nu) + [0] * (N - len(nu))
                if len(M) > N or any(M[i] > L[i] for i in range(N)):
                    sig.append(0)
                    continue
                sig.append(skew_at_pm1(L, M))
            shapes += 1
            if any(v not in (0, 1, -1) for v in sig):
                bi += 1
            if k % 2 == 1 and any(sig[i] and sig[i + 1] for i in range(len(sig) - 1)):
                bii += 1
            if any(sig[i] and sig[i + 1] and sig[i + 1] != -sig[i] for i in range(len(sig) - 1)):
                biii += 1
    print("  |lam| <= 14, <= 4 rows, k <= 7 : %d (lam,k) families" % shapes)
    print("     (i)   every sigma_m in {0,+-1}          : %d fail" % bi)
    print("     (ii)  k odd, no two consecutive nonzero : %d fail" % bii)
    print("     (iii) consecutive nonzero alternate     : %d fail" % biii)

    print()
    print(bar)
    print("CONTROLS -- each must be refuted")
    print(bar)
    ctl = []
    ctl.append(("the strip condition is load-bearing",
                all(skew_at_pm1(list(l) + [0] * (4 - len(l)), list(m) + [0] * (4 - len(m)))
                    == sgn_sigma(l, 4) * sgn_sigma(m, 4)
                    for l in all_parts(8, maxlen=4) for m in all_parts(sum(l), maxlen=4)
                    if len(m) <= len(l)
                    and all((list(m) + [0] * 4)[i] <= (list(l) + [0] * 4)[i] for i in range(4)))))
    ctl.append(("the core condition is load-bearing",
                skew_at_pm1([1, 0, 0, 0], [0, 0, 0, 0]) == sgn_sigma((1,), 4) * sgn_sigma((), 4)))
    ctl.append(("sgn(sigma) is not identically +1",
                all(sgn_sigma(l, 4) == 1 for l in all_parts(6, maxlen=4))))
    ctl.append(("the alternation is not vacuous: some run has length >= 2",
                not any(skew_at_pm1([4, 4, 0, 0], [m, m, 0, 0]) != 0
                        and skew_at_pm1([4, 4, 0, 0], [m + 1, m + 1, 0, 0]) != 0
                        for m in range(0, 5))))
    problems = 0
    for name, held in ctl:
        print("  %-52s %s" % (name, "correctly refuted" if not held else "*** HELD ***"))
        if held:
            problems += 1

    print()
    print(bar)
    print("TOTAL problems: %d" % (grandbad + bad2 + bi + bii + biii + problems))


if __name__ == "__main__":
    main()
