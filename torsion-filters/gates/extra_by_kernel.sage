# extra_by_kernel.sage
#
# Our Theorem 7.1: at the free part F = (z, 1/z), the independence criterion of [AK25, Thm 5.3]
# acquires exactly one extra family beyond the t-cores, for l(lambda) <= 2, and it is indexed by
# the order-two element r_B - r_A = t/2.
#
# Remark 7.2 explains WHY: on the reciprocal locus s_mu(z,1/z) = chi_{mu_1-mu_2}(z) depends only on
# mu_1 - mu_2, so the specialization has a KERNEL on the span of the universal characters, and the
# extra family measures it.  Problem 10.6 then asks whether the extra locus can always be read off
# that kernel.
#
# A GUESS OF MINE, KILLED BEFORE IT WAS RUN.  "The order-two element" invited a generalization
# indexed by the divisors of t: one extra family per d | t, ours being d = 2.  But the d-analogue
# of the free part would be (z, zeta_d z, ...), and at d = 2 that is (z, -z), which is NOT our
# alphabet.  The analogy breaks at its own base case, so it was never a conjecture -- it was a
# pattern read off a word.  What follows is the well-posed question instead.
#
# THE TEST.  Take several free parts F of size two, each with a different kernel, and compare the
# solution sets of  s_lambda(mu_t union F) = +- s_lambda(F)  for l(lambda) <= 2:
#
#     F = (z, 1/z)     inversion-closed;  s_mu(F) sees only mu_1 - mu_2      <- Theorem 7.1
#     F = (z, -z)      a zeta_2-orbit;    s_mu(F) is 0 or +- z^{|mu|}
#     F = (z, w)       free;              no kernel                          <- the control
#
# If the extra locus is a function of the kernel and not of the alphabet, the two collapsing
# specializations should behave alike in kind, and the free one should give the cores and nothing
# else.  That last is Remark 7.2's own control, so it must come out clean or the file is wrong.
#
# Authors: Carles Marin, Claude (AI assistant).

import sys
sys.stdout.reconfigure(line_buffering=True)


def parts_two(maxpart):
    out = [()]
    for a in range(0, maxpart + 1):
        out.append((a,))
        for b in range(1, a + 1):
            out.append((a, b))
    return [p for p in out if p == () or p[0] > 0]


def schur_from_p(lam, P):
    L = [k for k in lam if k > 0]
    n = len(L)
    if n == 0:
        return QQ(1)
    M = max(L) + n
    h = [QQ(1)]
    for m in range(1, M + 1):
        h.append(sum(P(i) * h[m - i] for i in range(1, m + 1)) / m)
    H = lambda d: QQ(0) if d < 0 else h[d]
    return matrix(QQ, [[H(L[i] - (i + 1) + (j + 1)) for j in range(n)] for i in range(n)]).det()


def is_core(lam, t):
    N = t + 2
    L = list(lam) + [0] * (N - len(lam))
    b = set(L[j] + N - (j + 1) for j in range(N))
    return not any(v >= t and (v - t) not in b for v in b)


def extra_family(lam, t):
    """the family of Theorem 7.1"""
    L = list(lam) + [0, 0]
    return (t % 2 == 0 and len([k for k in lam if k > 0]) == 2
            and L[0] == L[1] + 3 * t // 2 - 1 and t // 2 <= L[1] <= t - 1)


ZS = [(QQ(3) / 2, QQ(5) / 3), (QQ(7) / 4, QQ(2)), (QQ(5) / 2, QQ(3))]

FREE = {
    "(z, 1/z)": lambda z, w: (lambda k: z ** k + 1 / z ** k),
    "(z, -z)": lambda z, w: (lambda k: z ** k + (-z) ** k),
    "(z, w)  free": lambda z, w: (lambda k: z ** k + w ** k),
}

print("=" * 96)
print("Which free parts make the independence criterion acquire extra solutions?")
print("=" * 96)
print("%-14s %-4s %-9s %-8s %-8s %-9s %-22s"
      % ("free part", "t", "solutions", "cores", "Thm 7.1", "OTHER", "the others"))
print("-" * 96)

for name, mk in FREE.items():
    for t in range(2, 9):
        cores = extras = other = 0
        others = []
        for lam in parts_two(3 * t + 4):
            ok = True
            for (z, w) in ZS:
                Pf = mk(z, w)
                Pa = (lambda k, Pf=Pf: QQ(t if k % t == 0 else 0) + Pf(k))
                a = schur_from_p(lam, Pa)
                b = schur_from_p(lam, Pf)
                if abs(a) != abs(b):
                    ok = False
                    break
            if not ok:
                continue
            if is_core(lam, t):
                cores += 1
            elif extra_family(lam, t):
                extras += 1
            else:
                other += 1
                if len(others) < 3:
                    others.append(lam)
        print("%-14s %-4d %-9d %-8d %-8d %-9d %-22s"
              % (name, t, cores + extras + other, cores, extras, other, str(others) if others else ""))
    print("-" * 96)

print()
print("=" * 96)
print("READING")
print("=" * 96)
print("  CONTROL: the free row must be cores only (Thm 7.1 = 0, OTHER = 0) at every t.")
print("           That is Remark 7.2's own control; if it fails, this file is wrong.")
print("  (z,1/z): cores + the Theorem 7.1 family, OTHER = 0, extras only at even t.")
print("  (z,-z) : whatever it is, it is a different kernel.  If it also produces a family")
print("           indexed by a distinguished element of Z/t, Problem 10.6 has a pattern;")
print("           if it produces nothing or something unindexable, it does not.")
