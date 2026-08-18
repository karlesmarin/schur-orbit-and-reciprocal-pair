# sieving_1008.sage
#
# Problem 10.8 asks: is there a cyclic action on SSYT(lambda, t+2) and a q-analogue of Phi_t whose
# value at a t-th root of unity counts fixed points REFINED BY m_{t+1} - m_{t+2}?
#
# The paper poses it without an instrument.  There is one: Alexandersson-Amini [arXiv:1804.01447,
# Discrete Math. 2019] give a NECESSARY AND SUFFICIENT criterion for a cyclic action realising a CSP
# to exist, without constructing the action.  So the existence half of 10.8 is DECIDABLE.
#
# The candidate q-analogue is forced by the alphabet itself:
#
#     F(q,z) = s_lambda(1, q, ..., q^{t-1}, z, z^{-1}),          F(omega, z) = Phi_t(lambda; z),
#
# and the z-grading IS the refinement the problem asks for, since
#
#     s_lambda(x_1..x_t, z, 1/z) = sum_T x^{m_1..m_t} z^{m_{t+1} - m_{t+2}}.
#
# So write F(q,z) = sum_k f_k(q) z^k.  Then f_k(1) = #SSYT_k(lambda, t+2), the tableaux with
# m_{t+1} - m_{t+2} = k, and the refined statement of 10.8 is: for every k,
#
#     ( SSYT_k(lambda, t+2),  C_t,  f_k(q) )   is a CSP.
#
# One CSP per k.  Each is decidable by the AA criterion, which is what this file runs.
#
# NORMALISATION.  A CSP polynomial is only defined up to a power of q (RSW normalise by q^{-kappa}).
# Multiplying by q^a multiplies f(omega^d) by omega^{ad}, so the verdict depends on the shift.  We do
# NOT choose one: we test every shift a = 0..t-1 after stripping the lowest power of q, and report
# WHICH shifts work.  Reporting the shift is part of the answer -- if a single rule for a works
# across all (lambda, t, k), that rule is the kappa of the statement.
#
# CONTROLS, all of which must be able to fail:
#   C1  f_k(1) must equal the number of SSYT with that content difference, counted independently
#       by walking the tableaux.  Ties the polynomial to the set it is supposed to sieve.
#   C2  a deliberately WRONG refinement -- grading by m_{t+1} + m_{t+2} instead of the difference --
#       must not pass at the same rate.  If it passes just as often, the criterion is not seeing
#       our refinement and this whole file measures nothing.
#   C3  sum_k f_k(1) = #SSYT(lambda, t+2) = s_lambda(1,...,1).
#
# Run:  sage sieving_1008.sage        (writes everything to stdout; archive it)

R = PolynomialRing(ZZ, 'q'); q = R.gen()


def wt(T, N):
    """T.weight() padded to length N -- Sage truncates at the largest entry actually used."""
    w = list(T.weight())
    return w + [0] * (N - len(w))


def multiplicities(f, t):
    """the c_j of the Alexandersson-Amini criterion, by triangularity over the divisors of t.
    Returns dict d -> c_d, or None if some c_d is not a nonnegative integer."""
    K = CyclotomicField(t); w = K.gen()
    c = {}
    for d in sorted(ZZ(t).divisors()):
        val = f(w ** d)
        if val not in QQ:
            return None
        rem = QQ(val) - sum(j * c[j] for j in c if d % j == 0)
        if rem % d != 0:
            return None
        c[d] = ZZ(rem / d)
        if c[d] < 0:
            return None
    return c


def strip(f):
    """divide out the lowest power of q"""
    if f == 0:
        return f
    return R(f / q ** f.valuation())


def graded_pieces(lam, t, wrong=False):
    """f_k(q) for F(q,z) = s_lambda(1,q,..,q^{t-1}, z, 1/z), as a dict k -> polynomial in q.

    wrong=True grades by m_{t+1} + m_{t+2} instead of the difference -- control C2."""
    N = t + 2
    out = {}
    total = 0
    for T in SemistandardTableaux(lam, max_entry=N):
        w = wt(T, N)
        e = sum(w[i] * i for i in range(t))  # exponent of q from the first t letters
        k = (w[t] + w[t + 1]) if wrong else (w[t] - w[t + 1])
        out[k] = out.get(k, R(0)) + q ** e
        total += 1
    return out, total


def report(lam, t, wrong=False, verbose=True):
    pieces, total = graded_pieces(lam, t, wrong)
    N = t + 2
    # ---- control C3
    s_at_1 = SemistandardTableaux(lam, max_entry=N).cardinality()
    assert total == s_at_1, "C3 FAILED: %s vs %s" % (total, s_at_1)
    rows = []
    for k in sorted(pieces):
        f = strip(pieces[k])
        # ---- control C1: f_k(1) counts the tableaux in that graded piece
        n_k = sum(1 for T in SemistandardTableaux(lam, max_entry=N)
                  if ((wt(T, N)[t] + wt(T, N)[t + 1]) if wrong
                      else (wt(T, N)[t] - wt(T, N)[t + 1])) == k)
        assert f(1) == n_k, "C1 FAILED at k=%s: %s vs %s" % (k, f(1), n_k)
        good = [a for a in range(t) if multiplicities(R(q ** a * f), t) is not None]
        rows.append((k, n_k, good))
    if verbose:
        tag = "  [WRONG GRADING m_{t+1}+m_{t+2}]" if wrong else ""
        print("lambda = %-12s t = %d   N = %d   #SSYT = %d%s" % (lam, t, N, total, tag))
        for k, n_k, good in rows:
            print("    k = %-4d  #SSYT_k = %-6d  shifts a that pass: %s"
                  % (k, n_k, good if good else "NONE"))
        print()
    return rows


def sweep(t, maxsize, wrong=False):
    N = t + 2
    allrows = []
    for n in range(1, maxsize + 1):
        for lam in Partitions(n, max_length=N):
            allrows.extend([(list(lam), r) for r in report(list(lam), t, wrong, verbose=False)])
    npieces = len(allrows)
    passing = sum(1 for _, r in allrows if r[2])
    always0 = sum(1 for _, r in allrows if 0 in r[2])
    print("t = %d, |lambda| <= %d%s" % (t, maxsize, "   [WRONG GRADING]" if wrong else ""))
    print("    %d graded pieces, %d admit SOME shift (%.1f%%), %d admit a = 0"
          % (npieces, passing, 100.0 * passing / npieces, always0))
    bad = [(l, r) for l, r in allrows if not r[2]]
    for l, r in bad[:6]:
        print("        no shift: lambda = %-12s k = %d  #SSYT_k = %d" % (l, r[0], r[1]))
    if len(bad) > 6:
        print("        ... and %d more" % (len(bad) - 6))
    print()
    return npieces, passing


if __name__ == '__main__':
    print("=" * 96)
    print("Problem 10.8 by the Alexandersson-Amini criterion: one CSP per z-degree")
    print("=" * 96)
    print()
    report([2, 1], 2)
    report([3, 1], 3)
    print("=" * 96)
    print("SWEEPS -- the real grading")
    print("=" * 96)
    print()
    for t in [2, 3, 4]:
        sweep(t, 6 if t < 4 else 5)
    print("=" * 96)
    print("CONTROL C2 -- the WRONG grading, m_{t+1} + m_{t+2}.")
    print("If it passes at the same rate, this file measures nothing.")
    print("=" * 96)
    print()
    for t in [2, 3, 4]:
        sweep(t, 6 if t < 4 else 5, wrong=True)
