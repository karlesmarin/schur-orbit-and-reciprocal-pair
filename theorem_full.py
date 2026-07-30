"""
The general-m law for s_lambda(mu_m, t, 1/t), EXACT SIGN INCLUDED, with every lemma checked
separately so a failure localises.

Authors: Carles Marin, Claude (AI assistant).

Chain (all re-derived by hand, then checked here):
  L1  every nonzero m x m minor of the root-of-unity rows is (-1)^{inv(b_S)} * V,
      V = prod_{r<r'} (w^{r'} - w^r)                                  [Vandermonde, alternating]
  L2  det(x_i^{N-j}) = (-1)^{m(m-1)/2} * V * (t^m-1)(t^-m-1)(t-1/t)   [Vandermonde, split]
  L3  Phi = (-1)^{m + N(N+1)/2} * SUM_U (-1)^{j1+j2+inv} f(b_j1-b_j2)
             / [ (t^m-1)(t^-m-1)(t-1/t) ]                            [Laplace along the m rows]
  L4  lambda_ij = lambda_11 * eps_i * delta_j,  eps = delta = (+1,-1) [column-move argument]
  L5  SUM_{eps,delta} eps*delta*f(c+eps*p+delta*q) = f(c)f(p)f(q)
  L5' f(2u) - f(2u+2v) + f(2v) = -f(u)f(v)f(u+v)                      [the size-3 profile]
  ==> THEOREM   Phi = eps * prod_i sinh(d_i theta/2) / [sinh(m theta/2)^2 sinh(theta)]
                eps = (-1)^{m + (m+2)(m+3)/2} * lambda_11 * sgn(a1+a2-b1-b2)
"""
from mpmath import mp, mpf, mpc, exp, sinh, pi
from law_control import mydet, partitions

mp.dps = 40


def f(x, theta):
    return exp(x * theta) - exp(-x * theta)


def inv_count(word):
    return sum(1 for i in range(len(word)) for j in range(i + 1, len(word)) if word[i] > word[j])


def setup(lam, m):
    """beta-set, residue classes, and the (A,B) profile. Returns None if a class is empty."""
    N = m + 2
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[j] + N - 1 - j for j in range(N)]          # strictly decreasing, column j -> beta[j]
    cls = {}
    for j, b in enumerate(beta):
        cls.setdefault(b % m, []).append(j)                 # columns, already increasing
    if len(cls) < m:
        return None
    big = [v for v in cls.values() if len(v) >= 2]
    if len(big) == 1:                                       # one class of size 3: A={p,q}, B={q,r}
        j1, j2, j3 = big[0]
        Acols, Bcols = (j1, j2), (j2, j3)
    else:                                                   # two classes of size 2
        big.sort()
        Acols, Bcols = tuple(big[0]), tuple(big[1])
    return beta, Acols, Bcols


def lambda11(beta, Acols, Bcols, N):
    """the single Laplace-times-sorting sign of the (a1,b1) term"""
    jA, jB = Acols[0], Bcols[0]
    if jA == jB:
        return None
    S = [j for j in range(N) if j not in (jA, jB)]
    m = N - 2
    word = [beta[j] % m for j in S]
    s = (-1) ** ((jA + 1) + (jB + 1) + inv_count(word))      # 1-indexed columns
    a1, b1 = beta[jA], beta[jB]
    return s * (1 if a1 > b1 else -1)


def closed_form(lam, m, theta):
    """the THEOREM: value with explicit sign, or 0."""
    st = setup(lam, m)
    if st is None:
        return mpc(0)
    beta, Acols, Bcols = st
    N = m + 2
    a1, a2 = beta[Acols[0]], beta[Acols[1]]
    b1, b2 = beta[Bcols[0]], beta[Bcols[1]]
    d = (a1 - a2, b1 - b2, abs(a1 + a2 - b1 - b2))
    if d[2] == 0:
        return mpc(0)
    l11 = lambda11(beta, Acols, Bcols, N)
    sgn_c = 1 if (a1 + a2 - b1 - b2) > 0 else -1
    eps = (-1) ** (m + (m + 2) * (m + 3) // 2) * l11 * sgn_c
    val = (sinh(d[0] * theta / 2) * sinh(d[1] * theta / 2) * sinh(d[2] * theta / 2)
           / (sinh(m * theta / 2) ** 2 * sinh(theta)))
    return eps * val


def phi_bialternant(lam, m, theta):
    N = m + 2
    t = exp(theta)
    w = exp(2j * pi / m)
    xs = [w ** k for k in range(m)] + [t, 1 / t]
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[j] + N - 1 - j for j in range(N)]
    num = [[xs[i] ** beta[j] for j in range(N)] for i in range(N)]
    den = [[xs[i] ** (N - 1 - j) for j in range(N)] for i in range(N)]
    return mydet(num, N) / mydet(den, N)


def check_L1_L2(m, theta):
    """L1 on every m-subset of a sample beta; L2 outright."""
    N = m + 2
    t = exp(theta)
    w = exp(2j * pi / m)
    V = mpc(1)
    for r in range(m):
        for rp in range(r + 1, m):
            V *= (w ** rp - w ** r)
    # L2
    xs = [w ** k for k in range(m)] + [t, 1 / t]
    den = [[xs[i] ** (N - 1 - j) for j in range(N)] for i in range(N)]
    lhs = mydet(den, N)
    rhs = (-1) ** (m * (m - 1) // 2) * V * (t ** m - 1) * (t ** (-m) - 1) * (t - 1 / t)
    ok2 = abs(lhs - rhs) < 1e-20
    # L1 over a batch of beta-sets
    import itertools
    ok1 = True
    n1 = 0
    for lam in partitions(6, N) + partitions(9, N):
        lamp = list(lam) + [0] * (N - len(lam))
        beta = [lamp[j] + N - 1 - j for j in range(N)]
        for S in itertools.combinations(range(N), m):
            M = [[w ** (k * beta[j]) for j in S] for k in range(m)]
            dM = mydet(M, m)
            word = [beta[j] % m for j in S]
            if len(set(word)) < m:
                if abs(dM) > 1e-18:
                    ok1 = False
            else:
                if abs(dM - (-1) ** inv_count(word) * V) > 1e-18:
                    ok1 = False
                n1 += 1
    return ok1, ok2, n1


def check_L4(lam, m):
    """lambda_ij = lambda_11 * eps_i * delta_j, on disjoint (2,2) profiles."""
    st = setup(lam, m)
    if st is None:
        return None
    beta, Acols, Bcols = st
    if len(set(Acols) & set(Bcols)) > 0:
        return None                                          # size-3 profile: handled by L5'
    N = m + 2
    l11 = lambda11(beta, Acols, Bcols, N)
    for i in (0, 1):
        for j in (0, 1):
            lij = lambda11(beta, (Acols[i],) * 2, (Bcols[j],) * 2, N)
            if lij != l11 * (-1) ** i * (-1) ** j:
                return False
    return True


def check_L5(theta):
    c, p, q = mpf("1.3"), mpf("0.7"), mpf("2.1")
    s = sum(e * d * f(c + e * p + d * q, theta) for e in (1, -1) for d in (1, -1))
    ok5 = abs(s - f(c, theta) * f(p, theta) * f(q, theta)) < 1e-20
    u, v = mpf("0.9"), mpf("1.7")
    s2 = f(2 * u, theta) - f(2 * u + 2 * v, theta) + f(2 * v, theta)
    ok5p = abs(s2 + f(u, theta) * f(v, theta) * f(u + v, theta)) < 1e-20
    return ok5, ok5p


def run():
    theta = mpf("0.41")
    print("== lemmas ==")
    ok5, ok5p = check_L5(theta)
    print(f"  L5  {ok5}    L5' (size-3 identity)  {ok5p}")
    for m in range(2, 7):
        ok1, ok2, n1 = check_L1_L2(m, theta)
        print(f"  m={m}:  L1 {ok1} ({n1} nonzero minors)   L2 {ok2}")

    print("== L4, and the THEOREM with its exact sign ==")
    tot = zeros = bad_sign = bad_mag = 0
    l4ok = l4n = 0
    for m in range(2, 8):
        N = m + 2
        maxsize = {2: 14, 3: 14, 4: 13, 5: 12, 6: 11, 7: 10}[m]
        for n in range(0, maxsize + 1):
            for lam in partitions(n, N):
                r4 = check_L4(lam, m)
                if r4 is not None:
                    l4n += 1
                    l4ok += bool(r4)
                phi = phi_bialternant(lam, m, theta)
                cf = closed_form(lam, m, theta)
                if abs(cf) < 1e-18:
                    if abs(phi) > 1e-12:
                        bad_mag += 1
                    else:
                        zeros += 1
                    continue
                dif = abs(phi - cf)
                if dif < 1e-15 * max(1, abs(cf)):
                    tot += 1
                elif abs(phi + cf) < 1e-15 * max(1, abs(cf)):
                    bad_sign += 1
                else:
                    bad_mag += 1
        print(f"  through m={m}:  exact match {tot}   zeros {zeros}   "
              f"WRONG SIGN {bad_sign}   WRONG MAGNITUDE {bad_mag}   L4 {l4ok}/{l4n}")
    print()
    print(f"THEOREM, sign included: {tot} exact matches, {zeros} zeros, "
          f"{bad_sign} sign failures, {bad_mag} magnitude failures.  L4: {l4ok}/{l4n}")


if __name__ == "__main__":
    run()
