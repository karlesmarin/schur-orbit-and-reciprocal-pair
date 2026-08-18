# -*- coding: utf-8 -*-
# THE CONVERSE, attacked by its extremal term.  Does the leading monomial survive?
#
# WHERE THIS SITS.  (i)+(ii) => Phi_t = 0 is proved (the Laplace terms pair with ratio -1).  The
# converse Phi_t = 0 => (i)+(ii) is MEASURED (928 zeros, 0 exceptions) and contains conj:crit of
# arXiv:2608.09619.  The standard way to prove an alternating sum is nonzero is to exhibit ONE
# monomial with a unique preimage.  This gate asks whether that route is even open.
#
# THE OBJECT.  With beta_0 > .. > beta_{N-1}, N = t + 2r, and every residue class occupied,
# Laplace along the t root-of-unity columns keeps only the transversals P:
#
#     det M  =  V0 * (-1)^{t(t-1)/2} * sum_P  w(P) * A(T_P),
#     w(P)   =  (-1)^{sum_{i in P} i} * sgn(residues of P in row order),
#     T_P    =  the 2r leftover values, A(T) = det( col^{T_a} ) on the columns z_1,1/z_1,..,z_r,1/z_r.
#
# Every monomial of A(T) is prod_j z_j^{T_a - T_b} over a perfect matching of T into ordered pairs.
# For a fixed T the lex-largest dominant exponent is the NESTED matching
#     delta(T) = ( u_1-u_{2r}, u_2-u_{2r-1}, .., u_r-u_{r+1} ),   u_1 > .. > u_{2r} the values of T,
# because the top difference u_1-u_{2r} is attained only by pairing the extremes, then induct.  So
#     Delta(lambda) = lex-max over P of delta(T_P)
# is the leading dominant exponent of the whole sum, and
#     c_top = sum_P w(P) * [coefficient of z^Delta in A(T_P)]
# is a finite integer computable by hand.
#
# THE QUESTION, stated so it can fail:
#   H2   criterion FAILS  =>  c_top != 0.
# If H2 held, the converse -- hence conj:crit -- would reduce to a sign count over transversals.
# I expect it to break; the point is to find WHERE, and how far below Delta the survivor sits.
#
# CONTROLS, each able to fail:
#   K1  acceptance, fatal: the Laplace expansion must equal det M up to ONE constant, checked by
#       evaluating both over GF(p) at random z and across several shapes and several t.
#   K2  c_top from the fast recursion must equal the coefficient read off the full expansion.
#   K3  no monomial of the full expansion may have dominant exponent above Delta.
#   K4  criterion HOLDS => c_top = 0.  Forced by the theorem; a nonzero would kill the setup.
#   K5  non-vacuity: c_top != 0 must actually occur, or H2 says nothing.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

CONF = [(4, 1, 30), (4, 2, 24), (4, 3, 18), (6, 2, 18), (6, 3, 14), (8, 2, 16)]

L = lcm([c[0] for c in CONF])
p = next_prime(10 ** 9)
while (p - 1) % L != 0:
    p = next_prime(p)
F = GF(p)
G = F.multiplicative_generator()
print("field GF(%d); guard on t in %s -> PASS" % (p, sorted(set(c[0] for c in CONF))))


def zeta(t):
    z = G ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t))
    return z


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def perm_sign(q):
    """sign of the permutation a -> q[a]."""
    n = len(q)
    seen = [False] * n
    s = 1
    for i in range(n):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = q[j]
            ln += 1
        if ln % 2 == 0:
            s = -s
    return s


def classes(beta, t):
    cl = {}
    for b in beta:
        cl.setdefault(b % t, []).append(b)
    return cl


def terms(beta, t, r):
    """[(w(P), T values in row order)] over all transversals P."""
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    out = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        res = [beta[i] % t for i in P]
        w = perm_sign(res)
        if sum(P) % 2:
            w = -w
        rest = [i for i in range(N) if i not in set(P)]
        out.append((w, [beta[i] for i in rest]))
    return out


def nested(T):
    u = sorted(T, reverse=True)
    n = len(u)
    return tuple(u[j] - u[n - 1 - j] for j in range(n // 2))


def dominant(e):
    return tuple(sorted((abs(x) for x in e), reverse=True))


def coeff_of(T, d):
    """coefficient of prod_j z_j^{d_j} in A(T), by recursion over ordered pairs."""
    r = len(d)
    n = 2 * r
    used = [False] * n
    q = [0] * n

    def rec(u):
        if u == r:
            return perm_sign(q)
        tot = 0
        for a in range(n):
            if used[a]:
                continue
            for b in range(n):
                if used[b] or b == a:
                    continue
                if T[a] - T[b] != d[u]:
                    continue
                used[a] = used[b] = True
                q[a], q[b] = 2 * u, 2 * u + 1
                tot += rec(u + 1)
                used[a] = used[b] = False
        return tot

    return rec(0)


def full_expansion(tm, r):
    """the whole Laplace sum as {exponent tuple: integer coefficient}."""
    n = 2 * r
    D = {}
    for w, T in tm:
        for q in itertools.permutations(range(n)):
            e = [0] * r
            for a in range(n):
                u, sgnc = q[a] // 2, (1 if q[a] % 2 == 0 else -1)
                e[u] += sgnc * T[a]
            k = tuple(e)
            D[k] = D.get(k, 0) + w * perm_sign(list(q))
    return {k: v for k, v in D.items() if v != 0}


def criterion(beta, t):
    cl = classes(beta, t)
    if len(cl) < t:
        return None
    E = set(k for k in cl if len(cl[k]) >= 2)
    S = set(b for k in E for b in cl[k])
    if not S:
        return None
    C = min(S) + max(S)
    conc = set(C - b for b in S) == S
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    return conc and len(fixed) == 2 and all(k in E for k in fixed)


# ---------------------------------------------------------------- K1, acceptance -----------------
print("")
print("=" * 100)
print("K1  the Laplace expansion against the determinant, over GF(p) at random z")
print("=" * 100)
print("")
print("     t   r   lambda                  shapes  points  ratio constant   matches V0*(-1)^C(t,2)")
print("  " + "-" * 96)

k1fail = 0
for t, r, _ in CONF:
    N = t + 2 * r
    zt = zeta(t)
    V0 = prod([zt ** kp - zt ** k for k in range(t) for kp in range(k + 1, t)])
    want = V0 * (-1) ** (t * (t - 1) // 2)
    set_random_seed(97 + t + r)
    consts = set()
    npts = nsh = 0
    shown = None
    for size in range(0, 14):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            tm = terms(beta, t, r)
            if tm is None:
                continue
            nsh += 1
            if nsh > 6:
                break
            if shown is None:
                shown = list(l)
            for _ in range(3):
                zv = []
                while len(zv) < r:
                    x = F.random_element()
                    if x != 0 and x ** 2 != 1:
                        zv.append(x)
                cols = [zt ** k for k in range(t)]
                for u in range(r):
                    cols += [zv[u], 1 / zv[u]]
                M = matrix(F, N, N, lambda i, j: cols[j] ** beta[i])
                lhs = M.det()
                rhs = F(0)
                for w, T in tm:
                    A = matrix(F, 2 * r, 2 * r,
                               lambda a, b: (zv[b // 2] ** T[a]) if b % 2 == 0
                               else (zv[b // 2] ** (-T[a])))
                    rhs += w * A.det()
                npts += 1
                if rhs == 0:
                    if lhs != 0:
                        consts.add("BROKEN")
                else:
                    consts.add(lhs / rhs)
        if nsh > 6:
            break
    ok = (len(consts) == 1) and ("BROKEN" not in consts)
    match = ok and (list(consts)[0] == want)
    if not ok:
        k1fail += 1
    print("  %4d %3d   %-22s %6d %7d   %-15s %s"
          % (t, r, str(shown), min(nsh, 6), npts,
             "unique" if ok else "NOT CONSTANT", "yes" if match else "no"))

if k1fail:
    print("")
    print("  K1 FAILED -- the term signs are wrong, nothing below is meaningful.")
    raise SystemExit(1)
print("")
print("  K1 PASS: one constant per (t,r), and it is V0*(-1)^{t(t-1)/2}.  The relative signs are right.")

# ---------------------------------------------------------------- K2, K3 -------------------------
print("")
print("=" * 100)
print("K2  fast c_top against the full expansion    K3  nothing above Delta survives")
print("=" * 100)
print("")
print("     t   r  shapes  K2 bad  K3 bad   Delta example                 c_top")
print("  " + "-" * 96)

k23fail = 0
for t, r, _ in CONF:
    N = t + 2 * r
    nsh = k2b = k3b = 0
    ex = None
    for size in range(0, 13):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            tm = terms(beta, t, r)
            if tm is None:
                continue
            nsh += 1
            if nsh > 25:
                break
            Delta = max(nested(T) for _, T in tm)
            ct = sum(w * coeff_of(T, Delta) for w, T in tm)
            D = full_expansion(tm, r)
            if D.get(Delta, 0) != ct:
                k2b += 1
            for k in D:
                if dominant(k) > Delta:
                    k3b += 1
                    break
            if ex is None:
                ex = (Delta, ct)
        if nsh > 25:
            break
    k23fail += k2b + k3b
    print("  %4d %3d %7d %7d %7d   %-28s %6d"
          % (t, r, min(nsh, 25), k2b, k3b, str(ex[0]), ex[1]))

if k23fail:
    print("")
    print("  K2/K3 FAILED -- read the columns.")
    raise SystemExit(1)
print("")
print("  K2/K3 PASS: the recursion computes the right coefficient, and Delta really is the top.")

# ---------------------------------------------------------------- the sweep ----------------------
print("")
print("=" * 100)
print("H2  criterion FAILS => c_top != 0 ?     (K4: criterion HOLDS => c_top = 0, forced)")
print("=" * 100)
print("")
print("     t   r  |lam|<=  shapes   crit holds  K4 bad | crit fails  c_top != 0  H2 EXCEPTIONS")
print("  " + "-" * 96)

EXC = []
TOT = [0, 0, 0, 0, 0]
for t, r, MAX in CONF:
    N = t + 2 * r
    nsh = ncrit = k4 = nfail = nnz = nexc = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            tm = terms(beta, t, r)
            if tm is None:
                continue
            cr = criterion(beta, t)
            nsh += 1
            Delta = max(nested(T) for _, T in tm)
            ct = sum(w * coeff_of(T, Delta) for w, T in tm)
            if cr:
                ncrit += 1
                if ct != 0:
                    k4 += 1
            else:
                nfail += 1
                if ct != 0:
                    nnz += 1
                else:
                    nexc += 1
                    if len(EXC) < 400:
                        EXC.append((t, r, list(l), beta, Delta))
    print("  %4d %3d %8d %7d %12d %7d | %10d %11d %14d"
          % (t, r, MAX, nsh, ncrit, k4, nfail, nnz, nexc))
    TOT[0] += nsh
    TOT[1] += ncrit
    TOT[2] += k4
    TOT[3] += nfail
    TOT[4] += nnz

print("")
print("  totals: shapes %d | criterion holds %d, K4 bad %d | criterion fails %d, c_top != 0 in %d"
      % (TOT[0], TOT[1], TOT[2], TOT[3], TOT[4]))
print("  K5 non-vacuity: c_top != 0 occurs %d times -- %s"
      % (TOT[4], "PASS" if TOT[4] > 0 else "VACUOUS"))
print("  H2 exceptions (criterion fails but the leading term cancels): %d"
      % (TOT[3] - TOT[4]))

# ---------------------------------------------------------------- how deep -----------------------
if EXC:
    print("")
    print("=" * 100)
    print("THE EXCEPTIONS: how far below Delta does the first survivor sit?")
    print("=" * 100)
    print("")
    print("     t   r  lambda                     Delta            first surviving dominant  rank")
    print("  " + "-" * 96)
    shown = 0
    depth = {}
    for t, r, lam, beta, Delta in EXC:
        tm = terms(beta, t, r)
        D = full_expansion(tm, r)
        doms = sorted(set(dominant(k) for k in D), reverse=True)
        top = doms[0] if doms else None
        # rank of Delta among all dominant exponents that OCCUR in the raw term set
        levels = sorted(set(nested(T) for _, T in tm), reverse=True)
        rk = 1 + (levels.index(top) if top in levels else -1)
        depth[rk] = depth.get(rk, 0) + 1
        if shown < 14:
            print("  %4d %3d  %-26s %-16s %-24s %4s"
                  % (t, r, str(lam), str(Delta), str(top), str(rk) if rk else "off-list"))
            shown += 1
    print("")
    print("  distribution of the rank of the first surviving level (1 = Delta itself): %s"
          % (sorted(depth.items()),))

print("")
print("DONE")
