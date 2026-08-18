# -*- coding: utf-8 -*-
# THE TOP-DEGREE COMPONENT, which factorises into two GL(r) alternants -- and the anatomy of the
# shapes where the single leading monomial was not enough.
#
# WHAT topterm_gate.sage SETTLED, AND WHAT IT KILLED.  The lex-leading monomial z^Delta of
# Phi_t = sum_P w(P) A(T_P) is computable by hand, and H2 -- "criterion fails => c_top != 0" -- is
# FALSE: 261 exceptions in 12424 criterion-failing shapes (t=4 r=2: 189, t=4 r=3: 56, t=6 r=2: 11,
# t=6 r=3: 5; none at t=4 r=1 or t=8 r=2).  A single monomial is too crude.
#
# THE STRUCTURE THAT REPLACES IT.  Grade by total degree sum_j |e_j| instead of by lex.  For a
# 2r-set T with values u_1 > .. > u_{2r}, every monomial of A(T) is prod_j z_j^{u_a - u_b} over a
# perfect matching, so its degree is at most
#     deg(T) = (u_1 + .. + u_r) - (u_{r+1} + .. + u_{2r}),
# with equality exactly when every pair joins the top half H = {u_1..u_r} to the bottom half
# L = {u_{r+1}..u_{2r}}.  Those matchings split into a choice of bijection H -> the z_j columns and
# a choice of bijection L -> the 1/z_j columns, so the top-degree part FACTORISES:
#
#     [A(T)]_top  =  +- a_H(z) * a_L(1/z),      a_X(z) = det( z_j^{x_i} )_{r x r},
#
# ONE SECTOR AT A TIME -- and K1 caught me here on the first run.  A pair {high, low} can also be
# oriented the other way round, so the degree-Dmax slice splits into 2^r orientation sectors indexed
# by the set J of flipped variables, and Phi is alternating under z_j -> 1/z_j, so the sectors are
# copies of one another and never mix (a monomial has e_j > 0 exactly for j not in J, and e_j is
# never 0 at top degree because the values are distinct).  So [Phi]_top = 0 iff the all-positive
# sector vanishes, and that sector is the two-alternant product above.  Hence with
# Dmax = max_P deg(T_P),
#
#     [Phi_t]_top  =  sum_{P : deg(T_P) = Dmax}  w(P) * [A(T_P)]_top
#
# is a signed combination of PRODUCTS OF TWO SCHUR FUNCTIONS, one in z and one in 1/z.  That is a
# far bigger obstruction than one coefficient, and it is still finite and explicit.
#
# THE QUESTION, stated so it can fail:
#   H3   criterion FAILS  =>  [Phi_t]_top  is not identically zero.
#
# AND THE ANATOMY.  For every shape where the lex-leading monomial cancelled, record whether it is
# concentric (condition (i) holds, (ii) fails) or not, how many excess classes it has, and how many
# of the two fixed classes lie in E.  If the exceptions sit on the (i)-holds-(ii)-fails boundary,
# that boundary is where the converse is actually hard.
#
# CONTROLS, each able to fail:
#   K1  acceptance, fatal, three ways: sum_j |e_j| over the surviving monomials of the FULL
#       expansion must have maximum exactly Dmax; the all-positive sector of that slice must equal
#       the factorised [Phi_t]_top monomial for monomial; and the whole slice must be exactly 2^r
#       times as big as the sector, or the sectors are not what I claim.
#   K2  criterion HOLDS => [Phi_t]_top = 0.  Forced by the theorem.
#   K3  non-vacuity: [Phi_t]_top != 0 must actually occur.
#   K4  H3 must be tested on the SAME shapes where H2 failed, or it proves nothing new.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

CONF = [(4, 1, 30), (4, 2, 24), (4, 3, 18), (6, 2, 18), (6, 3, 14), (8, 2, 16)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def perm_sign(q):
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
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    out = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        Pset = set(P)
        out.append((w, [beta[i] for i in range(N) if i not in Pset]))
    return out


def nested(T):
    n = len(T)
    return tuple(T[j] - T[n - 1 - j] for j in range(n // 2))          # T already decreasing


def dominant(e):
    return tuple(sorted((abs(x) for x in e), reverse=True))


def degree(T, r):
    return sum(T[:r]) - sum(T[r:])


def topdeg_dict(T, r):
    """the degree-maximal part of A(T), as {exponent tuple: coefficient}."""
    D = {}
    n = 2 * r
    for a in itertools.permutations(range(r)):
        for b in itertools.permutations(range(r)):
            q = [0] * n
            e = [0] * r
            for i in range(r):
                q[i] = 2 * a[i]
                e[a[i]] += T[i]
            for i in range(r):
                q[r + i] = 2 * b[i] + 1
                e[b[i]] -= T[r + i]
            k = tuple(e)
            D[k] = D.get(k, 0) + perm_sign(q)
    return {k: v for k, v in D.items() if v != 0}


def full_expansion(tm, r):
    n = 2 * r
    D = {}
    for w, T in tm:
        for q in itertools.permutations(range(n)):
            e = [0] * r
            for a in range(n):
                e[q[a] // 2] += (T[a] if q[a] % 2 == 0 else -T[a])
            k = tuple(e)
            D[k] = D.get(k, 0) + w * perm_sign(list(q))
    return {k: v for k, v in D.items() if v != 0}


def phi_top(tm, r):
    Dmax = max(degree(T, r) for _, T in tm)
    out = {}
    for w, T in tm:
        if degree(T, r) != Dmax:
            continue
        for k, v in topdeg_dict(T, r).items():
            out[k] = out.get(k, 0) + w * v
    return Dmax, {k: v for k, v in out.items() if v != 0}


def anatomy(beta, t):
    """(criterion, concentric, |E|, number of fixed classes lying in E, sorted class sizes)."""
    cl = classes(beta, t)
    if len(cl) < t:
        return None
    E = set(k for k in cl if len(cl[k]) >= 2)
    S = set(b for k in E for b in cl[k])
    C = min(S) + max(S)
    conc = set(C - b for b in S) == S
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    s = sum(1 for k in fixed if k in E)
    crit = conc and len(fixed) == 2 and s == 2
    return crit, conc, len(E), s, tuple(sorted((len(cl[k]) for k in E), reverse=True))


# ---------------------------------------------------------------- K1, acceptance -----------------
print("=" * 100)
print("K1  the factorised top-degree part against the full expansion")
print("=" * 100)
print("")
print("     t   r  shapes  deg bad  part bad  sec bad   Dmax ex.   monomials in [Phi]_top")
print("  " + "-" * 96)

k1fail = 0
for t, r, _ in CONF:
    N = t + 2 * r
    nsh = degbad = partbad = secbad = 0
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
            Dmax, Ptop = phi_top(tm, r)
            FE = full_expansion(tm, r)
            if FE:
                if max(sum(abs(x) for x in k) for k in FE) > Dmax:
                    degbad += 1
            full = {k: v for k, v in FE.items() if sum(abs(x) for x in k) == Dmax}
            pos = {k: v for k, v in full.items() if all(x > 0 for x in k)}
            if pos != Ptop:
                partbad += 1
            if len(full) != 2 ** r * len(pos):
                secbad += 1
            if ex is None:
                ex = (Dmax, len(Ptop))
        if nsh > 25:
            break
    k1fail += degbad + partbad + secbad
    print("  %4d %3d %7d %8d %9d %8d %11d %10d"
          % (t, r, min(nsh, 25), degbad, partbad, secbad, ex[0], ex[1]))

if k1fail:
    print("")
    print("  K1 FAILED -- the factorisation is wrong, nothing below means anything.")
    raise SystemExit(1)
print("")
print("  K1 PASS: Dmax really is the top degree, and the two-alternant factorisation reproduces")
print("  the degree-Dmax slice of the full expansion monomial for monomial.")

# ---------------------------------------------------------------- the sweep ----------------------
print("")
print("=" * 100)
print("H3  criterion FAILS => [Phi_t]_top != 0 ?      (K2: criterion HOLDS => [Phi_t]_top = 0)")
print("=" * 100)
print("")
print("     t   r  |lam|<=  shapes   crit holds  K2 bad | crit fails  top != 0   H3 EXCEPTIONS")
print("  " + "-" * 96)

BAD = []
H2BAD = []
TOT = [0, 0, 0, 0, 0]
for t, r, MAX in CONF:
    N = t + 2 * r
    nsh = ncrit = k2 = nfail = nnz = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            tm = terms(beta, t, r)
            if tm is None:
                continue
            an = anatomy(beta, t)
            nsh += 1
            _, Ptop = phi_top(tm, r)
            Delta = max(nested(T) for _, T in tm)
            ctop = None
            if an[0]:
                ncrit += 1
                if Ptop:
                    k2 += 1
            else:
                nfail += 1
                if Ptop:
                    nnz += 1
                else:
                    BAD.append((t, r, list(l), beta, an))
            # H2 recomputation (cheap): the coefficient of the lex-leading monomial
            n2 = 2 * r
            tot = 0
            for w, T in tm:
                used = [False] * n2
                q = [0] * n2

                def rec(u, w=w, T=T, used=used, q=q):
                    if u == r:
                        return perm_sign(q)
                    acc = 0
                    for a in range(n2):
                        if used[a]:
                            continue
                        for b in range(n2):
                            if used[b] or b == a or T[a] - T[b] != Delta[u]:
                                continue
                            used[a] = used[b] = True
                            q[a], q[b] = 2 * u, 2 * u + 1
                            acc += rec(u + 1)
                            used[a] = used[b] = False
                    return acc
                tot += w * rec(0)
            if (not an[0]) and tot == 0:
                H2BAD.append((t, r, list(l), an, bool(Ptop)))
    print("  %4d %3d %8d %7d %12d %7d | %10d %9d %15d"
          % (t, r, MAX, nsh, ncrit, k2, nfail, nnz, nfail - nnz))
    TOT[0] += nsh
    TOT[1] += ncrit
    TOT[2] += k2
    TOT[3] += nfail
    TOT[4] += nnz

print("")
print("  totals: criterion holds %d (K2 bad %d) | criterion fails %d, top-degree part survives %d"
      % (TOT[1], TOT[2], TOT[3], TOT[4]))
print("  K3 non-vacuity: %s      H3 exceptions: %d"
      % ("PASS" if TOT[4] > 0 else "VACUOUS", TOT[3] - TOT[4]))

# ---------------------------------------------------------------- anatomy ------------------------
print("")
print("=" * 100)
print("K4  the anatomy of the shapes where the LEX-LEADING monomial cancelled (H2 exceptions)")
print("=" * 100)
print("")
print("  %d such shapes.  Does the top-degree part rescue them?" % len(H2BAD))
print("")
rescued = sum(1 for x in H2BAD if x[4])
print("     rescued by [Phi]_top (top-degree survives): %d of %d" % (rescued, len(H2BAD)))
print("")
print("  breakdown by (concentric?, |E|, fixed classes in E):")
tab = {}
for t, r, lam, an, resc in H2BAD:
    key = (an[1], an[2], an[3])
    tab[key] = tab.get(key, [0, 0])
    tab[key][0] += 1
    tab[key][1] += 1 if resc else 0
print("")
print("     concentric   |E|   fixed-in-E   count   of which rescued")
print("  " + "-" * 96)
for key in sorted(tab):
    print("     %-12s %4d %12d %8d %17d"
          % ("yes" if key[0] else "no", key[1], key[2], tab[key][0], tab[key][1]))

if BAD:
    print("")
    print("  the H3 exceptions themselves (criterion fails AND the whole top degree cancels):")
    print("")
    print("     t   r  lambda                     beta                             conc  |E|  fixed")
    print("  " + "-" * 96)
    for t, r, lam, beta, an in BAD[:20]:
        print("  %4d %3d  %-26s %-32s %-5s %3d %6d"
              % (t, r, str(lam), str(beta), "yes" if an[1] else "no", an[2], an[3]))
    print("")
    print("  total H3 exceptions: %d" % len(BAD))

print("")
print("DONE")
