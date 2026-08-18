# kt_calibrate.sage -- DERIVE the Koike-Terada determinant instead of recalling it.
#
# Twice now a remembered formula has produced a table of numbers that meant nothing:
#   - sp evaluated on r values instead of the 2r eigenvalues  (C0 caught it, 52/60)
#   - a wrong Jacobi-Trudi form for o                          (C0b caught it, 106/59)
# Both were caught by a control, but only after the fact.  So: do not recall the formula.  Search
# a family of candidate forms and keep the ones that agree with Sage's own o and sp bases.
#
# Candidate family, for kind in {sp, o}:
#     M_ij = h_{lam_i - i + j}  +  s * h_{lam_i - i - j + a}
# over s in {+1,-1} and a in {0,1,2,3}, optionally with the first column halved.
#
# The winner is then used to evaluate at a torus element in size l(lambda), which is what makes
# Example 4.9 of arXiv:2501.00275 (|lambda| = 55, l(lambda) = 11) reachable at all.
#
# Authors: Carles Marin, Claude (AI assistant).

import sys
sys.stdout.reconfigure(line_buffering=True)

Sym = SymmetricFunctions(QQ)
p, O, SP = Sym.p(), Sym.o(), Sym.sp()


def ev(f, P):
    out = 0
    for rho, c in p(f):
        term = c
        for k in rho:
            term *= P(k)
        out += term
    return out


def univ(basis, nu):
    nu = [k for k in nu if k > 0]
    return basis[nu] if nu else basis.one()


def h_from_P(P, M):
    h = [1]
    for m in range(1, M + 1):
        h.append(sum(P(i) * h[m - i] for i in range(1, m + 1)) / m)
    return h


def cand(lam, P, s, a, half_first):
    L = [k for k in lam if k > 0]
    l = len(L)
    if l == 0:
        return 1
    M = max(L) + l + 4
    h = h_from_P(P, M)
    H = lambda d: 0 if d < 0 or d > M else h[d]
    rows = []
    for i in range(l):
        row = []
        for j in range(l):
            v = H(L[i] - (i + 1) + (j + 1)) + s * H(L[i] - (i + 1) - (j + 1) + a)
            if half_first and j == 0:
                v = v / 2
            row.append(v)
        rows.append(row)
    return matrix(QQ, rows).det()


ZS = [[QQ(3) / 2], [QQ(3) / 2, QQ(5) / 3], [QQ(3) / 2, QQ(5) / 3, QQ(7) / 4]]

print("=" * 88)
print("Searching the candidate family against Sage's o and sp")
print("=" * 88)
print("%-6s %-4s %-4s %-12s %-10s" % ("kind", "s", "a", "half 1st col", "disagreements"))
print("-" * 88)
winners = {}
for kind, basis in (("sp", SP), ("o", O)):
    for s in (1, -1):
        for a in (0, 1, 2, 3):
            for half in (False, True):
                bad = tot = 0
                for zs in ZS:
                    r = len(zs)
                    P = lambda k: sum(z ** k + z ** (-k) for z in zs)
                    for size in range(0, 7):
                        for nu in Partitions(size, max_length=r):
                            tot += 1
                            if cand(list(nu), P, s, a, half) != ev(univ(basis, nu), P):
                                bad += 1
                if bad == 0:
                    winners.setdefault(kind, []).append((s, a, half))
                    print("%-6s %-4d %-4d %-12s %-10s   <== MATCHES on %d shapes"
                          % (kind, s, a, half, bad, tot))

print()
if 'sp' not in winners or 'o' not in winners:
    print("NO CANDIDATE MATCHES for %s -- the family is too small; do not proceed."
          % [k for k in ('sp', 'o') if k not in winners])
    raise SystemExit(1)
print("winners:  sp -> %s     o -> %s" % (winners['sp'], winners['o']))
s_sp, a_sp, h_sp = winners['sp'][0]
s_o, a_o, h_o = winners['o'][0]

print()
print("=" * 88)
print("INDEPENDENT CHECK of the winners: our Lemma 8.8, through the determinants only")
print("        o_nu(1,-1,z_1,1/z_1,...) = sp_nu(z_1,1/z_1,...)")
print("=" * 88)
bad = tot = 0
for zs in ZS:
    P_sp = lambda k: sum(z ** k + z ** (-k) for z in zs)
    P_o = lambda k: 1 + (-1) ** k + sum(z ** k + z ** (-k) for z in zs)
    for size in range(0, 8):
        for nu in Partitions(size):
            tot += 1
            if cand(list(nu), P_o, s_o, a_o, h_o) != cand(list(nu), P_sp, s_sp, a_sp, h_sp):
                bad += 1
print("    %d partitions, %d disagreements" % (tot, bad))
if bad:
    print("    *** the determinants do not reproduce Lemma 8.8 -- stop ***")
    raise SystemExit(1)
print("    the determinants reproduce Lemma 8.8 on their own.  Usable.")

print()
print("=" * 88)
print("EXAMPLE 4.9 of arXiv:2501.00275   (t=5, n=2, |mu|=44, |lambda|=55, l(lambda)=11)")
print("=" * 88)
mu49 = [9, 8, 6, 5, 4, 4, 3, 2, 2, 1]
lam49 = [10, 9, 7, 6, 5, 5, 4, 3, 3, 2, 1]
t49 = 5
Xs49 = [QQ(3) / 2, QQ(5) / 3]
P_W = lambda k: (t49 * sum(x ** k + x ** (-k) for x in Xs49)) if k % t49 == 0 else 0
P_A = lambda k: P_W(k) + 2
P_B = lambda k: P_W(k) + 1 + (-1) ** k
v_sp = cand(mu49, P_W, s_sp, a_sp, h_sp)
v_A = cand(lam49, P_A, s_o, a_o, h_o)
v_B = cand(lam49, P_B, s_o, a_o, h_o)
print("    sp_mu (X, zX, ..., z^4 X)          %s" % ("ZERO" if v_sp == 0 else "nonzero"))
print("    o_lam (X, zX, ..., z^4 X, 1)       %s     <- AK25 report ZERO" % ("ZERO" if v_A == 0 else "nonzero"))
print("    o_lam (X, zX, ..., z^4 X, 1, -1)   %s     <- the other component" % ("ZERO" if v_B == 0 else "nonzero"))
print()
if v_sp != 0 and v_A == 0:
    print("    >>> Example 4.9 REPRODUCED.  The model of their Section 4.3 is right.")
    print("    >>> and the other component is %s" % ("ALSO zero -- component is NOT the answer"
                                                     if v_B == 0 else "NONZERO -- component IS the answer"))
else:
    print("    >>> Example 4.9 NOT reproduced.  The model of their Section 4.3 is wrong,")
    print("        and no question should be sent on the strength of it.")
