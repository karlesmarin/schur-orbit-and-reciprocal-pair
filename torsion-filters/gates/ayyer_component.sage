# ayyer_component.sage
#
# THE QUESTION.  Ayyer-Kumari, arXiv:2501.00275, Section 4.3 compares Sp(2tn) with O(2tn+2) and
# finds the pair misbehaves: with lambda = mu_tn + (mu, 0),
#
#   "Suppose mu_tn != 0 and sp_mu(X, zX, ..., z^{t-1}X) is nonzero.  Then that does not imply
#    that o^even_lambda(X, zX, ..., z^{t-1}X, 1) is nonzero."
#
# closing with "... might be false".  Their alphabet adjoins the letter 1: that is a torus element
# of the IDENTITY component of O(2tn+2), tn+1 parameters.
#
# Our Lemma 8.8 is an identity in Lambda, for EVERY alphabet W:
#
#         o_nu(W, 1, -1) = sp_nu(W).
#
# W together with the two fixed points 1 and -1 is a torus element of the NON-IDENTITY component
# of O(2tn+2) -- tn parameters plus the two fixed points -- and there the even orthogonal
# character IS the symplectic character, no sign, no length hypothesis, because the orbit Lie
# algebra of D_{tn+1} is C_tn rather than the fixed subalgebra B_tn.
#
# So the test is: with W the Sp(2tn) torus of their own X, zX, ..., z^{t-1}X, compare
#
#   (A)  sp_mu(W)  against  o_lambda(W, 1)        their identity-component comparison
#   (B)  sp_mu(W)  against  o_lambda(W, 1, -1)    the other component
#
# and ask whether the vanishing loci agree in (B) where they do not in (A).
#
# CONVENTION, taken from our own verified anc/sec8_derivation.sage and not guessed:
#   sp_nu on the Sp(2r) torus       p_k -> sum_j (z_j^k + z_j^-k)          -- 2r eigenvalues
#   o_nu at A = {1,-1,z_j^{+-1}}    p_k -> 1 + (-1)^k + sum_j (z_j^k + z_j^-k)
#   the orbit torus                 p_k -> t*[t | k] * sum_i (x_i^k + x_i^-k)
#   adjoining the O-parameter y=1   p_k -> ... + 2
#
# CONTROLS, each of which must pass before the next means anything:
#   C0  Sage's o and sp must satisfy our Lemma 8.8 under this convention.  An earlier run of this
#       file evaluated sp on r values instead of 2r and C0 failed 52 of 60; the numbers it printed
#       were meaningless.  That is what C0 is for.
#   C1  the phenomenon must occur in range: there must BE pairs with sp_mu != 0 and o_lambda = 0
#       on the identity component, or the test is vacuous.
#
# Authors: Carles Marin, Claude (AI assistant).

import sys
sys.stdout.reconfigure(line_buffering=True)

Sym = SymmetricFunctions(QQ)
p, O, SP = Sym.p(), Sym.o(), Sym.sp()


def ev(f, P):
    """evaluate a universal character at an alphabet given by its power sums P(k)"""
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


print("=" * 94)
print("C0  CONTROL -- Sage's o and sp against our Lemma 8.8, in OUR convention")
print("        o_nu(1, -1, z_1, 1/z_1, ...) = sp_nu(z_1, 1/z_1, ...)")
print("=" * 94)
bad = tot = 0
for r in (1, 2):
    zs = [QQ(3) / 2, QQ(5) / 3][:r]
    P_sp = lambda k: sum(z ** k + z ** (-k) for z in zs)
    P_o = lambda k: 1 + (-1) ** k + sum(z ** k + z ** (-k) for z in zs)
    for size in range(0, 7):
        for nu in Partitions(size):
            tot += 1
            if ev(univ(O, nu), P_o) != ev(univ(SP, nu), P_sp):
                bad += 1
print("    %d partitions at r=1,2 ; %d disagreements" % (tot, bad))
if bad:
    print("    *** dictionary still wrong -- the rest of this file is meaningless ***")
    raise SystemExit(1)
print("    dictionary confirmed.  Proceeding.")

print()
print("=" * 94)
print("C0b CALIBRATION -- Koike-Terada determinants against Sage's bases")
print("        the character at a torus element is a determinant of size l(lambda) in the h's,")
print("        so |lambda| never enters the cost.  This is what makes Example 4.9 reachable.")
print("=" * 94)

def h_from_P(P, M):
    h = [1]
    for m in range(1, M + 1):
        h.append(sum(P(i) * h[m - i] for i in range(1, m + 1)) / m)
    return h

def kt(lam, P, kind):
    """sp_lambda or o_lambda at an alphabet given by its power sums, by the Koike-Terada
    determinant.  Size l(lambda); h computed once by Newton."""
    L = [k for k in lam if k > 0]
    l = len(L)
    if l == 0:
        return 1
    M = max(L) + l + 2
    h = h_from_P(P, M)
    H = lambda d: 0 if d < 0 or d > M else h[d]
    if kind == 'sp':
        rows = [[H(L[i] - (i + 1) + (j + 1)) - H(L[i] - (i + 1) - (j + 1)) for j in range(l)]
                for i in range(l)]
    else:
        rows = [[H(L[i] - (i + 1) + (j + 1)) - H(L[i] - (i + 1) - (j + 1) + 2) for j in range(l)]
                for i in range(l)]
    return matrix(QQ, rows).det()

calib_bad = calib_tot = 0
for r in (1, 2, 3):
    zs = [QQ(3) / 2, QQ(5) / 3, QQ(7) / 4][:r]
    P_sp = lambda k: sum(z ** k + z ** (-k) for z in zs)
    for size in range(0, 8):
        for nu in Partitions(size, max_length=r):
            calib_tot += 1
            if kt(list(nu), P_sp, 'sp') != ev(univ(SP, nu), P_sp):
                calib_bad += 1
            if kt(list(nu), P_sp, 'o') != ev(univ(O, nu), P_sp):
                calib_bad += 1
print("    %d shapes, %d determinant/basis disagreements" % (calib_tot, calib_bad))
if calib_bad:
    print("    *** Koike-Terada form wrong -- Example 4.9 below is not to be trusted ***")
else:
    print("    determinants agree with Sage on both series.")

print()
print("=" * 94)
print("EXAMPLE 4.9 of arXiv:2501.00275, computed  (t=5, n=2, |mu|=44, |lambda|=55)")
print("=" * 94)
mu49 = [9, 8, 6, 5, 4, 4, 3, 2, 2, 1]
lam49 = [10, 9, 7, 6, 5, 5, 4, 3, 3, 2, 1]
t49, n49 = 5, 2
Xs49 = [QQ(3) / 2, QQ(5) / 3]
P_W49 = lambda k: (t49 * sum(x ** k + x ** (-k) for x in Xs49)) if k % t49 == 0 else 0
P_A49 = lambda k: P_W49(k) + 2
P_B49 = lambda k: P_W49(k) + 1 + (-1) ** k
v_sp = kt(mu49, P_W49, 'sp')
v_A = kt(lam49, P_A49, 'o')
v_B = kt(lam49, P_B49, 'o')
print("    sp_mu(X, zX, ..., z^4 X)          = %s" % ("0" if v_sp == 0 else "nonzero"))
print("    o_lambda(X, zX, ..., z^4 X, 1)    = %s   <- they report 0" % ("0" if v_A == 0 else "nonzero"))
print("    o_lambda(X, zX, ..., z^4 X, 1,-1) = %s   <- the other component" % ("0" if v_B == 0 else "nonzero"))
print()
print("    their Example 4.9 is reproduced iff the first is nonzero and the second is 0.")

print()
print("=" * 94)
print("THE TEST -- lambda = mu_tn + (mu, 0), Section 4.3's construction")
print("=" * 94)
print("%-4s %-4s %-8s %-9s %-13s %-14s %-14s"
      % ("t", "n", "mu's", "sp != 0", "A mismatch", "B mismatch", "B fixes A"))
print("-" * 94)

for t in (2, 3, 4, 5):
    for n in (1, 2):
        tn = t * n
        if tn > 6:
            continue
        K = CyclotomicField(t) if t > 2 else QQ
        Xs = [QQ(3) / 2, QQ(5) / 3][:n]
        # the Sp(2tn) torus of X, zX, ..., z^{t-1}X : p_k = t*[t|k]*sum_i (x_i^k + x_i^-k)
        P_W = lambda k: (t * sum(x ** k + x ** (-k) for x in Xs)) if k % t == 0 else 0
        P_A = lambda k: P_W(k) + 2                      # adjoin the O-parameter y = 1
        P_B = lambda k: P_W(k) + 1 + (-1) ** k          # adjoin the two fixed points 1 and -1
        nmu = spnz = mA = mB = fixed = 0
        ex = []
        for size in range(0, 11):
            for mu in Partitions(size, max_length=tn):
                M = list(mu) + [0] * (tn - len(mu))
                last = M[-1]
                lam = [M[i] + last for i in range(tn)] + [last]
                nmu += 1
                v_sp = kt(M, P_W, 'sp')
                v_A = kt(lam, P_A, 'o')
                v_B = kt(lam, P_B, 'o')
                if v_sp != 0:
                    spnz += 1
                    if v_A == 0:
                        mA += 1
                        if v_B != 0:
                            fixed += 1
                        if len(ex) < 3:
                            ex.append((M, lam, v_B == 0))
                    if v_B == 0:
                        mB += 1
        print("%-4d %-4d %-8d %-9d %-13d %-14d %-14d" % (t, n, nmu, spnz, mA, mB, fixed))
        for M, lam, bzero in ex:
            print("      mu=%-16s lambda=%-18s  B also zero? %s" % (M, lam, bzero))

print()
print("=" * 94)
print("READING")
print("=" * 94)
print("  C1 : 'A mismatch' must be > 0 somewhere, or the phenomenon is out of range.")
print("  YES: 'B mismatch' is 0 wherever 'A mismatch' is not, i.e. 'B fixes A' = 'A mismatch'.")
print("  NO : 'B mismatch' tracks 'A mismatch' -- the component is not what is going on.")
