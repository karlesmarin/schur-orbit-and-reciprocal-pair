# -*- coding: utf-8 -*-
# NPP on the det = -1 component: the count runs over the LONG roots.
#
# WHERE THIS COMES FROM.  Nadimpalli-Pattanayak-Prasad, arXiv:2504.14684, Corollary 3.7, for a
# CONNECTED reductive G and the torsion element rho^(e^{2 pi i/m}):
#
#     chi_lambda != 0   <=>   #{a>0 : m | <rho^,a>}  =  #{a>0 : m | <lambda^+rho^,a>}.
#
# In type A the roots are e_i - e_j, so the count is #{i<j : m | beta_i - beta_j} and the criterion
# is Littlewood's: no residue class of beta is over-full.  That is branch (a).  It cannot reach our
# branch, which is a condition on beta_i + beta_j -- a root of type e_i + e_j, absent from A.
#
# The paper already folds the problem for us.  Lemma AtoSp: on the det = -1 coset the characters are
# SYMPLECTIC in the free variables, "the orbit Lie algebra of D_{r+1} being C_r rather than the fixed
# subalgebra B_r".  So the analogue of Corollary 3.7 for the non-identity component must be a count
# in C_R -- and C_R, unlike A, has TWO root lengths.  This script tests the first brick:
#
#     PROPOSAL (the long-root half).  Freezing a variable at a primitive t-th root of unity kills
#     sp_nu identically when a whole column of the symplectic bialternant dies, i.e. when
#         exists k in [1, t/2-1] :  t | 2k*m_j  for every j,      m_j = nu_j + R - j + 1,
#     and 2m_j = <nu + rho_C, 2e_j> is exactly the pairing with the LONG root 2e_j.
#
# So Littlewood/NPP is the short-root count e_a - e_b, and the interior's second branch is the
# long-root count 2e_a.  That is what "non-simply-laced" buys, and it is the shape KLP09 predicts.
#
# WHAT IS MEASURED, and the controls:
#   K1  LONG => sp_nu vanishes.  Forced by the column argument; a failure means the code is wrong.
#   K2  sp_nu vanishes => LONG.  NOT forced -- a determinant can die by column DEPENDENCE without
#       any column being zero.  This is the real content and it is allowed to fail; the shapes
#       where it does are printed, because they are the residue the long-root count does not see.
#   K3  a decoy that must do worse: "some m_j is even".
#   K4  the free part must be genuinely free: the count of nu tested and the number of random
#       points are printed, and vanishing is decided at 3 independent points.
#
# Authors: Carles Marin, Claude (AI assistant).

TS = [4, 6, 8, 10]
L = lcm(TS)
p = next_prime(10 ** 9)
while (p - 1) % L != 0:
    p = next_prime(p)
F = GF(p)
G = F.multiplicative_generator()
for t in TS:
    assert (p - 1) % t == 0, "GUARD: %d does not divide p-1" % t
print("field GF(%d); every t in %s divides p-1  ->  guard PASS" % (p, TS))


def zeta(t):
    z = G ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t))
    return z


def shifted(nu, R):
    """m_j = nu_j + R - j + 1, j = 1..R  (0-indexed: nu[j] + R - j).  m = (R,..,1) for nu = 0."""
    return [nu[j] + R - j for j in range(R)]


def sp_value(nu, xs, R):
    """sp_nu(x_1..x_R) = det(x^m - x^-m) / det(x^d - x^-d), d = shifted(0);
    ONE definition of the shift, used by both the numerator and the predicate."""
    m = shifted(nu, R)
    d = shifted([0] * R, R)
    num = matrix(F, R, R, lambda a, b: xs[a] ** m[b] - xs[a] ** (-m[b])).det()
    den = matrix(F, R, R, lambda a, b: xs[a] ** d[b] - xs[a] ** (-d[b])).det()
    return num, den


# ACCEPTANCE, fatal: sp_nu of the empty partition must be 1, and sp_(1) on R=2 must be x+1/x+y+1/y.
def _accept():
    R = 2
    xs = [F(7), F(11)]
    n0, d0 = sp_value([0, 0], xs, R)
    n1, d1 = sp_value([1, 0], xs, R)
    want = xs[0] + 1 / xs[0] + xs[1] + 1 / xs[1]
    ok = (n0 / d0 == 1) and (n1 / d1 == want)
    print("acceptance: sp_empty=1 and sp_(1)=sum(x+1/x)  ->  %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        raise SystemExit(1)


_accept()


def LONG(nu, R, t):
    m = shifted(nu, R)
    for k in range(1, t // 2):
        if all((2 * k * mj) % t == 0 for mj in m):
            return True
    return False


def DECOY(nu, R, t):
    return any(mj % 2 == 0 for mj in shifted(nu, R))


print("")
print("=" * 94)
print("THE LONG-ROOT COUNT AGAINST THE MEASURED VANISHING OF sp_nu AT A FROZEN t-th ROOT")
print("=" * 94)
print("")
print("   t   R  free  |nu|<=   nu tested   sp=0   LONG   K1 fail   K2 fail   decoy agrees")
print("  " + "-" * 88)

RESID = []
for t in TS:
    frozen_n = t // 2 - 1
    for R in (frozen_n + 1, frozen_n + 2):
        r = R - frozen_n
        if r < 1:
            continue
        MAX = 16 if R <= 3 else 12
        zt = zeta(t)
        frozen = [zt ** k for k in range(1, t // 2)]
        set_random_seed(4000 + 10 * t + R)
        PTS = []
        while len(PTS) < 3:
            zz = [F.random_element() for _ in range(r)]
            if any(x == 0 for x in zz):
                continue
            xs = list(zz) + frozen
            if len(set([x for x in xs] + [1 / x for x in xs])) == 2 * R:
                PTS.append(xs)
        ntest = nz = nlong = k1 = k2 = ndec = 0
        for size in range(MAX + 1):
            for l in Partitions(size, max_length=R):
                nu = list(l) + [0] * (R - len(l))
                m = shifted(nu, R)
                ntest += 1
                vals = []
                ok = True
                for xs in PTS:
                    num, den = sp_value(nu, xs, R)
                    if den == 0:
                        ok = False
                        break
                    vals.append(num / den)
                if not ok:
                    ntest -= 1
                    continue
                v = all(x == 0 for x in vals)
                lg = LONG(nu, R, t)
                if v:
                    nz += 1
                if lg:
                    nlong += 1
                if lg and not v:
                    k1 += 1
                    if k1 <= 2:
                        print("     K1 FAIL t=%d R=%d nu=%s m=%s" % (t, R, nu, m))
                if v and not lg:
                    k2 += 1
                    if len(RESID) < 14:
                        RESID.append((t, R, nu, m))
                if DECOY(nu, R, t) == v:
                    ndec += 1
        print("  %2d %3d %5d %7d %11d %6d %6d %9d %9d %10d/%d"
              % (t, R, r, MAX, ntest, nz, nlong, k1, k2, ndec, ntest))

print("")
print("  K1 must be 0: a zero column kills the determinant, so LONG forces vanishing.")
print("  K2 is the residue -- vanishing WITHOUT a dead column, i.e. genuine column dependence.")
print("  It is the part of the second branch the long-root count alone does not explain.")
print("")
if RESID:
    print("  the residue, by hand (nu with sp_nu = 0 and no dead column):")
    for t, R, nu, m in RESID:
        print("    t=%d R=%d  nu=%-18s shifted m=%s" % (t, R, nu, m))
else:
    print("  the residue is EMPTY: at these sizes the long-root count is the whole criterion.")
print("")
print("DONE")
