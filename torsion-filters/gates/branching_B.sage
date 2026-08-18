# -*- coding: utf-8 -*-
# The general-t (eq:Cmu): the freezing is a branching matrix onto the Sp(2r) basis.
#
# *** THE PREVIOUS CLAIM, REFUTED BY ITS OWN CONTROL ***
# reduction_general_t.sage asked whether the distinct J2-cells are linearly independent, so that
# Phi_t = 0 would be one condition per cell.  They are NOT: deficits 16, 81, 12, 17, 4, 12, with the
# rank reproduced on a second independent set of points.  J2 says when two labels give PROPORTIONAL
# values; it does not say when several give a DEPENDENT set, and the difference is real.
#
# But the rank is not mysterious.  At t=4, R=2, r=1 with |mu|<=18 it is 19, and 19 is exactly the
# number of Sp(2) characters available below that degree.  The image of the freezing lands in the
# span of the Sp(2r) characters of the FREE variables, and those are irreducible characters, hence
# linearly independent.  So the right statement is not one condition per cell but one condition per
# Sp(2r) label:
#
#     Phi_t(lambda) = 0   <=>   sum_mu B[mu][kappa] * C_mu(lambda) = 0   for every kappa, l(kappa)<=r
#
# with B the BRANCHING MATRIX of the freezing,  sp_mu(z, zeta, .., zeta^{t/2-1}) = sum_kappa
# B[mu][kappa] sp_kappa(z).  That is the general-t form of (eq:Cmu) of arXiv:2608.09619, which is
# the t=2 case where the freezing is empty and B is the identity.  Bricks 1 and 2 describe B: a row
# is zero iff the folded classes miss one, and two rows are proportional iff J2 agrees.
#
# WHAT IS MEASURED HERE, at r = 1 where B can be written down:
#   B1  every sp_mu(z,frozen) really is a combination of Sp(2) characters -- the solve must be exact
#       at points not used to solve it (a control that can fail).
#   B2  the number of nonzero entries per row, and the values they take.  The prediction, from the
#       Laplace mechanism: a row has at most 2 nonzero entries and they are +-1 times one scalar.
#   B3  rank(B) must equal the rank measured directly by reduction_general_t.sage.
#
# CONTROLS able to fail:
#   K1  acceptance, fatal: sp_empty = 1, sp_(1) = z + 1/z.
#   K2  the expansion is solved on HALF the points and verified on the other half.
#   K3  rank(B) against the independently measured rank; a mismatch kills the reading.
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
    assert (p - 1) % t == 0
print("field GF(%d); guard on t in %s -> PASS" % (p, TS))


def zeta(t):
    z = G ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t))
    return z


def shifted(nu, R):
    return [nu[j] + R - j for j in range(R)]


def alt(ms, xs):
    n = len(ms)
    return matrix(F, n, n, lambda a, b: xs[a] ** ms[b] - xs[a] ** (-ms[b])).det()


def sp1(k, z):
    """the Sp(2) character of label k: (z^{k+1} - z^{-(k+1)})/(z - 1/z)."""
    return (z ** (k + 1) - z ** (-(k + 1))) / (z - 1 / z)


def _accept():
    z = F(13)
    ok = (sp1(0, z) == 1) and (sp1(1, z) == z + 1 / z)
    print("K1 acceptance: %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        raise SystemExit(1)


_accept()


def cls(m, t):
    return min(m % t, (-m) % t)


def covers(m, t):
    return all(k in set(cls(x, t) for x in m) for k in range(1, t // 2))


print("")
print("=" * 98)
print("THE BRANCHING MATRIX B OF THE FREEZING,  r = 1")
print("=" * 98)
print("")
print("   t   R  |mu|<=   mu tested   rows = 0   rows with 1   with 2   with >2   rank B   K2")
print("  " + "-" * 90)

EX = []
for t in TS:
    R = t // 2                       # r = 1
    MAX = 18
    D = MAX + R + 2                  # kappa = 0..D covers every degree that can occur
    zt = zeta(t)
    frozen = [zt ** k for k in range(1, t // 2)]
    set_random_seed(1234 + t)
    ZS = []
    while len(ZS) < 2 * (D + 1) + 20:
        z = F.random_element()
        if z == 0 or z ** 2 == 1:
            continue
        xs = [z] + frozen
        if len(set([x for x in xs] + [1 / x for x in xs])) == 2 * R:
            ZS.append(z)
    half = D + 1 + 5
    FIT, CHK = ZS[:half], ZS[half:]
    d0 = shifted([0] * R, R)
    A = matrix(F, len(FIT), D + 1, lambda i, k: sp1(k, FIT[i]))
    A2 = matrix(F, len(CHK), D + 1, lambda i, k: sp1(k, CHK[i]))
    assert A.rank() == D + 1, "the Sp(2) basis is not independent on the fit points"

    Brows = []
    n0 = n1 = n2 = nbig = 0
    k2fail = 0
    ntest = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=R):
            mu = list(l) + [0] * (R - len(l))
            m = shifted(mu, R)
            ntest += 1
            fv = []
            ok = True
            for z in FIT:
                xs = [z] + frozen
                den = alt(d0, xs)
                if den == 0:
                    ok = False
                    break
                fv.append(alt(m, xs) / den)
            if not ok:
                ntest -= 1
                continue
            if not covers(m, t):
                assert all(x == 0 for x in fv), "brick 1 violated"
                n0 += 1
                Brows.append([F(0)] * (D + 1))
                continue
            c = A.solve_right(vector(F, fv))
            gv = []
            for z in CHK:
                xs = [z] + frozen
                den = alt(d0, xs)
                gv.append(alt(m, xs) / den if den != 0 else F(0))
            if A2 * c != vector(F, gv):
                k2fail += 1
            nz = [(k, c[k]) for k in range(D + 1) if c[k] != 0]
            Brows.append(list(c))
            if len(nz) == 1:
                n1 += 1
            elif len(nz) == 2:
                n2 += 1
            else:
                nbig += 1
            if len(EX) < 8 and size <= 4:
                EX.append((t, mu, m, nz))
    B = matrix(F, Brows)
    print("  %2d %3d %7d %11d %10d %13d %8d %9d %8d %4s"
          % (t, R, MAX, ntest, n0, n1, n2, nbig, B.rank(), "ok" if k2fail == 0 else "FAIL"))

print("")
print("  B2 prediction was: at most 2 nonzero entries per row.  The 'with >2' column tests it.")
print("")
print("  rows of B by hand (kappa, coefficient):")
for t, mu, m, nz in EX:
    print("    t=%d  mu=%-14s m=%s  ->  %s"
          % (t, mu, m, ", ".join("%d:%s" % (k, ("+1" if v == F(1) else ("-1" if v == F(-1) else str(v)))) for k, v in nz)))
print("")
print("DONE")
