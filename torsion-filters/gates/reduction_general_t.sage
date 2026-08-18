# -*- coding: utf-8 -*-
# Towards the converse: the general-t analogue of (eq:Cmu), and a defect of mine corrected.
#
# *** THE DEFECT, ON THE RECORD ***
# kernel_residual2.sage reported a "KERNEL dim" as (number of nu) minus (rank of the value matrix),
# with the value matrix built from TWELVE evaluation points.  The rank of a matrix with 12 columns
# is at most 12, so every "kernel dimension" I printed was capped by my own sampling and not by the
# mathematics: 68 of 80, 239 of 251 and the rest are artefacts.  The J2 result -- cells match
# distinct values, 0 split, 0 merged -- is unaffected, since that is a statement about
# proportionality of value vectors and 12 random points is a strong test of it.  The rank is not.
# Here the number of points is taken well above the number of cells, so the rank is real.
#
# WHY THE RANK IS THE POINT.  Section 8 of arXiv:2608.09619 reduces the t = 2 problem to
#     Psi_R(lambda) = sum_{l(mu)<=R} C_mu(lambda) o_mu(A),   Psi_R = 0 <=> C_mu = 0 for all mu,
# and the biconditional is legitimate for exactly one reason: {o_mu(A) : l(mu)<=R} = the irreducible
# Sp(2R) characters, which are LINEARLY INDEPENDENT.  Freezing t/2-1 of the pairs at roots of unity
# destroys part of that independence -- that is the whole content of the interior.  So the general-t
# analogue of (eq:Cmu) needs to know exactly what independence survives.
#
# THE CLAIM BEING TESTED.  Brick 1 says sp_mu(z,frozen) = 0 iff the folded classes miss one; brick 2
# says J2(mu) determines the survivor up to a scalar.  If, on top of that, the distinct J2-cells give
# LINEARLY INDEPENDENT functions, then
#
#     Phi_t(lambda) = 0   <=>   for every J2-cell K :   sum_{mu in K} s_mu * C_mu(lambda) = 0,
#
# which is the general-t (eq:Cmu): the vanishing is again a finite system of independent linear
# conditions, one per cell instead of one per label.  That is the platform the converse needs.
#
# CONTROLS able to fail:
#   K1  acceptance, fatal: sp_empty = 1, sp_(1) = sum(x + 1/x).
#   K2  points > cells, asserted, so the rank cannot be capped by the sampling again.
#   K3  rank == number of distinct nonzero cells is the claim; a shortfall is printed as a
#       DEPENDENCY DEFICIT and means the cells are NOT independent.
#   K4  the rank is recomputed on a second, independent set of points and must agree.
#
# Authors: Carles Marin, Claude (AI assistant).

TS = [4, 6, 8]
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


def _accept():
    xs = [F(7), F(11)]
    d = shifted([0, 0], 2)
    ok = (alt(d, xs) / alt(d, xs) == 1 and
          alt(shifted([1, 0], 2), xs) / alt(d, xs) == xs[0] + 1 / xs[0] + xs[1] + 1 / xs[1])
    print("K1 acceptance: %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        raise SystemExit(1)


_accept()


def cls(m, t):
    return min(m % t, (-m) % t)


def covers(m, t):
    h = set(cls(x, t) for x in m)
    return all(k in h for k in range(1, t // 2))


def J2(m, t):
    cs = [cls(x, t) for x in m]
    keep = [x for x in m if cls(x, t) in (0, t // 2) or cs.count(cls(x, t)) > 1]
    return (tuple(sorted(cs)), tuple(sorted(keep, reverse=True)))


def build(t, R, MAX, K, seed):
    r = R - (t // 2 - 1)
    zt = zeta(t)
    frozen = [zt ** k for k in range(1, t // 2)]
    set_random_seed(seed)
    PTS = []
    while len(PTS) < K:
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        xs = list(zz) + frozen
        if len(set([x for x in xs] + [1 / x for x in xs])) == 2 * R:
            PTS.append(xs)
    d = shifted([0] * R, R)
    rows = []
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=R):
            nu = list(l) + [0] * (R - len(l))
            m = shifted(nu, R)
            if not covers(m, t):
                continue
            vec = []
            ok = True
            for xs in PTS:
                den = alt(d, xs)
                if den == 0:
                    ok = False
                    break
                vec.append(alt(m, xs) / den)
            if not ok or all(v == 0 for v in vec):
                continue
            rows.append((nu, m, tuple(vec)))
    return rows, r


print("")
print("=" * 98)
print("IS THE VANISHING STILL A SYSTEM OF INDEPENDENT CONDITIONS AFTER FREEZING?")
print("=" * 98)
print("")
print("   t   R free |nu|<=  points  nonzero mu  J2 cells   RANK   rank2   deficit")
print("  " + "-" * 88)

for t in TS:
    fz = t // 2 - 1
    for R in (fz + 1, fz + 2):
        r = R - fz
        if r < 1:
            continue
        MAX = 18 if R <= 3 else 12
        K = 260
        rows, r = build(t, R, MAX, K, 8100 + 10 * t + R)
        cells = {}
        for nu, m, v in rows:
            cells.setdefault(J2(m, t), []).append((nu, m, v))
        assert K > len(cells), "K2 GUARD: points (%d) must exceed cells (%d)" % (K, len(cells))
        M = matrix(F, [list(v) for _, _, v in rows])
        rk = M.rank()
        rows2, _ = build(t, R, MAX, 120 + len(cells), 9900 + 10 * t + R)
        M2 = matrix(F, [list(v) for _, _, v in rows2])
        rk2 = M2.rank()
        print("  %2d %3d %4d %6d %7d %11d %9d %6d %7d %9d%s"
              % (t, R, r, MAX, K, len(rows), len(cells), rk, rk2, len(cells) - rk,
                 "   <-- DEPENDENCY DEFICIT" if len(cells) != rk else ""))

print("")
print("  K3: 'deficit' = cells minus rank.  Zero means the distinct J2-cells are linearly")
print("      INDEPENDENT, and then Phi_t = 0 is equivalent to one linear condition per cell --")
print("      the general-t form of (eq:Cmu).  A positive deficit means the freezing collapses")
print("      even further than J2 records, and the converse needs those extra relations too.")
print("  K4: 'rank2' is the same rank on an independent set of points; it must agree.")
print("")
print("DONE")
