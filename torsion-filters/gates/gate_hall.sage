# -*- coding: utf-8 -*-
# GATE on the Hall / matching framing.  Is it a matching condition, or just a rank drop in a
# fancy hat?
#
# WHAT I CLAIMED.  "The criterion is a MATCHING condition, not a count: the frozen block needs a
# system of distinct representatives of the non-fixed classes, and Hall's condition reduces to
# COVER."  Before gating that against the literature it has to be gated against itself, because
# there is a cheaper explanation available:
#
#   the determinant vanishes identically in the free variables  <=>  the FROZEN rows are linearly
#   dependent  <=>  rank(frozen block) < t/2 - 1,
#
# which is elementary (generic rows raise the rank to the maximum available).  Then everything
# hangs on ONE identity:
#
#   H  rank of the (t/2-1) x R frozen block  =  number of DISTINCT non-fixed folded classes present.
#
# If H holds, the "matching" language is a reformulation and adds nothing: each column carries one
# class, columns of a fixed class are ZERO, columns of equal class are PROPORTIONAL, and the
# class-vectors v(c) = (zeta^{kc} - zeta^{-kc})_{k=1..t/2-1} are linearly independent.  That last
# is the only real ingredient, and it is the invertibility of the discrete SINE matrix -- classical.
#
# So this script tests three things, in order of how much they would cost me:
#   S1  the sine matrix  S[k][c] = zeta^{kc} - zeta^{-kc},  k,c = 1..t/2-1,  is invertible.
#   S2  rank(frozen block) = #distinct non-fixed classes present, on every nu in range.
#   S3  hence COVER <=> full rank <=> nonvanishing, and the three agree with the determinant.
#
# CONTROLS able to fail:
#   K1  a decoy rank formula -- "number of distinct classes INCLUDING the fixed ones" -- must be
#       wrong, or the fixed classes are not doing what I say.
#   K2  S1 is checked for every even t up to 40, not just the ones used elsewhere.
#   K3  S3 is checked against the determinant at 3 independent points.
#
# Authors: Carles Marin, Claude (AI assistant).

TS = [4, 6, 8, 10, 12]
BIGTS = [t for t in range(4, 41, 2)]
L = lcm(BIGTS)
p = next_prime(10 ** 12)
while (p - 1) % L != 0:
    p = next_prime(p)
F = GF(p)
G = F.multiplicative_generator()
print("field GF(%d);  p-1 divisible by lcm(4,6,..,40) = %d" % (p, L))
for t in BIGTS:
    assert (p - 1) % t == 0, "GUARD: %d does not divide p-1" % t
print("guard: every even t in [4,40] divides p-1  ->  PASS")


def zeta(t):
    z = G ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t))
    return z


def shifted(nu, R):
    return [nu[j] + R - j for j in range(R)]


def alt(ms, xs):
    n = len(ms)
    return matrix(F, n, n, lambda a, b: xs[a] ** ms[b] - xs[a] ** (-ms[b])).det()


def cls(m, t):
    return min(m % t, (-m) % t)


# ---------------------------------------------------------------- S1: the sine matrix is invertible
print("")
print("=" * 92)
print("S1  the sine matrix  S[k][c] = zeta^{kc} - zeta^{-kc},  k,c = 1..t/2-1")
print("=" * 92)
bad = []
for t in BIGTS:
    n = t // 2 - 1
    if n < 1:
        continue
    zt = zeta(t)
    S = matrix(F, n, n, lambda k, c: zt ** ((k + 1) * (c + 1)) - zt ** (-(k + 1) * (c + 1)))
    if S.det() == 0:
        bad.append(t)
print("  even t from 4 to 40 : %d matrices, singular ones: %s"
      % (len(BIGTS), bad if bad else "NONE  ->  invertible, as the classical DST is"))
if bad:
    raise SystemExit(1)

# ------------------------------------------------- S2/S3: the rank identity and the whole criterion
print("")
print("=" * 92)
print("S2/S3  rank(frozen block) = #distinct non-fixed classes, and the criterion that follows")
print("=" * 92)
print("")
print("   t   R free |nu|<=  tested   S2 wrong   decoy wrong   S3 wrong (vs determinant)")
print("  " + "-" * 84)

for t in TS:
    fz = t // 2 - 1
    for r in (1, 2, 3):
        R = fz + r
        MAX = 14 if R <= 4 else 10
        zt = zeta(t)
        frozen = [zt ** k for k in range(1, t // 2)]
        rho0 = shifted([0] * R, R)
        set_random_seed(4400 + 10 * t + R)
        PTS = []
        tries = 0
        while len(PTS) < 3 and tries < 600:
            tries += 1
            zz = [F.random_element() for _ in range(r)]
            if any(x == 0 for x in zz):
                continue
            xs = list(zz) + frozen
            if len(set([x for x in xs] + [1 / x for x in xs])) == 2 * R:
                PTS.append(xs)
        if len(PTS) < 3:
            print("  %2d %3d %4d : no admissible points -- SKIPPED" % (t, R, r))
            continue
        n = w2 = wd = w3 = 0
        for size in range(MAX + 1):
            for l in Partitions(size, max_length=R):
                nu = list(l) + [0] * (R - len(l))
                m = shifted(nu, R)
                n += 1
                # the frozen block: rows k = 1..t/2-1, columns j = 1..R
                B = matrix(F, fz, R,
                           lambda a, b: frozen[a] ** m[b] - frozen[a] ** (-m[b]))
                rk = B.rank()
                nonfixed = set(cls(x, t) for x in m) - set([0, t // 2])
                allcls = set(cls(x, t) for x in m)
                if rk != len(nonfixed):
                    w2 += 1
                if rk == len(allcls):
                    wd += 1
                vals = []
                ok = True
                for xs in PTS:
                    den = alt(rho0, xs)
                    if den == 0:
                        ok = False
                        break
                    vals.append(alt(m, xs) / den)
                if not ok:
                    continue
                v = all(x == 0 for x in vals)
                if (rk == fz) == v:              # full rank should mean NONvanishing
                    w3 += 1
        print("  %2d %3d %4d %6d %7d %10d %13d %14d"
              % (t, R, r, MAX, n, w2, wd, w3))

print("")
print("  S2 wrong must be 0: then 'matching' is only bookkeeping for a rank, and the single real")
print("  ingredient is S1, the invertibility of the sine matrix -- which is classical.")
print("  'decoy wrong' counts the shapes where the rank ALSO equals the count including the fixed")
print("  classes; it must be far from the total, or the fixed classes are not being killed.")
print("  S3 wrong must be 0: full rank of the frozen block <=> the character does not vanish.")
print("")
print("DONE")
