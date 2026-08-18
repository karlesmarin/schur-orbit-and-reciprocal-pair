# -*- coding: utf-8 -*-
# The residual kernel, second attempt: measured against ratios, with no denominator to get wrong.
#
# *** WHAT FAILED IN kernel_residual.sage, ON THE RECORD ***
#
# H2 scored 0 out of 55 at every (t,R).  That was MY bug, not a refutation.  In the clean case the
# surviving Laplace term is [frozen minor] x [free alternant on the leftover columns], and the
# DENOMINATOR splits the same way -- for r = 1 the nu = 0 shifted vector is (t/2, .., 1), whose own
# leftover is [t/2].  So the ratio to compare against is alt(left_nu)/alt([t/2]), and I divided by
# alt(left_nu)/alt(staircase of r variables) instead.  For t=4 that is Z_b/Z_1 where the truth is
# Z_b/Z_2, and the ratio Z_1/Z_2 is not constant.  A denominator, not a mechanism.
#
# H1 was too coarse, and the printed cases say exactly why: J kept exact values only for the FIXED
# classes {0, t/2}.  But when a non-fixed class is hit TWICE, the frozen minor can only use one of
# the two columns and the other goes into the free alternant WITH ITS EXACT VALUE -- so m=[3,1] and
# m=[5,1] share J and must differ, and they do.
#
# THE FIX: DROP THE DENOMINATOR.  The kernel is a statement about ratios between different nu, so
# compare sp_nu against sp_mu directly and never form an absolute normalisation.
#
#   H3  restricted to the CLEAN nu -- every non-fixed class hit exactly once -- sp_nu(z,frozen) is
#       determined, up to a scalar, by the MULTISET OF LEFTOVER COLUMNS  {m_j : cls(m_j) in {0,t/2}}.
#       Both directions are scored: same leftover must give proportional values (that is the kernel),
#       and different leftovers must give non-proportional ones (or the kernel is bigger still).
#
#   H4  the general nu: the invariant must also carry the exact values of the duplicated columns.
#       Candidate: J2(nu) = ( sorted classes , sorted m_j over ALL columns whose class is either
#       fixed or repeated ).  Scored the same way, with H1's old J beside it as the retracted decoy.
#
# CONTROLS able to fail:
#   K1  acceptance, fatal.
#   K2  the kernel dimension, reported: zero would mean no second brick.
#   K3  both directions of H3 and H4, failures printed by hand.
#   K4  the retracted J is scored next to J2 and must do worse.
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
    ok = (alt(shifted([0, 0], 2), xs) / alt(d, xs) == 1 and
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


def clean(m, t):
    cs = [cls(x, t) for x in m]
    return all(cs.count(k) == 1 for k in range(1, t // 2))


def leftover(m, t):
    return tuple(sorted((x for x in m if cls(x, t) in (0, t // 2)), reverse=True))


def J_old(m, t):
    return (tuple(sorted(cls(x, t) for x in m)), leftover(m, t))


def J2(m, t):
    """classes, plus the exact values of every column in a fixed OR repeated class."""
    cs = [cls(x, t) for x in m]
    keep = [x for x in m if cls(x, t) in (0, t // 2) or cs.count(cls(x, t)) > 1]
    return (tuple(sorted(cs)), tuple(sorted(keep, reverse=True)))


def canon(vec):
    for v in vec:
        if v != 0:
            return tuple(x / v for x in vec)
    return None


def score(rows, keyf, t):
    """(cells carrying >1 value, distinct values that appear in >1 cell)"""
    cell = {}
    val = {}
    for nu, m, v in rows:
        c = canon(v)
        cell.setdefault(keyf(m, t), set()).add(c)
        val.setdefault(c, set()).add(keyf(m, t))
    return (sum(1 for s in cell.values() if len(s) > 1),
            sum(1 for s in val.values() if len(s) > 1),
            len(cell))


print("")
print("=" * 98)
print("THE RESIDUAL KERNEL, measured by ratios")
print("=" * 98)
print("")
print("   t   R free |nu|<=  COVER nu  values  KERNEL  | J2: cells  split  merged | J_old split")
print("  " + "-" * 94)

BAD3 = []
BAD4 = []
for t in TS:
    fz = t // 2 - 1
    for R in (fz + 1, fz + 2):
        r = R - fz
        if r < 1:
            continue
        MAX = 18 if R <= 3 else 12
        zt = zeta(t)
        frozen = [zt ** k for k in range(1, t // 2)]
        set_random_seed(7100 + 10 * t + R)
        PTS = []
        while len(PTS) < 12:
            zz = [F.random_element() for _ in range(r)]
            if any(x == 0 for x in zz):
                continue
            xs = list(zz) + frozen
            if len(set([x for x in xs] + [1 / x for x in xs])) == 2 * R:
                PTS.append((zz, xs))
        d = shifted([0] * R, R)
        rows = []
        for size in range(MAX + 1):
            for l in Partitions(size, max_length=R):
                nu = list(l) + [0] * (R - len(l))
                m = shifted(nu, R)
                if not covers(m, t):
                    continue
                vec = []
                bad = False
                for zz, xs in PTS:
                    den = alt(d, xs)
                    if den == 0:
                        bad = True
                        break
                    vec.append(alt(m, xs) / den)
                if bad or all(v == 0 for v in vec):
                    continue
                rows.append((nu, m, tuple(vec)))
        M = matrix(F, [list(v) for _, _, v in rows])
        kern = len(rows) - M.rank()
        nval = len(set(canon(v) for _, _, v in rows))
        s2, m2, c2 = score(rows, J2, t)
        s1, m1, c1 = score(rows, J_old, t)
        print("  %2d %3d %4d %6d %9d %7d %7d  | %9d %6d %7d | %11d"
              % (t, R, r, MAX, len(rows), nval, kern, c2, s2, m2, s1))
        cell = {}
        for nu, m, v in rows:
            cell.setdefault(J2(m, t), set()).add(canon(v))
        for k, s in cell.items():
            if len(s) > 1 and len(BAD4) < 6:
                ex = [(nu, m) for nu, m, v in rows if J2(m, t) == k][:2]
                BAD4.append((t, R, k, ex))

        # ---- H3 on the clean nu only
        cl = [(nu, m, v) for nu, m, v in rows if clean(m, t)]
        g = {}
        for nu, m, v in cl:
            g.setdefault(leftover(m, t), set()).add(canon(v))
        split = sum(1 for s in g.values() if len(s) > 1)
        vals = {}
        for nu, m, v in cl:
            vals.setdefault(canon(v), set()).add(leftover(m, t))
        merged = sum(1 for s in vals.values() if len(s) > 1)
        BAD3.append((t, R, r, len(cl), len(g), split, merged))

print("")
print("  'KERNEL' is the rank deficit: how much of the span the freezing collapses.")
print("  'split' = invariant cells carrying more than one value up to scalar (must be 0).")
print("  'merged' = values appearing in more than one cell (invariant too fine, informative).")
print("")
if BAD4:
    print("  J2 failures, by hand:")
    for t, R, k, ex in BAD4:
        print("    t=%d R=%d J2=%s" % (t, R, k))
        for nu, m in ex:
            print("        nu=%-20s m=%s" % (nu, m))
else:
    print("  J2: no cell carries two values.  The kernel is generated by the differences")
    print("      sp_nu - c*sp_mu over pairs with the same J2 -- classes, plus the exact values")
    print("      of every column in a fixed or repeated class.")

print("")
print("=" * 98)
print("H3  the CLEAN nu: is the leftover multiset the whole invariant?")
print("=" * 98)
print("")
print("   t   R free   clean nu   leftover cells   split (must be 0)   merged")
print("  " + "-" * 74)
for t, R, r, ncl, ng, split, merged in BAD3:
    print("  %2d %3d %4d %10d %16d %19d %8d%s"
          % (t, R, r, ncl, ng, split, merged, "   <-- H3 FAILS" if split else ""))
print("")
print("  If 'split' is 0 the clean part of the kernel is exactly: same leftover columns =>")
print("  proportional.  'merged' > 0 would mean even the leftover is finer than it needs to be.")
print("")
print("DONE")
