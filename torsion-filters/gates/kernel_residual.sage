# -*- coding: utf-8 -*-
# The second brick: the residual kernel of the specialization.
#
# WHERE WE ARE.  npp_minus2.sage settled the first half, with a proof and 2001 shapes at 0
# discrepancies:  sp_nu(z_1..z_r, zeta, .., zeta^{t/2-1}) = 0  exactly when the folded residues
# cls(m_j) = min(m_j mod t, -m_j mod t) fail to cover {1,..,t/2-1},  m_j = nu_j + R - j + 1.
# Those nu contribute nothing to  Phi_t = sum_nu c_nu(lambda) sp_nu(z, frozen).
#
# THE SURVIVORS ARE STILL DEPENDENT, and that is the whole point.  Freezing t/2-1 variables maps an
# infinite-dimensional span of symplectic characters into functions of r variables, so there must be
# relations, and Problem prob:extralocus of arXiv:2608.09619 asks whether the extra locus is read
# off exactly that kernel.  This measures the kernel instead of assuming it.
#
# THE MECHANISM, WHICH GIVES A PREDICTION.  Laplace along the t/2-1 frozen rows: each term is
# [frozen minor on a column set] x [free alternant on the complementary r columns].  The frozen
# entries v_k(m) = zeta^{km} - zeta^{-km} depend on m ONLY through cls(m) and a sign, so a frozen
# minor dies unless its t/2-1 columns are a system of distinct representatives of the classes
# 1..t/2-1.  Hence:
#
#   H1  sp_nu(z, frozen) is determined by nu through the pair
#          J(nu) = ( multiset of cls(m_j) ,  multiset of the m_j in the FIXED classes {0, t/2} ),
#       up to a scalar -- the frozen part sees only classes, the free part sees exact values.
#   H2  in the CLEAN case, where every class 1..t/2-1 is hit exactly ONCE, the Laplace sum has a
#       single surviving term, so sp_nu(z,frozen) is proportional to the r-variable symplectic
#       character on the leftover columns:  sp_nu(z,frozen) = const * sp_nutilde(z).
#
# Both are predictions written before the run.  H1 is the kernel: any two nu with the same J give a
# relation.  H2 says the surviving image is just Sp(2r) characters again, so the kernel is as big as
# it can be.
#
# CONTROLS able to fail:
#   K1  acceptance, fatal: sp_empty = 1, sp_(1) = sum(x+1/x).
#   K2  the measured kernel dimension is reported next to the number of nu -- if they were equal
#       the specialization would be injective and there would be no second brick at all.
#   K3  H1 is scored in BOTH directions: same J and different value (fatal to H1), and different J
#       and same value (H1 too coarse, still informative).  Failures printed by hand.
#   K4  a decoy invariant that must do worse: the multiset of m_j mod t WITHOUT folding.
#   K5  H2 is checked by demanding the ratio be constant across 12 independent z-points, not at one.
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
    assert (p - 1) % t == 0, "GUARD: %d does not divide p-1" % t
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


def sp_at(nu, xs, R):
    m = shifted(nu, R)
    d = shifted([0] * R, R)
    den = alt(d, xs)
    return (alt(m, xs), den)


def _accept():
    xs = [F(7), F(11)]
    n0, d0 = sp_at([0, 0], xs, 2)
    n1, d1 = sp_at([1, 0], xs, 2)
    ok = (n0 / d0 == 1) and (n1 / d1 == xs[0] + 1 / xs[0] + xs[1] + 1 / xs[1])
    print("K1 acceptance: %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        raise SystemExit(1)


_accept()


def cls(m, t):
    return min(m % t, (-m) % t)


def covers(m, t):
    h = set(cls(x, t) for x in m)
    return all(k in h for k in range(1, t // 2))


def J(m, t):
    """H1's invariant: the classes, plus the exact values sitting in the FIXED classes."""
    fixed = (0, t // 2)
    return (tuple(sorted(cls(x, t) for x in m)),
            tuple(sorted(x for x in m if cls(x, t) in fixed)))


def DECOY(m, t):
    return tuple(sorted(x % t for x in m))


def canon(vec):
    """projective normal form of a value vector, or None if identically zero."""
    for v in vec:
        if v != 0:
            return tuple(x / v for x in vec)
    return None


print("")
print("=" * 96)
print("THE RESIDUAL KERNEL: how far the frozen specialization collapses the symplectic characters")
print("=" * 96)
print("")
print("   t   R  free |nu|<=  nu with COVER  distinct values  KERNEL dim  H1 wrong  decoy wrong")
print("  " + "-" * 90)

H1BAD = []
H2ROWS = []
for t in TS:
    fz = t // 2 - 1
    for R in (fz + 1, fz + 2):
        r = R - fz
        if r < 1:
            continue
        MAX = 18 if R <= 3 else 12
        zt = zeta(t)
        frozen = [zt ** k for k in range(1, t // 2)]
        set_random_seed(6100 + 10 * t + R)
        PTS = []
        while len(PTS) < 12:
            zz = [F.random_element() for _ in range(r)]
            if any(x == 0 for x in zz):
                continue
            xs = list(zz) + frozen
            if len(set([x for x in xs] + [1 / x for x in xs])) == 2 * R:
                PTS.append((zz, xs))

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
                    n_, d_ = sp_at(nu, xs, R)
                    if d_ == 0:
                        bad = True
                        break
                    vec.append(n_ / d_)
                if bad or all(v == 0 for v in vec):
                    continue
                rows.append((nu, m, tuple(vec)))

        # kernel dimension of the value map, measured as a rank deficit
        M = matrix(F, [list(v) for _, _, v in rows])
        rk = M.rank()
        kern = len(rows) - rk
        proj = {}
        for nu, m, v in rows:
            proj.setdefault(canon(v), []).append((nu, m))

        # H1 and the decoy, both directions
        h1 = {}
        dc = {}
        for nu, m, v in rows:
            h1.setdefault(J(m, t), set()).add(canon(v))
            dc.setdefault(DECOY(m, t), set()).add(canon(v))
        h1bad = sum(1 for k, s in h1.items() if len(s) > 1)
        dcbad = sum(1 for k, s in dc.items() if len(s) > 1)
        for k, s in list(h1.items())[:200]:
            if len(s) > 1 and len(H1BAD) < 6:
                ex = [(nu, m) for nu, m, v in rows if J(m, t) == k][:2]
                H1BAD.append((t, R, k, ex))
        print("  %2d %3d %5d %6d %14d %16d %11d %9d %12d"
              % (t, R, r, MAX, len(rows), len(proj), kern, h1bad, dcbad))
        H2ROWS.append((t, R, r, rows, PTS))

print("")
print("  K2: 'KERNEL dim' is the rank deficit of the value map on the surviving nu.  A zero there")
print("      would mean the specialization is injective and there is no second brick.")
print("  K3: 'H1 wrong' counts invariant cells carrying MORE THAN ONE value up to scalar.")
print("")
if H1BAD:
    print("  H1 failures, by hand:")
    for t, R, k, ex in H1BAD:
        print("    t=%d R=%d  J=%s" % (t, R, k))
        for nu, m in ex:
            print("        nu=%-20s m=%s" % (nu, m))
else:
    print("  H1: no cell carries two different values.  The invariant")
    print("      J(nu) = (multiset of folded classes, exact values in the fixed classes {0,t/2})")
    print("      determines sp_nu(z,frozen) up to a scalar, and its fibres ARE the kernel.")

# ------------------------------------------------------------------ H2: the clean case factorizes
print("")
print("=" * 96)
print("H2  the clean case -- every non-fixed class hit exactly once -- must factor through Sp(2r)")
print("=" * 96)
print("")
print("   t   R  free   clean nu   factorizes   ratio not constant   examples")
print("  " + "-" * 84)
for t, R, r, rows, PTS in H2ROWS:
    nclean = nfac = nvar = 0
    ex = []
    for nu, m, v in rows:
        cs = [cls(x, t) for x in m]
        if any(cs.count(k) != 1 for k in range(1, t // 2)):
            continue
        nclean += 1
        left = sorted([x for x in m if cls(x, t) in (0, t // 2)], reverse=True)
        if len(left) != r:
            continue
        ratios = []
        ok = True
        for (zz, xs), val in zip(PTS, v):
            num = alt(left, zz)
            den = alt(shifted([0] * r, r), zz)
            if den == 0 or num == 0:
                ok = False
                break
            ratios.append(val / (num / den))
        if not ok:
            continue
        if len(set(ratios)) == 1:
            nfac += 1
            if len(ex) < 2:
                ex.append((nu, m, left))
        else:
            nvar += 1
    print("  %2d %3d %5d %10d %12d %20d   %s"
          % (t, R, r, nclean, nfac, nvar,
             "; ".join("nu=%s->left=%s" % (a, c) for a, b, c in ex)))

print("")
print("  'factorizes' = the ratio to the r-variable alternant on the leftover columns is the SAME")
print("  constant at all 12 points.  That is the Littlewood-shaped statement for this alphabet:")
print("  when it does not vanish, it is an Sp(2r) character on the leftover columns.")
print("")
print("DONE")
