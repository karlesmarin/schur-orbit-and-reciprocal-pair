# -*- coding: utf-8 -*-
# LITERATURE GATE on the folded lemma.  Is it Nadimpalli-Pattanayak-Prasad in type C?
#
# THE LEMMA.  sp_nu(z_1..z_r, zeta, .., zeta^{t/2-1}) != 0  <=>  the folded residues
# cls(m_j) = min(m_j mod t, -m_j mod t) of m_j = nu_j + R - j + 1 cover {1,..,t/2-1}.
#
# THE SUSPICION THAT MUST BE SETTLED FIRST.  For Sp(2n) the Weyl vector is rho = (n, n-1, .., 1), so
# the principal torsion element rho^(e^{2 pi i/t}) has eigenvalues zeta^{n}, .., zeta^{1} and their
# inverses.  With n = t/2 - 1 that multiset is exactly {zeta^k, zeta^{-k} : 1 <= k <= t/2-1}
# = mu_t minus {1,-1} -- OUR FROZEN PART.  So the frozen block is not an arbitrary choice: it is the
# principal element of order t of Sp(t-2), and at r = 0 the lemma is a statement NPP already make.
#
# NPP Corollary 3.7 (arXiv:2504.14684), for connected reductive G at C_m = rho^(e^{2 pi i/m}):
#     chi_nu(C_m) != 0  <=>  #{a>0 : m | <rho^,a>} = #{a>0 : m | <nu^+rho^,a>}.
# In type C_n the positive roots are e_i-e_j, e_i+e_j (i<j) and 2e_i, so with beta = nu+rho the
# right-hand count is
#     #{i<j : t | m_i-m_j} + #{i<j : t | m_i+m_j} + #{i : t | 2m_i}
#   = #{pairs with EQUAL folded class} + #{i with folded class in {0, t/2}},
# because t | m_i-m_j or t | m_i+m_j is exactly cls(m_i) = cls(m_j), and t | 2m_i is exactly
# cls(m_i) in {0, t/2}.  For rho itself the m are (n, .., 1) with n = t/2-1: all folded classes
# distinct, none fixed, so the left count is 0.
#
# WHAT IS TESTED:
#   G1  at r = 0 (n = t/2-1 variables, all frozen): does COVER agree with NPP's count equality,
#       and do BOTH agree with the determinant?  If yes, the r = 0 case is theirs, not ours.
#   G2  at r >= 1: NPP does not apply -- the element is not torsion -- but the COUNT is still
#       defined on the frozen block.  Is COVER still the NPP count equality?  This measures how
#       much of the lemma is the free-variable extension.
#   G3  a control that can fail: the type-A count (only e_i-e_j, i.e. forgetting e_i+e_j and 2e_i)
#       must do WORSE, or the run does not show that the folding is what matters.
#
# Authors: Carles Marin, Claude (AI assistant).

TS = [4, 6, 8, 10, 12]
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


def cls(m, t):
    return min(m % t, (-m) % t)


def COVER(m, t):
    h = set(cls(x, t) for x in m)
    return all(k in h for k in range(1, t // 2))


def npp_count_C(m, t):
    """#{i<j : t | m_i-m_j} + #{i<j : t | m_i+m_j} + #{i : t | 2 m_i}, the type-C count."""
    n = len(m)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (m[i] - m[j]) % t == 0:
                c += 1
            if (m[i] + m[j]) % t == 0:
                c += 1
        if (2 * m[i]) % t == 0:
            c += 1
    return c


def npp_count_A(m, t):
    """the type-A count only, the control."""
    n = len(m)
    return sum(1 for i in range(n) for j in range(i + 1, n) if (m[i] - m[j]) % t == 0)


# ---- the identification of the frozen block, checked and not assumed
print("")
print("G0  the frozen block IS the principal torsion element of Sp(t-2)")
for t in TS:
    n = t // 2 - 1
    zt = zeta(t)
    ours = sorted([(zt ** k, zt ** (-k)) for k in range(1, t // 2)])
    rho = shifted([0] * n, n)                      # (n, n-1, .., 1)
    theirs = sorted([(zt ** e, zt ** (-e)) for e in rho])
    same = set([zt ** k for k in range(1, t // 2)] + [zt ** (-k) for k in range(1, t // 2)]) \
        == set([zt ** e for e in rho] + [zt ** (-e) for e in rho])
    print("    t=%2d  n=%d  rho=%s  ->  %s" % (t, n, rho, "SAME multiset" if same else "DIFFERENT"))

print("")
print("=" * 96)
print("G1/G2  COVER  vs  NPP's type-C count equality  vs  the determinant")
print("=" * 96)
print("")
print("   t   R free |nu|<=  tested   sp=0   COVER wrong   NPP-C wrong   typeA wrong (control)")
print("  " + "-" * 88)

for t in TS:
    fz = t // 2 - 1
    for R in (fz, fz + 1, fz + 2):
        r = R - fz
        if R < 1:
            continue
        MAX = 16 if R <= 3 else 12
        zt = zeta(t)
        frozen = [zt ** k for k in range(1, t // 2)]
        rho0 = shifted([0] * R, R)
        base = npp_count_C(rho0, t)
        set_random_seed(2200 + 10 * t + R)
        PTS = []
        tries = 0
        while len(PTS) < 3 and tries < 400:
            tries += 1
            zz = [F.random_element() for _ in range(r)]
            if any(x == 0 for x in zz):
                continue
            xs = list(zz) + frozen
            if len(set([x for x in xs] + [1 / x for x in xs])) == 2 * R:
                PTS.append(xs)
        if len(PTS) < 3:
            print("  %2d %3d : no admissible points -- SKIPPED" % (t, R))
            continue
        n = nz = wc = wn = wa = 0
        baseA = npp_count_A(rho0, t)
        for size in range(MAX + 1):
            for l in Partitions(size, max_length=R):
                nu = list(l) + [0] * (R - len(l))
                m = shifted(nu, R)
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
                n += 1
                v = all(x == 0 for x in vals)          # True = vanishes
                if v:
                    nz += 1
                if COVER(m, t) == v:
                    wc += 1
                if (npp_count_C(m, t) == base) == v:
                    wn += 1
                if (npp_count_A(m, t) == baseA) == v:
                    wa += 1
        print("  %2d %3d %4d %6d %7d %6d %13d %13d %14d%s"
              % (t, R, r, MAX, n, nz, wc, wn, wa, "   <- r=0: NPP's own case" if r == 0 else ""))

print("")
print("  'wrong' = disagreements with the determinant.  If NPP-C is 0 at r = 0 the lemma's frozen")
print("  case is Corollary 3.7 of arXiv:2504.14684 in type C and must be cited as theirs.  If it is")
print("  also 0 for r >= 1, the whole lemma is their count and only the free-variable statement is")
print("  ours.  The type-A column is the control: it must be large, or the folding is not the point.")
print("")
print("DONE")
