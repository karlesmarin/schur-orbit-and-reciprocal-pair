# -*- coding: utf-8 -*-
# NPP on the det = -1 component, second attempt: the count runs over the FOLDED residues.
#
# *** THE FIRST PROPOSAL WAS THE DEGENERATE CASE, AND THE RUN SAID SO. ***
# npp_minus.sage proposed the LONG-ROOT count: sp_nu dies at the frozen roots when a whole row of
# the symplectic bialternant dies, i.e. t | 2k*m_j for every j.  K1 held (0 failures, it is forced)
# but K2 -- the converse -- failed massively for t >= 6: 90 of 92 at t=6, 85 of 87 at t=8, 128 of
# 128 at t=10.  It was exact only at t=4, and for a reason with no content: at t=4 there is only
# ONE frozen row, so "dependent" and "zero" coincide.  The residue printed there shows what really
# happens: the frozen rows do not die, they become LINEARLY DEPENDENT.
#
# THE CORRECTED STATEMENT.  The frozen rows are indexed by k = 1..t/2-1 and have entries
#     v_k(m_j) = zeta^{k m_j} - zeta^{-k m_j},
# which depend on m_j only through m_j mod t and are ODD in m_j.  The free rows are generic, so the
# determinant vanishes identically in the free variables exactly when the frozen rows are dependent,
# i.e. when some nonzero ODD function on Z/t vanishes at every m_j.  Odd functions on Z/t vanish at
# 0 and at t/2 automatically and are free elsewhere, so such a function exists iff the m_j miss one
# of the non-fixed classes.  Writing cls(m) = min(m mod t, -m mod t) for the FOLDED residue:
#
#     sp_nu(free, zeta, .., zeta^{t/2-1}) != 0   <=>   {cls(m_j)} contains all of {1, .., t/2-1}.
#
# This is Littlewood's shape -- "no residue class is empty" -- but on the residues FOLDED by
# negation, which is the orbit space the paper's Lemma AtoSp already points at when it says the
# orbit Lie algebra of D_{r+1} is C_r.  The long-root proposal is the special case where the m_j
# miss EVERY non-fixed class at once.
#
# CONTROLS able to fail:
#   C1  COVER must agree with the measured vanishing in BOTH directions, at every t and R.
#   C2  the old LONG predicate is scored beside it and must do strictly worse for t >= 6, or the
#       run does not distinguish the two proposals.
#   C3  a decoy that must do worse: "the m_j miss some class including the fixed ones {0, t/2}".
#   C4  acceptance, fatal: sp of the empty partition is 1 and sp_(1) is sum(x + 1/x).
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
    assert (p - 1) % t == 0, "GUARD: %d does not divide p-1" % t
print("field GF(%d); every t in %s divides p-1  ->  guard PASS" % (p, TS))


def zeta(t):
    z = G ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t))
    return z


def shifted(nu, R):
    return [nu[j] + R - j for j in range(R)]


def sp_value(nu, xs, R):
    m = shifted(nu, R)
    d = shifted([0] * R, R)
    num = matrix(F, R, R, lambda a, b: xs[a] ** m[b] - xs[a] ** (-m[b])).det()
    den = matrix(F, R, R, lambda a, b: xs[a] ** d[b] - xs[a] ** (-d[b])).det()
    return num, den


def _accept():
    xs = [F(7), F(11)]
    n0, d0 = sp_value([0, 0], xs, 2)
    n1, d1 = sp_value([1, 0], xs, 2)
    ok = (n0 / d0 == 1) and (n1 / d1 == xs[0] + 1 / xs[0] + xs[1] + 1 / xs[1])
    print("acceptance: sp_empty = 1 and sp_(1) = sum(x+1/x)  ->  %s" % ("PASS" if ok else "FAIL"))
    if not ok:
        raise SystemExit(1)


_accept()


def cls(m, t):
    a = m % t
    return min(a, (-m) % t)


def COVER(nu, R, t):
    """True when sp_nu is predicted NONZERO: every non-fixed folded class is hit."""
    hit = set(cls(m, t) for m in shifted(nu, R))
    return all(k in hit for k in range(1, t // 2))


def LONG(nu, R, t):
    m = shifted(nu, R)
    return any(all((2 * k * mj) % t == 0 for mj in m) for k in range(1, t // 2))


def DECOY(nu, R, t):
    """miss some class counting the fixed ones too -- must do worse."""
    hit = set(cls(m, t) for m in shifted(nu, R))
    return all(k in hit for k in range(0, t // 2 + 1))


print("")
print("=" * 96)
print("THE FOLDED-RESIDUE COUNT AGAINST THE MEASURED VANISHING OF sp_nu")
print("=" * 96)
print("")
print("   t   R  free |nu|<=   tested   sp=0   COVER wrong   LONG wrong   decoy wrong")
print("  " + "-" * 84)

BAD = []
for t in TS:
    fz = t // 2 - 1
    for R in (fz + 1, fz + 2):
        r = R - fz
        if r < 1:
            continue
        MAX = 16 if R <= 3 else (14 if R <= 5 else 10)
        zt = zeta(t)
        frozen = [zt ** k for k in range(1, t // 2)]
        set_random_seed(5000 + 10 * t + R)
        PTS = []
        tries = 0
        while len(PTS) < 3 and tries < 300:
            tries += 1
            zz = [F.random_element() for _ in range(r)]
            if any(x == 0 for x in zz):
                continue
            xs = list(zz) + frozen
            if len(set([x for x in xs] + [1 / x for x in xs])) == 2 * R:
                PTS.append(xs)
        if len(PTS) < 3:
            print("  %2d %3d : could not build points -- SKIPPED" % (t, R))
            continue
        n = nz = wc = wl = wd = 0
        for size in range(MAX + 1):
            for l in Partitions(size, max_length=R):
                nu = list(l) + [0] * (R - len(l))
                vals = []
                ok = True
                for xs in PTS:
                    num, den = sp_value(nu, xs, R)
                    if den == 0:
                        ok = False
                        break
                    vals.append(num / den)
                if not ok:
                    continue
                n += 1
                v = all(x == 0 for x in vals)       # True = vanishes
                if v:
                    nz += 1
                if COVER(nu, R, t) == v:            # COVER predicts NONzero
                    wc += 1
                    if len(BAD) < 10:
                        BAD.append((t, R, nu, shifted(nu, R), v))
                if (not LONG(nu, R, t)) == v:
                    wl += 1
                if DECOY(nu, R, t) == v:
                    wd += 1
        print("  %2d %3d %5d %6d %8d %6d %13d %12d %12d"
              % (t, R, r, MAX, n, nz, wc, wl, wd))

print("")
print("  'wrong' counts disagreements with the determinant.  COVER must be 0 in every row.")
print("  LONG is the retracted proposal and must be large for t >= 6.  The decoy must be large.")
print("")
if BAD:
    print("  COVER failures, by hand:")
    for t, R, nu, m, v in BAD:
        print("    t=%d R=%d nu=%s  m=%s  vanishes=%s  classes=%s"
              % (t, R, nu, m, v, sorted(set(cls(x, t) for x in m))))
else:
    print("  COVER: no failure anywhere.  The criterion for the frozen part is")
    print("     sp_nu != 0  <=>  the folded residues of nu + rho_C cover {1, .., t/2-1}.")
print("")
print("DONE")
