# -*- coding: utf-8 -*-
# Why the root count dies at r = 2 and the covering condition does not.
#
# THE DIAGNOSIS, WRITTEN BEFORE THE RUN.
#
# NPP's count is the right criterion for a TORSION element, where every row of the Weyl numerator is
# frozen: two weights in the same folded class then make two rows agree up to sign and the character
# dies.  Our element has r FREE rows, and a free row sees the exact value of a column, not its class.
# So a collision of classes is fatal only if the FROZEN block has to absorb it.  The frozen block is
# t/2-1 rows and it needs one column from each non-fixed class 1..t/2-1; the other r columns go to
# the free block, which is generic and does not care what classes they carry.
#
# Hence the correct condition is a MATCHING condition, not a counting one:
#
#     sp_nu(free, frozen) != 0  <=>  the columns admit a system of distinct representatives of the
#                                    classes 1..t/2-1   <=>   every such class is hit   <=>  COVER,
#
# and Hall's condition here degenerates to "no class is missed" because each column carries exactly
# one class.  NPP's count is the shadow of this at r = 0, where the matching is forced to be a
# bijection and "missing a class" and "colliding in a class" are the same event.
#
# THE ARITHMETIC OF WHY r = 2 IS THE FIRST PLACE THEY PART.  With R = r + t/2 - 1 columns and
# t/2 - 1 non-fixed classes to cover, there are r SPARE columns.
#   r = 0 : no spare.  A missed class and a collision are the same thing.  Count = matching.
#   r = 1 : one spare.  It can sit in a fixed class (count contribution 1, same as rho) or double a
#           non-fixed class (contribution C(2,2) = 1, again the same).  Both give N = N_0, so the
#           count STILL agrees -- by an accident of small numbers, not by a mechanism.
#   r = 2 : two spares.  Put BOTH in class 0.  COVER still holds, the frozen block still finds its
#           system of representatives, the free 2x2 block absorbs the two colliding columns and is
#           generically nonzero -- so sp_nu != 0.  But NPP's count sees n_0 = 2, which contributes
#           2*C(2,2) + 2 = 4 against N_0 = 2.  The count says zero; the determinant says nonzero.
#
# PREDICTIONS, all falsifiable:
#   P1  every disagreement between the NPP count and the determinant has a class with n_k >= 2.
#   P2  the minimal witness at each t has two columns in the SAME FIXED class (0 or t/2), the
#       cheapest way to spend two spares.
#   P3  restricted to nu with all n_k <= 1, the NPP count and COVER agree at EVERY r.  This is the
#       control that isolates the collisions as the whole difference.
#   P4  the number of disagreements must be 0 at r <= 1 and grow with r.
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


def cls(m, t):
    return min(m % t, (-m) % t)


def profile(m, t):
    n = {}
    for x in m:
        n[cls(x, t)] = n.get(cls(x, t), 0) + 1
    return n


def COVER(m, t):
    n = profile(m, t)
    return all(n.get(k, 0) >= 1 for k in range(1, t // 2))


def npp(m, t):
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


def collide(m, t):
    return any(v >= 2 for v in profile(m, t).values())


print("")
print("=" * 98)
print("P1/P4  where the count and the determinant part company, and what those shapes look like")
print("=" * 98)
print("")
print("   t   R free |nu|<=  tested   NPP wrong   of those, some n_k>=2   COVER wrong   P3 wrong")
print("  " + "-" * 92)

MINW = []
for t in TS:
    fz = t // 2 - 1
    for r in (0, 1, 2, 3):
        R = fz + r
        if R < 1:
            continue
        MAX = 16 if R <= 3 else 12
        zt = zeta(t)
        frozen = [zt ** k for k in range(1, t // 2)]
        rho0 = shifted([0] * R, R)
        N0 = npp(rho0, t)
        set_random_seed(3300 + 10 * t + R)
        PTS = []
        tries = 0
        while len(PTS) < 3 and tries < 500:
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
        n = wn = wn_coll = wc = p3 = n_p3 = 0
        best = None
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
                v = all(x == 0 for x in vals)           # True = vanishes
                if COVER(m, t) == v:
                    wc += 1
                if (npp(m, t) == N0) == v:
                    wn += 1
                    if collide(m, t):
                        wn_coll += 1
                    if best is None or size < best[0]:
                        best = (size, nu, m, profile(m, t))
                if not collide(m, t):                   # P3: no collisions at all
                    n_p3 += 1
                    if (npp(m, t) == N0) == v:
                        p3 += 1
        print("  %2d %3d %4d %6d %7d %11d %22d %13d %10d/%d"
              % (t, R, r, MAX, n, wn, wn_coll, wc, p3, n_p3))
        if best:
            MINW.append((t, r, best))

print("")
print("  'NPP wrong' counts disagreements with the determinant.  P1 says the next column must")
print("  equal it: every disagreement has a repeated folded class.  P4 says the column is 0 for")
print("  r <= 1.  'P3 wrong' is over the nu with NO repeated class: it must be 0/N there.")
print("")
print("=" * 98)
print("P2  the smallest shape on which the count is wrong, at each (t, r)")
print("=" * 98)
for t, r, (size, nu, m, prof) in MINW:
    fixed = [k for k in prof if k in (0, t // 2) and prof[k] >= 2]
    rep = [k for k, v in prof.items() if v >= 2]
    print("  t=%2d r=%d  |nu|=%d  nu=%-18s m=%s  profile=%s  repeated=%s%s"
          % (t, r, size, nu, m, dict(sorted(prof.items())), rep,
             "   <-- repeat in a FIXED class" if fixed else ""))
print("")
print("DONE")
