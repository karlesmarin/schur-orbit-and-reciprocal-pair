# -*- coding: utf-8 -*-
# The interior is the t = 2 problem with pairs specialized to roots of unity.
#
# READING THE ALPHABET AGAIN.  For t EVEN,
#     mu_t = {1, zeta, ..., zeta^{t-1}} = {1, -1} u {zeta^k, zeta^{-k} : 1 <= k <= t/2 - 1},
# because zeta^{t/2} = -1 and the remaining powers pair off under inversion.  So mu_t is {1,-1}
# together with (t-2)/2 RECIPROCAL PAIRS, and therefore
#
#     Phi_t(lambda; z_1..z_r)  =  Psi_R(lambda) evaluated at (z_1,..,z_r, zeta, zeta^2, .., zeta^{t/2-1}),
#     with  R = r + (t-2)/2   and   N = t + 2r = 2R + 2.
#
# The whole (t,r) interior is one object -- Psi_R of Section 9 of arXiv:2608.09619 -- with some of
# its free pairs frozen at roots of unity.  Three consequences, all testable:
#
#   P1  the identity itself: Phi_t at (z) must EQUAL Psi_R at the specialized point, coefficient by
#       coefficient.  Fatal if it fails, because everything below rests on it.
#   P2  specializing can only CREATE zeros, never destroy them, so at fixed N the zero sets must be
#       NESTED INCREASING in t:   Z(t=2) subset Z(t=4) subset Z(t=6) subset ...
#       A single lambda vanishing at t=2 and not at t=4 kills the reading.
#   P3  for t ODD, prod(mu_t) = +1: the element sits in SO(N), the det-twist is trivial, and the
#       second branch cannot exist -- which is the representation-theoretic reason behind the
#       measured "odd t has no second branch", not a separate combinatorial accident.
#       Tested as: zeros at odd t, on shapes with every residue class occupied, must be 0.
#
# CONTROL able to fail: P2 is a containment, so the run also reports the STRICT growth.  If the
# sets were equal the reading would be true but empty -- the interior would add no zeros at all.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

PRIME = 2013265921
F = GF(PRIME)
G = F.multiplicative_generator()


def roots_of(t):
    z = G ** ((PRIME - 1) // t)
    return [z ** k for k in range(t)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def schur(beta, alph, N):
    return matrix(F, [[alph[b] ** beta[a] for b in range(N)] for a in range(N)]).det()


# =============================================================== P1: the alphabets are the same set
print("=" * 92)
print("P1  mu_t = {1,-1} u (t-2)/2 reciprocal pairs, and the two evaluations agree")
print("=" * 92)
bad = 0
for t in (2, 4, 6, 8, 10):
    RT = roots_of(t)
    z = RT[1]
    rebuilt = [F(1), F(-1)] + [y for k in range(1, t // 2) for y in (z ** k, z ** (-k))]
    same = set(RT) == set(rebuilt)
    print("  t=%2d : mu_t rebuilt from {1,-1} and %d pairs -> %s"
          % (t, t // 2 - 1, "SAME SET" if same else "DIFFERENT -- FATAL"))
    if not same:
        bad += 1
if bad:
    raise SystemExit(1)

# the evaluation identity on real shapes
print("")
print("  the evaluation identity Phi_t(z) = Psi_R(z, zeta, .., zeta^{t/2-1}), on 200 shapes each:")
for t, r in ((4, 2), (6, 2), (8, 1)):
    N = t + 2 * r
    RT = roots_of(t)
    zeta = RT[1]
    set_random_seed(70 + t)
    n = ok = 0
    for size in range(0, 14):
        for l in Partitions(size, max_length=N):
            if n >= 200:
                break
            beta = beta_of(list(l), N)
            zz = [F.random_element() for _ in range(r)]
            if any(x == 0 for x in zz):
                continue
            A_phi = RT + [y for x in zz for y in (x, 1 / x)]
            A_psi = [F(1), F(-1)] + [y for x in zz for y in (x, 1 / x)] \
                + [y for k in range(1, t // 2) for y in (zeta ** k, zeta ** (-k))]
            if len(set(A_phi)) < N or len(set(A_psi)) < N:
                continue
            n += 1
            # both are bialternants over the SAME multiset, so they agree up to the sign of the
            # column permutation; compare the Schur value, i.e. divide by the Vandermonde
            v1 = schur(beta, A_phi, N) / schur(beta_of([], N), A_phi, N)
            v2 = schur(beta, A_psi, N) / schur(beta_of([], N), A_psi, N)
            if v1 == v2:
                ok += 1
    print("    t=%d r=%d N=%2d : %3d / %3d shapes agree%s" % (t, r, N, ok, n,
          "" if ok == n else "   <-- FATAL"))
    if ok != n:
        raise SystemExit(1)


# ================================================================ P2: the zero sets nest in t
def zeros_at(t, N, seed):
    r = (N - t) // 2
    RT = roots_of(t)
    set_random_seed(seed)
    pts = []
    while len(pts) < 3:
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        al = RT + [y for x in zz for y in (x, 1 / x)]
        if len(set(al)) == N:
            pts.append(al)
    return pts, r


print("")
print("=" * 92)
print("P2  at fixed N the zero sets must be NESTED INCREASING in t")
print("=" * 92)

for N, MAX in ((10, 26), (12, 22)):
    print("")
    print("  N = %d,  |lambda| <= %d" % (N, MAX))
    print("     t   r   zeros   new vs previous t   LOST (must be 0)   all-classes-occupied shapes")
    print("    " + "-" * 84)
    prev = None
    prev_t = None
    for t in range(2, N - 1, 2):
        r = (N - t) // 2
        if r < 1:
            continue
        pts, r = zeros_at(t, N, 500 + t + N)
        EMAX = MAX + N
        POW = [[[a ** e for e in range(EMAX + 1)] for a in al] for al in pts]
        Z = set()
        nsh = 0
        for size in range(MAX + 1):
            for l in Partitions(size, max_length=N):
                lam = tuple(l)
                beta = beta_of(list(l), N)
                cls = {}
                for k, b in enumerate(beta):
                    cls.setdefault(b % t, []).append(k)
                if len(cls) == t:
                    nsh += 1
                v = True
                for p in range(3):
                    if matrix(F, [[POW[p][b][beta[a]] for b in range(N)]
                                  for a in range(N)]).det() != 0:
                        v = False
                        break
                if v:
                    Z.add(lam)
        if prev is None:
            print("    %2d %3d %7d %19s %18s %14d" % (t, r, len(Z), "-", "-", nsh))
        else:
            lost = prev - Z
            print("    %2d %3d %7d %19d %18d %14d"
                  % (t, r, len(Z), len(Z - prev), len(lost), nsh))
            for lam in sorted(lost)[:4]:
                print("        LOST at t=%d (was a zero at t=%d): lam=%s" % (t, prev_t, list(lam)))
        prev, prev_t = Z, t

# ================================================================ P3: odd t sits in SO(N)
print("")
print("=" * 92)
print("P3  for t odd  prod(mu_t) = +1, the element is in SO(N), and the second branch cannot exist")
print("=" * 92)
print("")
print("     t   prod(mu_t)   component     r   N  |lam|<=  shapes  ZEROS (all classes occupied)")
print("    " + "-" * 82)
for t, r, MAX in ((3, 2, 24), (5, 2, 20), (7, 2, 18), (4, 2, 24), (6, 2, 20)):
    RT = roots_of(t)
    pr = prod(RT)
    N = t + 2 * r
    pts, _ = zeros_at(t, N, 800 + t)
    EMAX = MAX + N
    POW = [[[a ** e for e in range(EMAX + 1)] for a in al] for al in pts]
    nsh = z = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            cls = {}
            for k, b in enumerate(beta):
                cls.setdefault(b % t, []).append(k)
            if len(cls) < t:
                continue
            nsh += 1
            v = True
            for p in range(3):
                if matrix(F, [[POW[p][b][beta[a]] for b in range(N)]
                              for a in range(N)]).det() != 0:
                    v = False
                    break
            if v:
                z += 1
    print("    %2d %12s   %-11s %3d %3d %7d %8d %6d"
          % (t, "+1" if pr == F(1) else "-1", "SO(N)" if pr == F(1) else "O(N) minus",
             r, N, MAX, nsh, z))

print("")
print("  The odd rows must show 0 zeros and the even rows must not.  That is the det-twist")
print("  argument of Section 9 of the paper, now as the reason the interior has one branch")
print("  for even t and none for odd t.")
print("")
print("DONE")
