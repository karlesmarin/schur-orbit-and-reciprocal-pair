# -*- coding: utf-8 -*-
# The specialization reading, with the two defects of gate_specialization.sage fixed.
#
# *** WHAT WENT WRONG THERE, LEFT ON THE RECORD ***
#
# (1) NO GUARD ON THE FIELD.  It used PRIME = 2013265921, whose p-1 = 3*5*2^27 is NOT divisible by
#     7, so G^((p-1)/7) is not a primitive 7th root of unity -- integer division truncated and the
#     row for t = 7 was computed on a fake alphabet.  It printed "prod = -1, O(N) minus, 0 zeros"
#     with perfect manners.  Here the field is chosen so that every t used divides p-1, and the
#     script ASSERTS it for each t before using it.
#
# (2) THE NESTING WAS STATED WRONG BY ME.  mu_t is contained in mu_{t'} exactly when t divides t',
#     not when t' is the next even number.  The old script compared consecutive t, so it tested
#     t=4 vs t=6 -- a containment nothing predicts -- and the 183 "LOST" shapes it reported refute
#     my sentence, not the reading.  The comparisons here are the ones divisibility actually
#     licenses.
#
# THE READING.  For t even, mu_t = {1,-1} u {zeta^k, zeta^{-k} : 1 <= k <= t/2 - 1}: the alphabet of
# Section 9 of arXiv:2608.09619 with (t-2)/2 of its free pairs frozen at roots of unity.  So
#     Phi_t(lambda; z_1..z_r) = Psi_R(lambda) at (z_1..z_r, zeta, .., zeta^{t/2-1}),  R = r+(t-2)/2,
# and the whole (t,r) interior is ONE object with more or fewer pairs frozen.  Freezing can only
# create zeros, so at fixed N:
#
#     P2   t | t'   ==>   Z(t, N)  subset  Z(t', N).      A lost shape kills the reading.
#
# and, since prod(mu_t) = +1 for t odd and -1 for t even, the element is in SO(N) for odd t and in
# the other component of O(N) for even t:
#
#     P3   t odd  ==>  no second branch at all, for any r -- the det-twist of Section 9 is trivial.
#
# CONTROL able to fail, both ways: P2 also reports the STRICT growth.  Containment with zero growth
# would make the reading true and empty.
#
# Authors: Carles Marin, Claude (AI assistant).

TS = [2, 3, 4, 5, 6, 7, 8, 9, 10]
L = lcm(TS)
p = next_prime(10 ** 9)
while (p - 1) % L != 0:
    p = next_prime(p)
F = GF(p)
G = F.multiplicative_generator()
print("field: GF(%d),  p-1 divisible by lcm%s = %d" % (p, tuple(TS), L))
for t in TS:
    assert (p - 1) % t == 0, "GUARD: %d does not divide p-1" % t
print("guard: every t in %s divides p-1  ->  PASS" % TS)


def roots_of(t):
    assert (p - 1) % t == 0
    z = G ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t)), "not primitive"
    return [z ** k for k in range(t)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def points(t, N, seed, howmany=3):
    r = (N - t) // 2
    RT = roots_of(t)
    set_random_seed(seed)
    pts = []
    tries = 0
    while len(pts) < howmany and tries < 500:
        tries += 1
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        al = RT + [y for x in zz for y in (x, 1 / x)]
        if len(set(al)) == N:
            pts.append(al)
    return pts, r


def zeros_at(t, N, MAX, seed):
    pts, r = points(t, N, seed)
    POW = [[[a ** e for e in range(MAX + N + 1)] for a in al] for al in pts]
    Z = set()
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            v = True
            for q in range(3):
                if matrix(F, [[POW[q][b][beta[a]] for b in range(N)]
                              for a in range(N)]).det() != 0:
                    v = False
                    break
            if v:
                Z.add(tuple(l))
    return Z, r


print("")
print("=" * 92)
print("P2  Z(t,N) subset Z(t',N) whenever t divides t'   [the containment divisibility licenses]")
print("=" * 92)
for N, MAX in ((10, 24), (12, 20)):
    print("")
    print("  N = %d,  |lambda| <= %d" % (N, MAX))
    Zs = {}
    for t in range(2, N - 1):
        if (N - t) % 2 or (N - t) // 2 < 1:
            continue
        Zs[t], _ = zeros_at(t, N, MAX, 600 + t + N)
    print("    zeros by t: %s" % ", ".join("t=%d:%d" % (t, len(Zs[t])) for t in sorted(Zs)))
    print("")
    print("      t   t'   |Z(t)|  |Z(t')|   LOST (must be 0)   GAINED (must be > 0)")
    print("    " + "-" * 72)
    for t in sorted(Zs):
        for t2 in sorted(Zs):
            if t2 > t and t2 % t == 0:
                lost = Zs[t] - Zs[t2]
                print("    %3d %4d %8d %8d %14d %20d%s"
                      % (t, t2, len(Zs[t]), len(Zs[t2]), len(lost), len(Zs[t2] - Zs[t]),
                         "   <-- P2 VIOLATED" if lost else ""))
                for lam in sorted(lost)[:3]:
                    print("        LOST: lam=%s" % list(lam))

print("")
print("=" * 92)
print("P3  prod(mu_t) decides the component, and the component decides whether a branch exists")
print("=" * 92)
print("")
print("     t  prod(mu_t)  component      r   N |lam|<=   shapes  ZEROS with every class occupied")
print("    " + "-" * 84)
for t, r, MAX in ((3, 2, 22), (5, 2, 18), (7, 2, 16), (9, 1, 16),
                  (4, 2, 22), (6, 2, 18), (8, 2, 16), (10, 1, 16)):
    RT = roots_of(t)
    pr = prod(RT)
    N = t + 2 * r
    pts, _ = points(t, N, 900 + t)
    POW = [[[a ** e for e in range(MAX + N + 1)] for a in al] for al in pts]
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
            for q in range(3):
                if matrix(F, [[POW[q][b][beta[a]] for b in range(N)]
                              for a in range(N)]).det() != 0:
                    v = False
                    break
            if v:
                z += 1
    print("    %3d %11s  %-12s %3d %3d %7d %8d %7d%s"
          % (t, "+1" if pr == F(1) else "-1", "SO(N)" if pr == F(1) else "O(N)-minus",
             r, N, MAX, nsh, z, "   <-- odd t with a zero!" if (t % 2 and z) else ""))

print("")
print("  Odd t must show 0 and even t must not.  prod(mu_t) = zeta^{t(t-1)/2} is +1 for t odd")
print("  (t divides t(t-1)/2) and -1 for t even (t(t-1)/2 = t/2 mod t), so the twist by det is")
print("  trivial exactly for odd t.  That is the reason, not a coincidence of the sweep.")
print("")
print("DONE")
