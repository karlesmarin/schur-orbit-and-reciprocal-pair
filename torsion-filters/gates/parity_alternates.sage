# -*- coding: utf-8 -*-
# Does the pairwise cancellation alternate with the parity of r?
#
# The prediction comes from crossing two things we already have.  A(C-T) = (-1)^r A(T), so the
# reflection pairs terms that cancel only when the Laplace signs cooperate, and what "cooperate"
# means FLIPS with r:
#
#     r even :  A(C-T) = +A(T)   ->  cancels pairwise iff  eps = sgn(C-T)*sgn(T) = -1
#     r odd  :  A(C-T) = -A(T)   ->  cancels pairwise iff  eps = +1
#
# At r = 2, eps = -1 is rare -- 2 of 30 vanishing shapes -- which is why pairwise cancellation
# explains almost none of the interior zeros.  If the reading is right, at r = 3 the SAME mechanism
# should come back with the opposite sign, and most vanishing shapes should have eps = +1.  If
# instead r = 3 looks like r = 2, the parity reading is wrong and the two lines of the map differ
# for some other reason.
#
# ACCEPTANCE TEST first, fatal: lambda = (5,4,3) at r = 2 must give 8 terms, eps = -1, total zero.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

t = 4
PRIME = 998244353
F = GF(PRIME)
I4 = F(-1).sqrt()
ROOTS = [F(1), I4, F(-1), -I4]


def make_points(r, N, seed):
    set_random_seed(seed)
    pts = []
    while len(pts) < 3:
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        alph = ROOTS + [y for x in zz for y in (x, 1 / x)]
        if len(set(alph)) == N:
            pts.append(alph)
    return pts


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def split(beta):
    cls = {}
    for k, b in enumerate(beta):
        cls.setdefault(b % t, []).append(k)
    E = sorted(i for v in cls.values() if len(v) >= 2 for i in v)
    return cls, E


def inv_sign(pick, T):
    return (-1) ** sum(1 for i in pick for j in T if i > j)


def eps_of(beta, cls, C, N):
    pos = dict((b, i) for i, b in enumerate(beta))
    eps = None
    for pick in itertools.product(*[cls[k] for k in sorted(cls)]):
        pick = tuple(sorted(pick))
        T = tuple(sorted(set(range(N)) - set(pick)))
        try:
            T2 = tuple(sorted(pos[C - beta[i]] for i in T))
        except KeyError:
            return None
        p2 = tuple(sorted(set(range(N)) - set(T2)))
        e = inv_sign(pick, T) * inv_sign(p2, T2)
        if eps is None:
            eps = e
        elif eps != e:
            return 0
    return eps


def centres(beta, E):
    S = set(beta[i] for i in E)
    if not S:
        return []
    return [C for C in range(2 * max(S) + 1) if set(C - b for b in S) == S]


# ---- acceptance test ----------------------------------------------------------------------
b0 = beta_of([5, 4, 3], 8)
cls0, E0 = split(b0)
e0 = eps_of(b0, cls0, 12, 8)
n0 = 1
for v in cls0.values():
    n0 *= len(v)
print("ACCEPTANCE: beta=%s  terms=%d (want 8)  eps=%s (want -1)  ->  %s"
      % (b0, n0, e0, "PASS" if (n0 == 8 and e0 == -1) else "FAIL"))
if not (n0 == 8 and e0 == -1):
    raise SystemExit(1)

print("")
print("=" * 78)
print("Pairwise cancellation against the parity of r,  t = 4")
print("=" * 78)
print("")
print("  r  N  |lam|<=  shapes  zeros   eps=-1  eps=+1  eps=0   pairwise ones")
print("  " + "-" * 72)

for r, MAX in ((2, 24), (3, 20)):
    N = t + 2 * r
    pts = make_points(r, N, 4242 + r)

    def vanishes(beta):
        for alph in pts:
            if matrix(F, N, N, lambda a, b: alph[b] ** beta[a]).det() != 0:
                return False
        return True

    nsh = 0
    zeros = []
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            cls, E = split(beta)
            if len(cls) < t:
                continue
            nsh += 1
            if vanishes(beta):
                zeros.append((list(l), beta, cls, E))

    cm = cp = cz = cnone = 0
    for lam, beta, cls, E in zeros:
        got = None
        for C in centres(beta, E):
            e = eps_of(beta, cls, C, N)
            if e is not None:
                got = e
                break
        if got is None:
            cnone += 1
        elif got == -1:
            cm += 1
        elif got == 1:
            cp += 1
        else:
            cz += 1
    # the prediction: pairwise cancellation needs eps = -1 for r even, eps = +1 for r odd
    pw = cm if r % 2 == 0 else cp
    print("  %d %2d %7d %7d %6d %8d %7d %6d %14d"
          % (r, N, MAX, nsh, len(zeros), cm, cp, cz, pw))

print("")
print("  'pairwise ones' counts the zeros whose reflection cancels TERM BY TERM:")
print("     r even needs eps = -1 because A(C-T) = +A(T);")
print("     r odd  needs eps = +1 because A(C-T) = -A(T).")
print("  The prediction is that this column is small for r=2 and large for r=3.")
print("")
print("DONE")
