"""Does the t-core decide the vanishing?  This is the question behind Remark 8.8; the
acknowledgements say whose it is.

Three questions, all decided by exhaustive enumeration rather than by argument:

  (Q1) At fixed N, does the residue profile (n_0,...,n_{t-1}) determine core_t(lambda)?
       If yes, then EVERY condition on the profile -- in particular "n_i = 0 for some i",
       which is branch (a) / Theorem 3.1(i) -- IS a condition on the t-core.
  (Q2) Which core does a degenerate profile give? (Answer: never the empty one.)
  (Q3) Does the core determine whether Psi_r(lambda) = 0?
       (Answer: no. The second branch -- d_3 = 0, i.e. self-complementary of odd width
       at t = 2 -- is invisible to the profile, hence to the core.)

Q3 is the point: core -> vanishing is not a function, so no core-based restatement of the
criterion can exist, however the first branch is dressed up.

Exact arithmetic throughout (Fraction); Schur values by the bialternant at three generic
rational points, so "vanishes" means "vanishes at all three".
Authors: Carles Marin, Claude (AI assistant)."""
from fractions import Fraction as F
from collections import defaultdict


def beta(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - 1 - j for j in range(N)]


def profile(b, t):
    return tuple(sum(1 for x in b if x % t == i) for i in range(t))


def core_from_beta(b, t):
    """t-core by lowering beta numbers into free slots; independent of anc/ak53_consistency.py,
    which packs residue classes instead. The two agree on everything checked below."""
    b = sorted(b, reverse=True)
    moved = True
    while moved:
        moved = False
        occupied = set(b)
        for i, x in enumerate(b):
            if x - t >= 0 and (x - t) not in occupied:
                b[i] = x - t
                b = sorted(b, reverse=True)
                moved = True
                break
    N = len(b)
    return tuple(p for p in (b[j] - (N - 1 - j) for j in range(N)) if p > 0)


def det(M):
    M = [row[:] for row in M]
    n = len(M)
    d = F(1)
    for i in range(n):
        p = next((k for k in range(i, n) if M[k][i] != 0), None)
        if p is None:
            return F(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        d *= M[i][i]
        inv = F(1) / M[i][i]
        for k in range(i + 1, n):
            f = M[k][i] * inv
            if f:
                for j in range(i, n):
                    M[k][j] -= f * M[i][j]
    return d


def schur(lam, alph):
    N = len(alph)
    b = beta(lam, N)
    return det([[a ** b[j] for j in range(N)] for a in alph]) / \
           det([[a ** (N - 1 - j) for j in range(N)] for a in alph])


def branch_a(b, t=2):
    return len({x % t for x in b}) == 1


def branch_b(lam, N):
    l = list(lam) + [0] * (N - len(lam))
    w = {l[i] + l[N - 1 - i] for i in range(N)}
    return len(w) == 1 and w.pop() % 2 == 1


def parts(nmax, N):
    out = [()]

    def rec(pref, rem, mx):
        if pref:
            out.append(tuple(pref))
        if len(pref) >= N:
            return
        for k in range(min(mx, rem), 0, -1):
            rec(pref + [k], rem - k, k)

    rec([], nmax, nmax)
    return sorted(set(out), key=lambda p: (sum(p), p))


# ---------------------------------------------------------------- Q1/Q2, all t
print("== Q1, Q2: the residue profile determines the t-core (N = t+2, |lambda| <= 24) ==")
for t in range(2, 6):
    N = t + 2
    seen, clash = {}, 0
    for lam in parts(24, N):
        b = beta(lam, N)
        p, c = profile(b, t), core_from_beta(b[:], t)
        if p in seen and seen[p] != c:
            clash += 1
        seen.setdefault(p, c)
    deg = [(p, c) for p, c in sorted(seen.items()) if 0 in p]
    empty = [p for p, c in seen.items() if c == ()]
    print(f"  t={t} N={N}: {len(seen)} profiles, profile->core clashes {clash}; "
          f"empty core <=> profile {empty}; "
          f"{len(deg)} degenerate profiles, all with NON-empty core: {all(c != () for _, c in deg)}")
    if t == 2:
        for p, c in sorted(seen.items()):
            print(f"      {p} -> core {c}" + ("   <- some n_i = 0" if 0 in p else ""))

# ---------------------------------------------------------------- Q3, t = 2
print()
print("== Q3: the core does NOT determine the vanishing of Psi_r (t = 2) ==")
ZS = [F(3, 2), F(5, 3), F(7, 4)]
WS = [F(11, 5), F(13, 6), F(17, 7)]
US = [F(19, 8), F(23, 9), F(29, 10)]
for r, NMAX in ((1, 14), (2, 16), (3, 14)):
    N = 2 * r + 2
    alphs = []
    for k in range(3):
        a = [ZS[k], 1 / ZS[k]]
        if r >= 2:
            a += [WS[k], 1 / WS[k]]
        if r >= 3:
            a += [US[k], 1 / US[k]]
        alphs.append(a + [F(1), F(-1)])
    rows = []
    for lam in parts(NMAX, N):
        b = beta(lam, N)
        rows.append((lam, profile(b, 2), core_from_beta(b[:], 2), branch_a(b), branch_b(lam, N),
                     all(schur(lam, a) == 0 for a in alphs)))
    bad_a = sum(1 for _, p, _, a, _, _ in rows if a != (0 in p))
    bad_c = [x for x in rows if x[5] != (x[3] or x[4])]
    g = defaultdict(list)
    for lam, p, c, a, bb, v in rows:
        g[(p, c)].append((lam, v))
    split = [(k, vs) for k, vs in g.items() if len({v for _, v in vs}) > 1]
    print(f"  r={r} N={N} |lambda|<={NMAX}: {len(rows)} shapes")
    print(f"     branch (a) <=> some n_i = 0            : violations {bad_a}")
    print(f"     (a) or (b) <=> Psi_r = 0               : violations {len(bad_c)}")
    print(f"     core classes carrying BOTH behaviours  : {len(split)}")
    for k, vs in split:
        van = [lam for lam, v in vs if v][:3]
        non = [lam for lam, v in vs if not v][:3]
        print(f"        profile {k[0]} core {k[1] if k[1] else '(empty)'}: "
              f"vanish {van} ... | non-vanish {non} ...")
print()
print("Conclusion: branch (a) is a core condition (Q1), and never the empty core (Q2);")
print("the criterion as a whole is not (Q3), because branch (b) does not move the core.")
