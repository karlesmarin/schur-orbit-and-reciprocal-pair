# -*- coding: utf-8 -*-
# Does the alternation survive thirty zeros?
#
# parity_alternates.sage measured, at t = 4:  r = 2 -> 3 of 30 zeros cancel pairwise,
# r = 3 -> 6 of 6.  Six is an anecdote.  This widens r = 3 until the zero count is well past
# thirty, and adds the controls the first pass did not have.
#
# THE STATEMENT BEING TESTED.  Psi ~ sum_T sgn(T) A(T) over complements of transversals, and the
# reflection T -> C - T obeys A(C-T) = (-1)^r A(T).  Pairing the terms,
#
#     sgn(T)A(T) + sgn(C-T)A(C-T)  =  sgn(T)A(T) * (1 + eps*(-1)^r),      eps = sgn(C-T)sgn(T),
#
# so the pair kills itself iff eps = (-1)^{r+1}.  That is a PROOF, not a measurement: whenever eps
# is constant over the terms and equals (-1)^{r+1}, Psi vanishes identically.  What is NOT proved,
# and is the whole content of the sweep, is the CONVERSE -- whether every zero arises this way.
# At r = 2 it does not (27 of 30 unexplained).  At r = 3 the six were all explained.  The question
# here is whether "r odd => the pairing explains every zero" survives a real sample.
#
# TWO INDEPENDENT ROUTES, forced to agree.  eps is computed
#   (A) brute force, enumerating every transversal and counting inversions;
#   (B) by a closed form derived separately: sgn(pick|rest) = (-1)^{sum(pick) - t(t-1)/2}, hence
#         eps = (-1)^{ const + sum_k [p_k - iota(p_k)] },   iota(i) = position of C - beta[i],
#       so eps is constant over terms iff i - iota(i) has one parity per excess class.
# Route B is O(N) and drives the wide sweep; route A audits it on every shape of the acceptance
# range.  If they ever disagree the run dies.
#
# CONTROLS, each able to fail:
#   C1  eps = (-1)^{r+1}  =>  vanishes.  This is forced by the algebra; a single exception means
#       the mechanism as written is wrong.  Must be 100%.
#   C2  concentric-and-non-vanishing must be NON-EMPTY, or concentricity alone decides the locus
#       and eps is idle decoration.
#   C3  no zero may be non-concentric (necessity, re-checked in the wider range).
#   C4  every zero re-verified at 6 further random points AND over a second prime.
#   C5  eps non-constant over terms is counted, not swept under the rug.
#
# ACCEPTANCE TESTS, fatal, run first:
#   A1  r=2, lambda=(5,4,3): 8 terms, eps = -1        (the case printed by hand in parity2)
#   A2  r=3, |lambda| <= 20: 2387 shapes, 6 zeros, all eps = +1  (parity_alternates_OUT.txt)
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools
import time

t = 4
P1 = 998244353
P2 = 1004535809          # second prime, also 1 mod 4, for the independent re-check of the zeros


def field_roots(prime):
    F = GF(prime)
    i = F(-1).sqrt()
    return F, [F(1), i, F(-1), -i]


def make_points(F, roots, r, N, seed, howmany):
    set_random_seed(seed)
    pts = []
    while len(pts) < howmany:
        zz = [F.random_element() for _ in range(r)]
        if any(x == 0 for x in zz):
            continue
        alph = roots + [y for x in zz for y in (x, 1 / x)]
        if len(set(alph)) == N:
            pts.append(alph)
    return pts


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def classes(beta):
    cls = {}
    for k, b in enumerate(beta):
        cls.setdefault(b % t, []).append(k)
    return cls


def centre(beta, cls):
    """The unique C with C - S = S on the excess values S, or None."""
    S = set(beta[i] for v in cls.values() if len(v) >= 2 for i in v)
    if not S:
        return None
    C = min(S) + max(S)
    return C if set(C - b for b in S) == S else None


# ---------------------------------------------------------------------------- route A: brute force
def shuffle_sign(pick, rest):
    perm = list(pick) + list(rest)
    sg = 1
    for a in range(len(perm)):
        for b in range(a + 1, len(perm)):
            if perm[a] > perm[b]:
                sg = -sg
    return sg


def eps_brute(beta, cls, C, N):
    pos = dict((b, i) for i, b in enumerate(beta))
    val = None
    for pick in itertools.product(*[cls[k] for k in sorted(cls)]):
        pick = tuple(sorted(pick))
        rest = tuple(sorted(set(range(N)) - set(pick)))
        try:
            rest2 = tuple(sorted(pos[C - beta[i]] for i in rest))
        except KeyError:
            return None
        pick2 = tuple(sorted(set(range(N)) - set(rest2)))
        e = shuffle_sign(pick, rest) * shuffle_sign(pick2, rest2)
        if val is None:
            val = e
        elif val != e:
            return 0
    return val


# ------------------------------------------------------------------- route B: the O(N) closed form
def eps_fast(beta, cls, C, N):
    pos = dict((b, i) for i, b in enumerate(beta))
    exc = [k for k in cls if len(cls[k]) >= 2]
    sing = [k for k in cls if len(cls[k]) == 1]
    iota = {}
    for k in exc:
        for i in cls[k]:
            j = pos.get(C - beta[i])
            if j is None:
                return None
            iota[i] = j
    const = sum(cls[k][0] for k in sing) + sum(iota[i] for k in exc for i in cls[k]) \
        + N * (N - 1) // 2
    tot = const
    for k in exc:
        par = set((i - iota[i]) % 2 for i in cls[k])
        if len(par) > 1:
            return 0                      # eps not constant over the terms
        tot += cls[k][0] - iota[cls[k][0]]
    return (-1) ** (tot % 2)


# ---------------------------------------------------------------------------------- acceptance A1
b0 = beta_of([5, 4, 3], 8)
c0 = classes(b0)
C0 = centre(b0, c0)
n0 = prod(len(v) for v in c0.values())
eA = eps_brute(b0, c0, C0, 8)
eB = eps_fast(b0, c0, C0, 8)
print("A1  beta=%s  C=%s  terms=%d(want 8)  epsA=%s epsB=%s (want -1)  ->  %s"
      % (b0, C0, n0, eA, eB, "PASS" if (n0 == 8 and eA == -1 and eB == -1) else "FAIL"))
if not (n0 == 8 and eA == -1 and eB == -1):
    raise SystemExit(1)


# ------------------------------------------------------------------------------------ the sweeper
def sweep(r, MAX, seed, brute_upto, label):
    N = t + 2 * r
    EMAX = MAX + N
    want = (-1) ** (r + 1)                        # the eps that makes the pairs cancel

    F1, R1 = field_roots(P1)
    pts = make_points(F1, R1, r, N, seed, 3)
    POW = [[[a ** e for e in range(EMAX + 1)] for a in al] for al in pts]

    def vanishes(beta):
        for p in range(len(POW)):
            M = matrix(F1, [[POW[p][b][beta[a]] for b in range(N)] for a in range(N)])
            if M.det() != 0:
                return False
        return True

    nsh = ncon = nnc = nagree = ndis = 0
    tab = {}                                       # (eps, vanishes) -> count
    zeros = []
    t0 = time.time()
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            lam = list(l)
            beta = beta_of(lam, N)
            cls = classes(beta)
            if len(cls) < t:
                continue
            nsh += 1
            C = centre(beta, cls)
            e = None if C is None else eps_fast(beta, cls, C, N)
            if size <= brute_upto:
                eb = None if C is None else eps_brute(beta, cls, C, N)
                if eb == e:
                    nagree += 1
                else:
                    ndis += 1
                    print("   ROUTE MISMATCH lam=%s C=%s  brute=%s fast=%s" % (lam, C, eb, e))
            if C is None or e is None:
                nnc += 1
            else:
                ncon += 1
            v = vanishes(beta)
            if v:
                zeros.append((lam, beta, cls, C, e))
            if C is not None and e is not None:
                tab[(e, v)] = tab.get((e, v), 0) + 1
    el = time.time() - t0

    expl = [z for z in zeros if z[4] == want]
    unex = [z for z in zeros if z[4] != want]
    return dict(r=r, N=N, MAX=MAX, want=want, nsh=nsh, ncon=ncon, nnc=nnc, tab=tab,
                zeros=zeros, expl=expl, unex=unex, el=el, nagree=nagree, ndis=ndis,
                label=label, pts=pts, POW=POW)


def report(S):
    print("")
    print("=" * 88)
    print("%s :  t = 4,  r = %d,  N = %d,  |lambda| <= %d      [%.1fs]"
          % (S['label'], S['r'], S['N'], S['MAX'], S['el']))
    print("=" * 88)
    print("  shapes with all %d classes occupied : %d" % (t, S['nsh']))
    print("  of them concentric (eps defined)    : %d" % S['ncon'])
    print("  not concentric / eps undefined      : %d" % S['nnc'])
    if S['nagree'] + S['ndis']:
        print("  route A vs route B on |lambda|<=%2d   : %d agree, %d DISAGREE"
              % (BRUTE, S['nagree'], S['ndis']))
    print("")
    print("  the contingency table over concentric shapes (eps x vanishes):")
    print("     eps      vanishes   count")
    for e in (1, -1, 0):
        for v in (True, False):
            n = S['tab'].get((e, v), 0)
            mark = ""
            if e == S['want'] and v is False:
                mark = "   <-- C1 VIOLATION, the mechanism is wrong"
            print("    %+3s        %-5s    %6d%s" % (e if e else "0(nc)", v, n, mark))
    print("")
    print("  ZEROS                               : %d" % len(S['zeros']))
    print("  explained by the pairing (eps=%+d)   : %d" % (S['want'], len(S['expl'])))
    print("  NOT explained                       : %d" % len(S['unex']))
    return S


# ---------------------------------------------------------------------------------- acceptance A2
BRUTE = 20
A2 = sweep(3, 20, 4242 + 3, 20, "A2 acceptance")
ok = (A2['nsh'] == 2387 and len(A2['zeros']) == 6 and len(A2['expl']) == 6 and A2['ndis'] == 0)
print("A2  r=3 |lam|<=20: shapes=%d(want 2387) zeros=%d(want 6) explained=%d(want 6) "
      "routes disagree=%d(want 0)  ->  %s"
      % (A2['nsh'], len(A2['zeros']), len(A2['expl']), A2['ndis'], "PASS" if ok else "FAIL"))
if not ok:
    raise SystemExit(1)


# ------------------------------------------------------------------------------- the wide r = 3
BRUTE = 20
S3 = report(sweep(3, 34, 4245, 20, "THE WIDE r = 3"))

print("")
print("  the unexplained zeros, by hand (empty is the claim):")
if not S3['unex']:
    print("    (none)")
for lam, beta, cls, C, e in S3['unex']:
    print("    lam=%-26s beta=%s  C=%s  sizes=%s  eps=%s"
          % (lam, beta, C, sorted(len(v) for v in cls.values()), e))

print("")
print("  five explained zeros, by hand, so the reader can check one:")
for lam, beta, cls, C, e in S3['zeros'][:5]:
    print("    lam=%-26s beta=%s  C=%s  sizes=%s  eps=%+d"
          % (lam, beta, C, sorted(len(v) for v in cls.values()), e))

# ---- C3: necessity of concentricity in the wider range
bad = [z for z in S3['zeros'] if z[3] is None or z[4] is None]
print("")
print("  C3  zeros that are NOT concentric      : %d  (must be 0)" % len(bad))

# ---- C2: is eps doing any work?
cnv = sum(n for (e, v), n in S3['tab'].items() if v is False)
print("  C2  concentric shapes that do NOT vanish: %d  (must be > 0, else eps is idle)" % cnv)

# ---- C4: the zeros re-verified, 6 more points and a second prime
N3 = S3['N']
F1, R1 = field_roots(P1)
extra = make_points(F1, R1, 3, N3, 9999, 6)
F2, R2 = field_roots(P2)
pts2 = make_points(F2, R2, 3, N3, 31337, 3)
surv1 = surv2 = 0
for lam, beta, cls, C, e in S3['zeros']:
    if all(matrix(F1, [[al[b] ** beta[a] for b in range(N3)] for a in range(N3)]).det() == 0
           for al in extra):
        surv1 += 1
    if all(matrix(F2, [[al[b] ** beta[a] for b in range(N3)] for a in range(N3)]).det() == 0
           for al in pts2):
        surv2 += 1
print("  C4  zeros surviving 6 further points   : %d of %d" % (surv1, len(S3['zeros'])))
print("      zeros surviving the second prime   : %d of %d" % (surv2, len(S3['zeros'])))

# ---- C5 is already in the table as eps = 0(nc)

# ----------------------------------------------------------------- the alternation, r = 1 .. 4
print("")
print("=" * 88)
print("THE ALTERNATION, one row per r,  t = 4")
print("=" * 88)
print("")
print("   r   N  |lam|<=  shapes  concentric  ZEROS  explained  unexplained   eps that cancels")
print("  " + "-" * 84)
rows = [S3]
for r, MAX, seed in ((1, 30, 7001), (2, 26, 7002), (4, 18, 7004)):
    BRUTE = 14
    rows.append(sweep(r, MAX, seed, 14, "r=%d" % r))
for S in sorted(rows, key=lambda d: d['r']):
    print("  %2d %3d %7d %8d %11d %6d %10d %12d %10s%d"
          % (S['r'], S['N'], S['MAX'], S['nsh'], S['ncon'], len(S['zeros']),
             len(S['expl']), len(S['unex']), "", S['want']))
print("")
print("  'explained' = the reflection cancels the sum TERM BY TERM.  The prediction, written")
print("  before the wide sweep: the unexplained column is ZERO for r odd and large for r even.")
print("")
for S in sorted(rows, key=lambda d: d['r']):
    if S['r'] == 3:
        continue
    if S['unex']:
        print("  r=%d, three unexplained zeros by hand: %s"
              % (S['r'], "; ".join("lam=%s" % z[0] for z in S['unex'][:3])))
print("")
print("DONE")
