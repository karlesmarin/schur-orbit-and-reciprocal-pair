# -*- coding: utf-8 -*-
# The Laplace parity, rewritten from scratch, with the printed case as an acceptance test.
#
# The previous attempt (quadrant_parity.sage, retracted) enumerated the terms and reported numbers
# that contradicted a case computed by hand.  Rather than debug it, this recomputes eps by a SECOND
# and independent route and demands the two agree on every shape, so that the failure mode that
# caught me -- one method quietly wrong, nothing to compare it against -- cannot repeat.
#
# Route A: enumerate the transversals, take eps = sgn(C-T) * sgn(T) per term.
#
# Route B: from the structure.  Write E for the excess positions in increasing order, |E| = m, and
# S for the singleton positions.  The reflection b -> C-b REVERSES E, because beta decreases with
# the index, and fixes S pointwise.  The Laplace sign is (-1)^inv with inv counting pairs (i in
# pick, j in T) with i > j, and it splits:
#   * the part with i in E \ T flips to its complement, changing by (-1)^{(m-2r)*2r} = +1, since
#     2r is even -- so the whole E-internal contribution cancels, for every r;
#   * the part with i in S changes by Delta = sum over singletons of
#         #(T meets the top a_i of E)  -  #(T meets the bottom a_i of E),
#     where a_i is the number of excess positions below the singleton i.
# So eps = (-1)^Delta, and its being constant is the statement that Delta has constant parity.
#
# ACCEPTANCE TEST, run first and fatal: lambda = (5,4,3), beta = [12,10,8,4,3,2,1,0], C = 12 must
# give 8 terms, every pair of opposite sign, eps = -1, and a signed total of exactly zero.  Those
# numbers are on the record in quadrant_debug_OUT.txt.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

t, r = 4, 2
N = t + 2 * r
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
zs = list(R.gens())


def beta_of(lam):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def split(beta):
    cls = {}
    for k, b in enumerate(beta):
        cls.setdefault(b % t, []).append(k)
    E = sorted(i for v in cls.values() if len(v) >= 2 for i in v)
    S = sorted(i for v in cls.values() if len(v) == 1 for i in v)
    return cls, E, S


def inv_sign(pick, T):
    return (-1) ** sum(1 for i in pick for j in T if i > j)


def alt(vals):
    return matrix(R, len(vals), len(vals),
                  lambda a, b: (zs[b // 2] ** vals[a]) if b % 2 == 0
                  else (zs[b // 2] ** (-vals[a]))).det()


def route_A(beta, cls, C, want_total=False):
    """eps per term; returns (eps or 0 if not constant, signed total or None)."""
    pos = dict((b, i) for i, b in enumerate(beta))
    eps, total = None, R(0)
    for pick in itertools.product(*[cls[k] for k in sorted(cls)]):
        pick = tuple(sorted(pick))
        T = tuple(sorted(set(range(N)) - set(pick)))
        s1 = inv_sign(pick, T)
        if want_total:
            total += s1 * alt([beta[i] for i in T])
        try:
            T2 = tuple(sorted(pos[C - beta[i]] for i in T))
        except KeyError:
            return None, None
        p2 = tuple(sorted(set(range(N)) - set(T2)))
        e = s1 * inv_sign(p2, T2)
        if eps is None:
            eps = e
        elif eps != e:
            return 0, (total if want_total else None)
    return eps, (total if want_total else None)


def route_B(beta, cls, C):
    """eps from the structure, and 0 if Delta does not have constant parity."""
    _, E, S = split(beta)
    m = len(E)
    pos = dict((b, i) for i, b in enumerate(beta))
    if any(C - beta[i] not in pos for i in E):
        return None
    seen = set()
    for pick in itertools.product(*[cls[k] for k in sorted(cls)]):
        T = sorted(set(range(N)) - set(pick))
        d = 0
        for i in S:
            a = sum(1 for j in E if j < i)
            top, bot = set(E[m - a:]), set(E[:a])
            d += len(set(T) & top) - len(set(T) & bot)
        seen.add(d % 2)
    if len(seen) > 1:
        return 0
    return (-1) ** seen.pop()


print("=" * 78)
print("ACCEPTANCE TEST on the case printed in quadrant_debug_OUT.txt")
print("=" * 78)
lam0 = [5, 4, 3]
b0 = beta_of(lam0)
cls0, E0, S0 = split(b0)
epsA, tot = route_A(b0, cls0, 12, want_total=True)
epsB = route_B(b0, cls0, 12)
nterms = 1
for v in cls0.values():
    nterms *= len(v)
print("  beta            : %s" % b0)
print("  terms           : %d   (expected 8)" % nterms)
print("  eps by route A  : %s   (expected -1)" % epsA)
print("  eps by route B  : %s   (expected -1)" % epsB)
print("  signed total    : %s   (expected zero)" % ("zero" if tot == 0 else "NOT ZERO"))
ok = (nterms == 8 and epsA == -1 and epsB == -1 and tot == 0)
print("  ACCEPTANCE: %s" % ("PASS" if ok else "FAIL -- everything below is void"))
if not ok:
    raise SystemExit(1)

print("")
print("=" * 78)
print("eps over the reflective shapes, both routes,  |lambda| <= 22")
print("=" * 78)

agree = disagree = 0
rows = []
for size in range(23):
    for l in Partitions(size, max_length=N):
        beta = beta_of(list(l))
        cls, E, S = split(beta)
        if len(cls) < t or not E:
            continue
        vals = set(beta[i] for i in E)
        Cs = [C for C in range(2 * max(vals) + 1) if set(C - b for b in vals) == vals]
        for C in Cs:
            a = route_A(beta, cls, C)[0]
            b = route_B(beta, cls, C)
            if a is None or b is None:
                continue
            if a == b:
                agree += 1
            else:
                disagree += 1
                if disagree <= 5:
                    print("    DISAGREE lam=%s C=%d  A=%s B=%s" % (list(l), C, a, b))
            rows.append(dict(lam=list(l), beta=beta, C=C, eps=a, m=len(E), S=S, E=E,
                             sizes=tuple(sorted(len(v) for v in cls.values()))))
            break

print("")
print("  CONTROL, the two routes agree : %d   disagree: %d" % (agree, disagree))
print("  shapes measured               : %d" % len(rows))
print("  eps = -1 (the sum cancels)    : %d" % sum(1 for d in rows if d['eps'] == -1))
print("  eps = 0  (not constant)       : %d" % sum(1 for d in rows if d['eps'] == 0))

print("")
print("  eps against the singleton geometry, which is what route B says decides it:")
tab = {}
for d in rows:
    a_list = tuple(sorted(sum(1 for j in d['E'] if j < i) for i in d['S']))
    tab.setdefault((d['m'], a_list), set()).add(d['eps'])
for key in sorted(tab):
    print("    m=%d  singletons at depths %-14s ->  eps in %s"
          % (key[0], str(key[1]), sorted(tab[key])))
print("  cells carrying BOTH signs: %d of %d"
      % (sum(1 for v in tab.values() if len(v) > 1), len(tab)))

print("")
print("DONE")
