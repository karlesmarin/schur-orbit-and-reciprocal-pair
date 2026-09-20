# -*- coding: utf-8 -*-
r"""momentos_D12.py -- the residual cases of the theorem of moments for n <= 12, enumerated and closed.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Standard library only (polynomials in symbolic coefficients are dictionaries of monomials).

Setting of the proof: A(Y) = 1 + sum_{g in supp A} c_g Y^g of degree D <= e/2, AB = 1 - Y^e,
Gamma_T = < supp A > (e in <supp A>, D in supp A), multiplicity kappa.  Closed without computation:
kappa = 2 or 3, and two minimal generators.  Left: kappa >= 4, at least three minimal generators,
all <= D, D in Gamma_T, and
      e >= 2D,    D(D-1) <= (e-1)(D-kappa),    c(Gamma_T) > e.
(1) enumerate every such triple (D, Gamma_T, e) with D <= 12;
(2) for each, with symbolic c_g on the elements of Gamma_T in [1, D] and B = [A^{-1}]_{<= e-D},
    find a coefficient of AB in degree 0 < k < e that is a single nonzero monomial in the
    coefficients of the minimal generators: it cannot vanish, contradicting AB = 1 - Y^e
    (the integer factor is 1, 2 or 3 and p > e >= 23).
Public domain (CC0).
"""
import sys
from itertools import combinations
from math import gcd, ceil
from functools import reduce

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def miembros(gens, top):
    S = [False] * (top + 1)
    S[0] = True
    for x in range(1, top + 1):
        S[x] = any(x >= g and S[x - g] for g in gens)
    return S


def conductor(gens):
    top = max(gens) ** 2 + 5
    S = miembros(gens, top)
    c = top
    while c > 0 and S[c - 1]:
        c -= 1
    return c


def minimales(gens):
    gens = sorted(set(gens))
    S = miembros(gens, gens[-1])
    return [a for a in gens if not any(S[x] and S[a - x] for x in range(1, a))]


# (1) the residual triples
casos, vistos = [], set()
for r in range(1, 12):
    for gens in combinations(range(2, 13), r):
        if reduce(gcd, gens) != 1:
            continue
        mins = tuple(minimales(gens))
        if mins in vistos:
            continue
        vistos.add(mins)
        kappa = mins[0]
        if kappa < 4 or len(mins) < 3:
            continue
        c = conductor(list(mins))
        for D in range(max(mins), 13):
            if not miembros(list(mins), D)[D]:
                continue
            emin = max(2 * D, 1 + ceil(D * (D - 1) / (D - kappa)))
            for e in range(emin, c):
                casos.append((D, mins, c, e))
print("(1) residual triples (D, Gamma_T, e) with D <= 12: %d   [shown as (D, Gamma_T, c(Gamma_T), e)]" % len(casos))
for x in casos:
    print("    ", x)
print("     semigroups involved:", sorted({x[1] for x in casos}))


# (2) polynomial arithmetic in the symbols c_g
def pmul(P, Q):
    R = {}
    for m1, a in P.items():
        for m2, b in Q.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            R[m] = R.get(m, 0) + a * b
    return {m: v for m, v in R.items() if v}


def padd(P, Q, s=1):
    R = dict(P)
    for m, v in Q.items():
        R[m] = R.get(m, 0) + s * v
    return {m: v for m, v in R.items() if v}


def texto(P, soporte):
    (m, a), = P.items()
    f = "*".join("c%d%s" % (g, "^%d" % x if x > 1 else "") for g, x in zip(soporte, m) if x)
    return ("%d*" % a if a != 1 else "") + f if a not in (1,) else f


todos = True
for D, mins, c, e in casos:
    S = miembros(list(mins), D)
    soporte = [g for g in range(1, D + 1) if S[g]]
    nv = len(soporte)
    cero = tuple([0] * nv)
    var = {g: {tuple(int(i == j) for i in range(nv)): 1} for j, g in enumerate(soporte)}
    A = {0: {cero: 1}}
    for g in soporte:
        A[g] = var[g]
    B = {0: {cero: 1}}
    for k in range(1, e - D + 1):
        acc = {}
        for g in soporte:
            if g <= k:
                acc = padd(acc, pmul(A[g], B[k - g]), -1)
        B[k] = acc
    idx_min = {soporte.index(g) for g in mins}
    hallado = None
    for k in range(1, e):
        acc = {}
        for g, a in A.items():
            if 0 <= k - g <= e - D:
                acc = padd(acc, pmul(a, B[k - g]))
        if len(acc) == 1:
            (m, v), = acc.items()
            if all(x == 0 or i in idx_min for i, x in enumerate(m)) and abs(v) in (1, 2, 3):
                hallado = (k, texto(acc, soporte))
                break
    todos &= hallado is not None
    print("(2) D=%2d Gamma=%s e=%d: %s" % (D, mins, e, "[Y^%d] of AB = %s" % hallado if hallado
                                           else "*** NOT EXCLUDED ***"))
print("every residual triple excluded:", todos, "| smallest e:", min(x[3] for x in casos))
