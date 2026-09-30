# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# COMPUERTA (28-sep-2026): el menor de Kucera frente al mcd de los menores maximales de M_m.
#
# QUE: para cada m <= 90, m no = 2 (mod 4), m >= 3, se construye M_m = [s_m(u^-1 c)] (filas u en el
#   semisistema Y de Kucera, columnas TODOS los c mod m; s_m(0)=0, s_m(c)=2c-m) y se compara:
#     - el menor d x d de M_m sobre las columnas M_- de Kucera (d = phi(m)/2),
#     - Delta = mcd de los menores maximales = indice del reticulo de columnas en Z^d.
# POR QUE: la auditoria de prior art (28-sep-2026) encontro que Kucera, J. Number Theory 40 (1992)
#   284-316, Teoremas 5.1 (p.301) y 5.2 (p.313), evalua |det(1/2 - <tk/m>)_{k en M_-, t en Y}| =
#   2^b h^- / w (b = 0 si m es potencia de primo, b = 2^{n-2}-1 si tiene n >= 2 primos), que por el
#   factor (-2m)^d es un menor maximal de M_m con el valor de la Proposicion prop:sinnott.  Esta
#   compuerta mide que ese menor coincide con Delta (lo que Kucera NO enuncia: que sea el mcd).
# FUENTES: M_+- y N_+- de Kuc92 p.293; Y = {a: (a,m)=1, <a/q_n> < 1/2} de Kuc92 Lemma 4.2 (p.298),
#   q_n = potencia del mayor primo (primos en orden creciente).  Lectura: 0 esta en N_- (se toma
#   "a | m" como cierto para a = 0, que es m mod m); con esa lectura card M_- = phi(m)/2 en todo el rango.
# CONTROLES que pueden fallar:
#   (1) card M_- = d; (2) menor(M_-) != 0 y menor(M_-) == Delta; (3) todo menor calculado es multiplo
#   de Delta; (4) control negativo: M_- con su mayor elemento cambiado por el menor residuo fuera de
#   M_- u {0}; se cuenta cuantas veces su menor != Delta (si nunca difiere, la compuerta no discrimina);
#   (5) el bloque de unidades (columnas = Y): donde es singular su menor es 0 != Delta; donde no, se
#   imprime cociente; (6) para los 17 m de la auditoria, Delta tambien por SNF de sympy.
#
# Authors: Carles Marin, Claude (AI assistant).
from fractions import Fraction as F
from math import gcd
from sympy import factorint, Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

AUDITORIA = [5, 7, 9, 12, 15, 16, 20, 21, 23, 24, 28, 35, 39, 40, 56, 60, 84]


def frac(x):
    return x - (x.numerator // x.denominator)


def m_menos(m):
    fa = sorted(factorint(m).items())
    ps = [p for p, _ in fa]
    qs = [p ** e for p, e in fa]
    n = len(ps)
    X = [a for a in range(m) if all((a % p != 0) or (a % q == 0) for p, q in zip(ps, qs))]
    X0 = sorted(set(X) | {0})
    Xs = set(X)
    cnd = lambda a: sum(1 for q in qs if a % q != 0)
    out = []
    for a in X0:
        g = gcd(m, a)
        c1 = a in Xs and any(a % q != 0 and ((a // g) + 1) % q == 0 for q in qs)
        nm = (a == 0 or m % a == 0) and (-1) != (-1) ** cnd(a)
        c3 = any(a in Xs and a % qs[k] != 0 and frac(F(a, qs[k] * g)) > F(1, 2)
                 and all((a - g) % qs[i] == 0 for i in range(k + 1, n)) for k in range(n))
        if not (c1 or nm or c3):
            out.append(a)
    return out, qs


def s(m, c):
    c %= m
    return 0 if c == 0 else 2 * c - m


def det_bareiss(A):
    A = [row[:] for row in A]
    k = len(A)
    sign, prev = 1, 1
    for i in range(k - 1):
        if A[i][i] == 0:
            sw = next((r for r in range(i + 1, k) if A[r][i] != 0), None)
            if sw is None:
                return 0
            A[i], A[sw] = A[sw], A[i]
            sign = -sign
        for r in range(i + 1, k):
            for c in range(i + 1, k):
                A[r][c] = (A[r][c] * A[i][i] - A[r][i] * A[i][c]) // prev
        prev = A[i][i]
    return sign * A[k - 1][k - 1]


def egcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def indice_reticulo(vecs, d, D):
    # indice en Z^d del reticulo generado por vecs; D != 0 es un multiplo conocido del indice,
    # asi que el reticulo contiene D.Z^d y las entradas fuera del pivote se reducen mod D.
    base = [None] * d
    todos = [list(v) for v in vecs] + [[D if j == i else 0 for j in range(d)] for i in range(d)]
    for v in todos:
        for i in range(d):
            if v[i] == 0:
                continue
            if base[i] is None:
                base[i] = [v[j] if j <= i else v[j] % D for j in range(d)]
                break
            b = base[i]
            g, x, y = egcd(b[i], v[i])
            nb = [x * bb + y * vv for bb, vv in zip(b, v)]
            nv = [(v[i] // g) * bb - (b[i] // g) * vv for bb, vv in zip(b, v)]
            base[i] = [nb[j] if j <= i else nb[j] % D for j in range(d)]
            v = [nv[j] if j <= i else nv[j] % D for j in range(d)]
    if any(b is None for b in base):
        return 0
    r = 1
    for i in range(d):
        r *= base[i][i]
    return abs(r)


def phi(m):
    return sum(1 for a in range(1, m + 1) if gcd(a, m) == 1)


fallos = 0
neg_distinto = 0
unid_sing = 0
casos = 0
print("m  d  card(M-)  menor(M-)  Delta  cociente  unidades  control_neg(menor==Delta?)")
for m in range(3, 91):
    if m % 4 == 2:
        continue
    casos += 1
    Mm, qs = m_menos(m)
    qn = qs[-1]
    Y = [t for t in range(1, m) if gcd(t, m) == 1 and frac(F(t, qn)) < F(1, 2)]
    d = len(Y)
    ok = True
    if d != phi(m) // 2 or len(Mm) != d:
        ok = False
    inv = [pow(u, -1, m) for u in Y]
    col = lambda c: [s(m, ui * c) for ui in inv]
    sub = lambda cs: [[s(m, ui * c) for c in cs] for ui in inv]
    menor = abs(det_bareiss(sub(Mm))) if len(Mm) == d else 0
    if menor == 0:
        ok = False
        D = abs(det_bareiss(sub(Y))) or 1
        Delta = indice_reticulo([col(c) for c in range(m)], d, D) if D != 1 else None
    else:
        Delta = indice_reticulo([col(c) for c in range(m)], d, menor)
    if not Delta or menor != Delta:
        ok = False
    # (5) bloque de unidades
    du = abs(det_bareiss(sub(Y)))
    if du == 0:
        unid_sing += 1
        txt_u = "singular"
    else:
        if Delta and du % Delta != 0:
            ok = False
        txt_u = "no sing, det/Delta=" + str(F(du, Delta) if Delta else "?")
    # (4) control negativo
    fuera = [c for c in range(1, m) if c not in Mm]
    alt = sorted(Mm[:-1] + [fuera[0]]) if fuera else Mm
    dalt = abs(det_bareiss(sub(alt)))
    if Delta and dalt % Delta != 0:
        ok = False
    if dalt != Delta:
        neg_distinto += 1
    # (6) SNF de sympy en los 17 m de la auditoria
    txt_s = ""
    if m in AUDITORIA:
        snf = smith_normal_form(Matrix([col(c) for c in range(m)]).T, domain=ZZ)
        P = 1
        for i in range(d):
            P *= snf[i, i]
        if abs(P) != Delta:
            ok = False
        txt_s = "  [SNF sympy = Delta: " + str(abs(P) == Delta) + "]"
    if not ok:
        fallos += 1
    print(m, d, len(Mm), menor, Delta, (F(menor, Delta) if Delta else "?"), txt_u,
          "alt=" + str(alt) + " igual=" + str(dalt == Delta), "OK" if ok else "FALLO", txt_s)
print()
print("casos:", casos, " fallos:", fallos)
print("bloque de unidades singular en", unid_sing, "de", casos, "(ahi su menor 0 != Delta)")
print("control negativo: menor != Delta en", neg_distinto, "de", casos,
      "(si fuese 0, la compuerta no discriminaria)")
