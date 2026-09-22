# -*- coding: utf-8 -*-
u"""P3d_grado_Qk.py -- formula cerrada para el grado de Q_k.

Con Q_k el factor propio de s_k escrito en u = e2^3/e3^2,

        deg Q_k  =  floor(k/6) - [ k = 1 (mod 6) ]

Las seis clases de k modulo 6 son seis rectas de pendiente 1/6, y la de k = 1 va un escalon
por debajo.  Eso explica el k = 16: 16 = 4 (mod 6), no tiene el escalon, luego
floor(16/6) = 2, y no 3 como daria una extrapolacion lineal ingenua.

El grado se calcula de dos maneras independientes --- factorizando s_k simbolicamente
(sympy) y contando, y por la formula --- y se comparan para 3 <= k <= 90.  La primera
comprobacion (todo factor propio es homogeneo en u) es la que puede fallar si la
reduccion a u no fuera legitima; k = 16 y la escalera k = 7, 13, 19, 25 son los casos que
separan la formula de sus vecinas sin escalon.

Uso:  python P3d_grado_Qk.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

import sympy as sp

e2, e3, u = sp.symbols('e2 e3 u')

OK = 0
MAL = 0
L = []


def comprueba(et, cond, det=u''):
    global OK, MAL
    if cond:
        OK += 1
        L.append(u'  ok   %s%s' % (et, (u'   ' + det) if det else u''))
    else:
        MAL += 1
        L.append(u'  MAL  %s%s' % (et, (u'   ' + det) if det else u''))


KMAX = 90
S = {0: sp.Integer(3), 1: sp.Integer(0), 2: -2 * e2, 3: 3 * e3}
for k in range(4, KMAX + 1):
    S[k] = sp.expand(-e2 * S[k - 2] + e3 * S[k - 3])


def grado_medido(k):
    u"""Factoriza s_k, quita las potencias de e2 y e3, y devuelve el grado en u."""
    facs = sp.factor_list(sp.factor(S[k]))[1]
    resto = sp.Integer(1)
    for base, m in facs:
        b = sp.expand(base)
        if b not in (e2, e3):
            resto *= b ** m
    if resto == 1:
        return 0
    P = sp.Poly(sp.expand(resto), e2, e3)
    g = 0
    for (a, b), _ in zip(P.monoms(), P.coeffs()):
        if a % 3 or b % 2:
            return None                       # no homogeneo en u: seria un fallo real
        g = max(g, a // 3)
    return g


def grado_formula(k):
    return k // 6 - (1 if k % 6 == 1 else 0)


L.append(u'  k    deg medido   deg formula')
L.append(u'  ' + u'-' * 34)
malos = []
no_homog = []
for k in range(3, KMAX + 1):
    gm = grado_medido(k)
    if gm is None:
        no_homog.append(k)
        continue
    gf = grado_formula(k)
    if gm != gf:
        malos.append((k, gm, gf))
    if k <= 26:
        L.append(u'  %-4d %-12d %d' % (k, gm, gf))
L.append(u'  ...')

comprueba(u'todo factor propio de s_k es homogeneo en u = e2^3/e3^2', not no_homog,
      u'fallos: %s' % (no_homog,))
comprueba(u'deg Q_k = floor(k/6) - [k = 1 mod 6], para 3 <= k <= %d' % KMAX, not malos,
      u'discrepancias: %s' % (malos[:5],))

# el caso que separa la formula de una extrapolacion lineal
comprueba(u'k = 16: la formula da 2, no 3 (16 = 4 mod 6, sin escalon)',
      grado_formula(16) == 2 and grado_medido(16) == 2)
comprueba(u'k = 7, 13, 19, 25: son los del escalon y dan 0, 1, 2, 3',
      [grado_medido(k) for k in (7, 13, 19, 25)] == [0, 1, 2, 3])

# Q_k constante <=> grado 0 <=> k en {3,4,5,7}
ceros = [k for k in range(3, KMAX + 1) if grado_medido(k) == 0]
comprueba(u'Q_k es constante exactamente para k en {3,4,5,7}', ceros == [3, 4, 5, 7],
      u'%s' % (ceros,))
comprueba(u'y la formula lo predice sin calcular nada',
      [k for k in range(3, KMAX + 1) if grado_formula(k) == 0] == [3, 4, 5, 7])

L.append(u'')
L.append(u'VEREDICTO')
L.append(u'  deg Q_k = floor(k/6) - [k = 1 (mod 6)].  Seis rectas de pendiente 1/6, una de')
L.append(u'  ellas un escalon por debajo.  De ahi salen los cuatro ceros {3,4,5,7} sin')
L.append(u'  factorizar nada, y de ahi sale tambien que NO hay mas: el grado es creciente')
L.append(u'  en cada clase, asi que una vez pasa de cero no vuelve.')
L.append(u'  Queda abierta la recursion en u: no hay Q_k = A(u) Q_{k-a} + B(u) Q_{k-b}')
L.append(u'  con A, B fijos para (a,b) en')
L.append(u'  {(3,6), (2,4), (1,2), (6,12)}.')

print(u'\n'.join(L))
print(u'')
print(u'=' * 70)
print(u'VERIFICACION GRADO DE Q_k:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 70)
sys.exit(1 if MAL else 0)
