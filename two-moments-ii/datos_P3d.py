# -*- coding: utf-8 -*-
u"""Produce los datos que dibujan las figuras de P3d, y los deja en disco.

REGLA: la figura no calcula, lee.  Si el dato se recalcula dentro del guion de dibujo, no
hay forma de auditarlo por separado ni de saber si la figura que esta en el PDF corresponde
al numero que dice la prosa.

SALIDA:  datos_P3d.json  con
   grado[k]      deg Q_k, medido factorizando s_k
   Qk[k]         los coeficientes de Q_k en u, de mayor a menor
   valores[k]    { valor de 6*C_k(3,p)/(p-1) : lista de primos que lo dan }, p <= PMAX

    python datos_P3d.py
"""
from __future__ import print_function

import json
import os

import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
KMAX_GRADO = 60
KMAX_VAL = 14
PMAX = 257

e2, e3, u = sp.symbols('e2 e3 u')


def primos(lo, hi):
    crib = [True] * (hi + 1)
    crib[0] = crib[1] = False
    for i in range(2, int(hi ** 0.5) + 1):
        if crib[i]:
            for j in range(i * i, hi + 1, i):
                crib[j] = False
    return [n for n in range(max(2, lo), hi + 1) if crib[n]]


# ---------------------------------------------------------------- s_k y Q_k
S = {0: sp.Integer(3), 1: sp.Integer(0), 2: -2 * e2, 3: 3 * e3}
for k in range(4, KMAX_GRADO + 1):
    S[k] = sp.expand(-e2 * S[k - 2] + e3 * S[k - 3])


def Qk(k):
    u"""Coeficientes de Q_k en u (mayor a menor); [1] si es constante."""
    resto = sp.Integer(1)
    for base, m in sp.factor_list(sp.factor(S[k]))[1]:
        b = sp.expand(base)
        if b not in (e2, e3):
            resto *= b ** m
    if resto == 1:
        return [1]
    P = sp.Poly(sp.expand(resto), e2, e3)
    expr = sp.Integer(0)
    for (a, b), c in zip(P.monoms(), P.coeffs()):
        assert a % 3 == 0 and b % 2 == 0, u'k=%d no homogeneo en u' % k
        expr += c * u ** (a // 3)
    return [int(c) for c in sp.Poly(sp.expand(expr), u).all_coeffs()]


# ---------------------------------------------------------------- el recuento
def cuenta(k, p):
    u"""#{T de 3 elementos distintos de F_p con suma 0 y suma de k-esimas 0}, en O(p^2)."""
    c = 0
    pk = [pow(x, k, p) for x in range(p)]
    for a in range(p):
        for b in range(a + 1, p):
            d = (-a - b) % p
            if d <= b:
                continue
            if (pk[a] + pk[b] + pk[d]) % p == 0:
                c += 1
    return c


if __name__ == '__main__':
    print(u'grados y Q_k hasta k = %d' % KMAX_GRADO)
    grado, coefs = {}, {}
    for k in range(3, KMAX_GRADO + 1):
        q = Qk(k)
        coefs[k] = q
        grado[k] = len(q) - 1
    print(u'   ceros (Q_k constante): %s'
          % [k for k in sorted(grado) if grado[k] == 0])

    print(u'valores medidos de 6C/(p-1) hasta k = %d, p <= %d' % (KMAX_VAL, PMAX))
    PS = primos(5, PMAX)
    valores = {}
    for k in range(3, KMAX_VAL + 1):
        d = {}
        for p in PS:
            if k % p == 0:
                continue
            r = sp.Rational(cuenta(k, p) * 6, p - 1)
            d.setdefault(str(r), []).append(p)
        valores[k] = d
        print(u'   k = %-2d  %s' % (k, sorted(d, key=lambda s: sp.Rational(s))))

    ruta = os.path.join(AQUI, u'datos_P3d.json')
    with open(ruta, 'w') as f:
        json.dump({u'grado': {str(k): v for k, v in grado.items()},
                   u'Qk': {str(k): v for k, v in coefs.items()},
                   u'valores': {str(k): v for k, v in valores.items()},
                   u'PMAX': PMAX, u'n_primos': len(PS)}, f, indent=1)
    print(u'escrito %s' % ruta)
