# -*- coding: utf-8 -*-
u"""P3d_weil_D.py -- los polinomios de Weil de la curva de genero dos D en p = 11 y p = 19,
desde los recuentos sobre F_p Y sobre F_{p^2}.

D: v^2 = t^6 + 12t^4 + 40t^3 + 36t^2 + 48t + 16 (P3d, la curva de genero cinco de k = 5, n = 6).
Para una curva de genero dos con polinomio de Weil T^4 + a1 T^3 + a2 T^2 + a1 p T + p^2,
    #D(F_p)     = p + 1 + a1,
    #D(F_{p^2}) = p^2 + 1 - (a1^2 - 2 a2),
asi que #D(F_p) solo fija a1, y a2 pide el recuento sobre F_{p^2}.  Los dos recuentos se hacen
por fuerza bruta, F_{p^2} construido como F_p[i]/(i^2 - r) con r no residuo.  El sextico es
monico, asi que en el infinito hay dos puntos racionales sobre todo cuerpo.

Control: con SOLO #D(F_p) hay decenas de valores de a2 compatibles con la cota de Weil (las dos
raices de X^2 + a1 X + (a2 - 2p) reales y de modulo <= 2 sqrt p), de modo que el recuento sobre
F_{p^2} no sobra: es el que fija a2.

Uso:  python P3d_weil_D.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

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


F = [16, 48, 36, 40, 12, 0, 1]          # coeficientes de t^0 .. t^6


def no_residuo(p):
    return next(r for r in range(2, p) if pow(r, (p - 1) // 2, p) == p - 1)


def cuenta_Fp(p):
    n = 2                                 # los dos puntos del infinito
    for t in range(p):
        y = sum(c * pow(t, i, p) for i, c in enumerate(F)) % p
        n += 1 if y == 0 else (2 if pow(y, (p - 1) // 2, p) == 1 else 0)
    return n


def cuenta_Fp2(p):
    r = no_residuo(p)

    def mul(a, b):
        return ((a[0] * b[0] + r * a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    def pot(a, e):
        x = (1, 0)
        while e:
            if e & 1:
                x = mul(x, a)
            a = mul(a, a)
            e >>= 1
        return x

    q = p * p
    n = 2
    for t0 in range(p):
        for t1 in range(p):
            t = (t0, t1)
            y = (0, 0)
            tp = (1, 0)
            for c in F:
                y = ((y[0] + c * tp[0]) % p, (y[1] + c * tp[1]) % p)
                tp = mul(tp, t)
            if y == (0, 0):
                n += 1
            else:
                n += 2 if pot(y, (q - 1) // 2) == (1, 0) else 0
    return n


def howe_zhu(a1, a2, p):
    """las condiciones del Teorema 6 de Howe--Zhu que se cumplen (vacia = absolutamente simple)."""
    c = [(u'a1 = 0', a1 == 0), (u'a1^2 = p + a2', a1 * a1 == p + a2),
         (u'a1^2 = 2 a2', a1 * a1 == 2 * a2), (u'a1^2 = 3 a2 - 3p', a1 * a1 == 3 * a2 - 3 * p)]
    return [n for n, v in c if v]


IMPRESO = {11: (-1, -5), 19: (5, 18)}     # (a1, a2) que imprime la nota

for p in (11, 19):
    N1, N2 = cuenta_Fp(p), cuenta_Fp2(p)
    a1 = N1 - p - 1
    a2 = (N2 - p * p - 1 + a1 * a1) // 2
    comprueba(u'p = %d: #D(F_p) = %d' % (p, N1), N1 == {11: 11, 19: 25}[p])
    comprueba(u'p = %d: #D(F_{p^2}) = %d da a2 entero' % (p, N2),
              (N2 - p * p - 1 + a1 * a1) % 2 == 0)
    comprueba(u'p = %d: (a1, a2) = (%d, %d), el par impreso' % (p, a1, a2),
              (a1, a2) == IMPRESO[p], u'impreso %s' % (IMPRESO[p],))
    # con solo N1: cuantos a2 respeta la cota de Weil (|a2| <= 2p + a1^2/4 basta de sobra)
    candidatos = [b for b in range(-6 * p, 6 * p + 1)
                  if all(abs(x) <= 2 * p ** 0.5 + 1e-9 for x in
                         __import__('numpy').roots([1, a1, b - 2 * p]).real)
                  and all(abs(__import__('numpy').roots([1, a1, b - 2 * p]).imag) < 1e-9)]
    comprueba(u'control: con solo #D(F_%d) caben %d valores de a2 --- F_{p^2} no sobra'
              % (p, len(candidatos)), len(candidatos) > 1)
    L.append(u'       p = %d: N1 = %d, N2 = %d' % (p, N1, N2))
    # Howe--Zhu, Teorema 6: para una superficie simple ORDINARIA (p no divide a2) la simplicidad
    # absoluta falla solo si a1 = 0, a1^2 = p + a2, a1^2 = 2 a2 o a1^2 = 3 a2 - 3p.
    comprueba(u'p = %d: ordinaria, p no divide a2 = %d' % (p, a2), a2 % p != 0)
    comprueba(u'p = %d: ninguna de las cuatro condiciones de Howe--Zhu: %d no esta en {%d, %d, %d}'
              % (p, a1 * a1, p + a2, 2 * a2, 3 * a2 - 3 * p), not howe_zhu(a1, a2, p))
    # control negativo: la restriccion de Weil T^4 - T^2 + 121 en p = 11, (a1, a2) = (0, -1),
    # ordinaria y con polinomio irreducible, tiene que caer en a1 = 0 por la MISMA funcion
    if p == 11:
        comprueba(u'control negativo: (a1, a2) = (0, -1) en p = 11 dispara %s'
                  % (howe_zhu(0, -1, 11),), howe_zhu(0, -1, 11) == [u'a1 = 0'])
    # los resolventes X^2 + a1 X + (a2 - 2p) y sus discriminantes
    d = a1 * a1 - 4 * (a2 - 2 * p)
    comprueba(u'p = %d: resolvente X^2 %+dX %+d, discriminante %d' % (p, a1, a2 - 2 * p, d),
              d == {11: 109, 19: 105}[p])

print(u'\n'.join(L))
print(u'=' * 70)
print(u'VERIFICACION POLINOMIOS DE WEIL DE D:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 70)
sys.exit(1 if MAL else 0)
