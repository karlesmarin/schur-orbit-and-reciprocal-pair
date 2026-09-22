# -*- coding: utf-8 -*-
u"""P3b_pofs_vs_porc.py -- las dos medidas del bloque medido de obs:dosfallos.

obs:dosfallos separa dos fallos de congruencia de la clasificacion:

  n = 3, k = 6  ->  C_6(3,p) toma dos valores polinomicos, p-1 y 0, y la particion de los
                    primos es por escision completa de 2Y^3 + 3Y - 3 (grupo S_3): no PORC,
                    pero POFS (Stanojkovski-Voll, Math. Ann. 381 (2021) 593-629, §1.5.2).
  n >= 4        ->  sobrevive una traza de Frobenius de genero positivo.

Que el segundo NO sea POFS no se mide aqui: lo da el peso uno (thm:clasificacion), y la nota
lo dice asi.  Este guion solo mide lo que el bloque medido afirma:

  1. n = 5, k = 3: el numero de trazas a_p distintas de la curva de conductor 20
     Y^2 = X^3 + X^2 - X, en las ventanas p <= 200, 500, 1000, 2000 (p no | 20), es
     15, 24, 31, 52 (14, 23, 30, 51 sin p = 3, cuyo a_3 = -2 no lo repite ningun otro
     primo); y chi_p(-3) toma dos valores en la misma ventana.  Cuatro ventanas no
     fijan ningun exponente: el crecimiento se imprime, no se convierte en conclusion.
  2. n = 3, k = 6: C_6(3,p), contado por fuerza bruta en todo primo p <= 2000 con p no | 66,
     toma exactamente los dos valores p-1 y 0, y vale p-1 exactamente cuando 2Y^3 + 3Y - 3
     escinde del todo mod p (prop:n3k6); ningun modulo de la lista
     4, 6, 8, 12, 24, 33, 36, 44, 66, 72, 99, 132, 264, 396, 1188 asigna un solo valor a cada
     clase de residuo.  Control negativo: la misma prueba de modulos SI encuentra un modulo
     que separa una cantidad PORC de verdad ([chi_p(-3) = 1], separada por 6).

Uso:  python P3b_pofs_vs_porc.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

import numpy as np

OK = 0
MAL = 0
LINEAS = []


def comprueba(etiqueta, cond, detalle=u''):
    global OK, MAL
    if cond:
        OK += 1
        LINEAS.append(u'  ok   %s%s' % (etiqueta, (u'   ' + detalle) if detalle else u''))
    else:
        MAL += 1
        LINEAS.append(u'  MAL  %s%s' % (etiqueta, (u'   ' + detalle) if detalle else u''))


def seccion(t):
    LINEAS.append(u'')
    LINEAS.append(t)
    LINEAS.append(u'-' * len(t))


def primos(lo, hi):
    crib = [True] * (hi + 1)
    crib[0] = crib[1] = False
    for i in range(2, int(hi ** 0.5) + 1):
        if crib[i]:
            for j in range(i * i, hi + 1, i):
                crib[j] = False
    return [n for n in range(max(2, lo), hi + 1) if crib[n]]


def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def a_p(p):
    u"""Traza de Frobenius de E: Y^2 = X^3 + X^2 - X (conductor 20), por conteo directo."""
    n = 0
    for x in range(p):
        r = (x * x * x + x * x - x) % p
        if r == 0:
            n += 1
        elif pow(r, (p - 1) // 2, p) == 1:
            n += 2
    # #E(F_p) = n + 1 (punto del infinito) = p + 1 - a_p
    return p + 1 - (n + 1)


def C6(p):
    u"""#{T de F_p, |T| = 3, suma 0, suma de sextas 0}, por fuerza bruta (vectorizada)."""
    A = np.arange(p, dtype=np.int64)
    pw = np.array([pow(int(t), 6, p) for t in range(p)], dtype=np.int64)
    a = A[:, None]
    b = A[None, :]
    c = (-a - b) % p
    s = (pw[:, None] + pw[None, :] + pw[c]) % p
    distintos = (a != b) & (c != a) & (c != b)
    n = int(np.count_nonzero((s == 0) & distintos))
    assert n % 6 == 0
    return n // 6


def escinde(p):
    return len([y for y in range(p) if (2 * y ** 3 + 3 * y - 3) % p == 0]) == 3


def modulos_que_separan(ps, clase_de, modulos):
    u"""Los N de la lista para los que toda clase de residuo mod N lleva un solo valor."""
    sep = []
    for N in modulos:
        clases = {}
        mezclada = False
        for p in ps:
            c = p % N
            v = clase_de[p]
            if c in clases and clases[c] != v:
                mezclada = True
                break
            clases[c] = v
        if not mezclada:
            sep.append(N)
    return sep


# ---------------------------------------------------------------------------
seccion(u'1.  n = 5, k = 3: trazas de Frobenius de la curva de conductor 20')
# ---------------------------------------------------------------------------
cuentas, cuentas7 = [], []
for X in (200, 500, 1000, 2000):
    ps = [p for p in primos(3, X) if 20 % p]
    vals = set(a_p(p) for p in ps)
    vals7 = set(a_p(p) for p in ps if p >= 7)
    cuentas.append(len(vals))
    cuentas7.append(len(vals7))
    LINEAS.append(u'       p <= %-5d  p no | 20: %-4d primos, %-3d valores de a_p;'
                  u'  sin p = 3: %d valores' % (X, len(ps), len(vals), len(vals7)))
LINEAS.append(u'   a_3 = %d: con torsion Z/6, a_p = p + 1 mod 6 para p >= 7, asi que a_3 es un valor'
              % a_p(3))
LINEAS.append(u'   que ningun otro primo repite, y la ventana con o sin p = 3 difiere en uno.')
buenos = [p for p in primos(7, 2000)]
comprueba(u'#E(F_p) = p + 1 - a_p es multiplo de 6 en todo 7 <= p <= 2000',
          all((p + 1 - a_p(p)) % 6 == 0 for p in buenos))
comprueba(u'a_3 no lo repite ningun 7 <= p <= 2000', a_p(3) not in set(a_p(p) for p in buenos))
comprueba(u'valores distintos de a_p, p no | 20, en p <= 200, 500, 1000, 2000: 15, 24, 31, 52',
          cuentas == [15, 24, 31, 52], u'%s' % (cuentas,))
comprueba(u'lo mismo con 7 <= p (sin p = 3): 14, 23, 30, 51',
          cuentas7 == [14, 23, 30, 51], u'%s' % (cuentas7,))
comprueba(u'el soporte observado crece en cada ventana',
          all(x < y for x, y in zip(cuentas, cuentas[1:])))
chis = set(chi(-3, p) for p in primos(5, 2000))
comprueba(u'control: chi_p(-3) toma exactamente 2 valores en la misma ventana',
          chis == {1, -1}, u'%s' % (sorted(chis),))
LINEAS.append(u'   (observacion, no comprobacion) cuatro ventanas no fijan ningun exponente;')
LINEAS.append(u'   que este fallo no sea POFS lo da el peso uno, no este crecimiento.')

# ---------------------------------------------------------------------------
seccion(u'2.  n = 3, k = 6: dos valores, y ninguna congruencia los separa')
# ---------------------------------------------------------------------------
ps = [p for p in primos(5, 2000) if 66 % p]
LINEAS.append(u'   primos p <= 2000 con p no | 66: %d  (de %d a %d)' % (len(ps), ps[0], ps[-1]))
C = dict((p, C6(p)) for p in ps)
valores = set(u'p-1' if C[p] == p - 1 else (u'0' if C[p] == 0 else u'otro') for p in ps)
comprueba(u'C_6(3,p), por fuerza bruta, toma exactamente los dos valores p-1 y 0',
          valores == {u'p-1', u'0'}, u'%s' % (sorted(valores),))
malos = [p for p in ps if C[p] != (p - 1) * (1 if escinde(p) else 0)]
comprueba(u'C_6(3,p) = (p-1) [2Y^3+3Y-3 escinde] en todos ellos (prop:n3k6)', not malos,
          u'fallos: %s' % malos[:5])
nesc = sum(1 for p in ps if C[p] == p - 1)
LINEAS.append(u'   primos con C_6 = p-1: %d de %d' % (nesc, len(ps)))

MODS = (4, 6, 8, 12, 24, 33, 36, 44, 66, 72, 99, 132, 264, 396, 1188)
sep = modulos_que_separan(ps, dict((p, C[p] == p - 1) for p in ps), MODS)
comprueba(u'ningun modulo de la lista asigna un solo valor a cada clase de residuo',
          not sep, u'modulos que separarian: %s' % (sep,))
sep_ctrl = modulos_que_separan(ps, dict((p, chi(-3, p) == 1) for p in ps), MODS)
comprueba(u'control negativo: la misma prueba separa [chi_p(-3) = 1], que es PORC',
          6 in sep_ctrl, u'modulos que la separan: %s' % (sep_ctrl,))
LINEAS.append(u'   (observacion, no comprobacion) la particion es por escision completa en un')
LINEAS.append(u'   cuerpo de grupo S_3, que es un conjunto de Frobenius: POFS por definicion')
LINEAS.append(u'   (Stanojkovski-Voll §1.5.2).')

print(u'\n'.join(LINEAS))
print(u'')
print(u'=' * 70)
print(u'VERIFICACION POFS vs PORC:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 70)
sys.exit(1 if MAL else 0)
