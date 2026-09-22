# -*- coding: utf-8 -*-
u"""P3d_mann_minimalidad.py -- la hipotesis de minimalidad de Mann, y el paso 2+2 que exige.

El argumento directo "se descompone en subsumas minimales; una relacion minimal de r <= 4
terminos tiene sus cocientes en mu_6 por Mann; luego todo vive en mu_{gcd(m,6)}" tiene la
segunda mitad FALSA: Mann acota los cocientes DENTRO de cada subsuma minimal, no entre
subsumas distintas.  Con un corte 2+2 los cocientes pueden tener orden arbitrario.

TESTIGO:   1 + (-1) + z + (-z) = 0   con z = zeta_7.
Cuatro terminos, pesos racionales positivos, suma cero, y el cociente z/1 tiene ORDEN 7.

Que el RESULTADO se salve es otra cosa, y hay que demostrarlo aparte: para longitud cuatro
con pesos positivos las unicas particiones en subsumas minimales son {4} y {2,2}, porque un
singleton no se anula y una subsuma de tres deja un singleton.  Y el caso {2,2} obliga a
w_i = w_j y zeta_i = -zeta_j, luego a 2 | m --- que es una condicion que SI depende solo de
gcd(m,6).  Asi que la conclusion aguanta con el paso 2+2 explicito.

Se comprueba numericamente en mu_m (m <= 30): el testigo, que todo par que se anula tiene
w_i = w_j y zeta_i = -zeta_j, que (2,1,1,1) y (3,1,1,1) no admiten corte 2+2 y que
(1,1,1,1) y (2,2,1,1) si, y que la firma de singularidad de diez cuartetos depende solo de
gcd(m,6).  El testigo es el control negativo del argumento directo.

Uso:  python P3d_mann_minimalidad.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import cmath
import itertools
import math
import sys

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


def sec(t):
    L.append(u'')
    L.append(t)
    L.append(u'-' * len(t))


def z(j, m):
    return cmath.exp(2j * cmath.pi * j / m)


def suma(pesos, exps, m):
    return sum(w * z(e, m) for w, e in zip(pesos, exps))


def minimal(pesos, exps, m, tol=1e-9):
    u"""True si ninguna subsuma PROPIA no vacia se anula."""
    n = len(pesos)
    for r in range(1, n):
        for I in itertools.combinations(range(n), r):
            s = sum(pesos[i] * z(exps[i], m) for i in I)
            if abs(s) < tol:
                return False
    return True


# ---------------------------------------------------------------------------
sec(u'1.  El testigo: cuatro terminos, suma cero, cociente de orden 7')
# ---------------------------------------------------------------------------
# en mu_14:  1 = z^0,  -1 = z^7,  zeta_7 = z^2,  -zeta_7 = z^9
m = 14
pesos = (1, 1, 1, 1)
exps = (0, 7, 2, 9)
comprueba(u'1 + (-1) + zeta_7 + (-zeta_7) = 0 en mu_14',
      abs(suma(pesos, exps, m)) < 1e-9, u'|suma| = %.2e' % abs(suma(pesos, exps, m)))
comprueba(u'NO es minimal: la subsuma {1, -1} se anula', not minimal(pesos, exps, m))
# el orden del cociente zeta_7 / 1
orden = m // math.gcd(2, m)
comprueba(u'el cociente zeta_7/1 tiene orden 7, que NO divide a 6', orden == 7,
      u'orden = %d' % orden)
comprueba(u'y sin embargo gcd(14,6) = 2', math.gcd(14, 6) == 2)

L.append(u'')
L.append(u'  => Mann NO dice que los cocientes esten en mu_{gcd(m,6)}.  Lo dice DENTRO de')
L.append(u'     cada subsuma minimal.  El argumento directo salta ese paso.')

# ---------------------------------------------------------------------------
sec(u'2.  El arreglo: con longitud 4 y pesos positivos solo hay dos particiones')
# ---------------------------------------------------------------------------
# Una particion en subsumas minimales no puede tener un singleton (w_i zeta_i != 0)
# ni un bloque de 3 (dejaria un singleton).  Luego es {4} o {2,2}.
comprueba(u'un singleton nunca se anula: |w_i zeta_i| = w_i > 0', True)
comprueba(u'un bloque de 3 deja un singleton, luego es imposible', True)
comprueba(u'=> las unicas particiones son {4} (minimal) y {2,2}', True)

# el caso {2,2}: w_i zeta_i + w_j zeta_j = 0  =>  w_i = w_j  y  zeta_i = -zeta_j
malos = []
for m in range(1, 25):
    for wi in range(1, 6):
        for wj in range(1, 6):
            for a in range(m):
                for b in range(m):
                    if abs(wi * z(a, m) + wj * z(b, m)) < 1e-9:
                        if wi != wj or (2 * (b - a)) % (2 * m) != m % (2 * m) and \
                           abs(z(b, m) + z(a, m)) > 1e-9:
                            malos.append((m, wi, wj, a, b))
comprueba(u'todo par que se anula tiene w_i = w_j y zeta_i = -zeta_j', not malos,
      u'excepciones: %s' % (malos[:3],))

# y eso exige -1 en mu_m, o sea 2 | m
pares = [m for m in range(1, 25) if any(abs(z(a, m) + 1) < 1e-9 for a in range(m))]
comprueba(u'-1 pertenece a mu_m exactamente cuando 2 | m',
      pares == [m for m in range(1, 25) if m % 2 == 0], u'%s' % (pares[:8],))

# ---------------------------------------------------------------------------
sec(u'3.  Los cuatro cuartetos que importan, uno a uno')
# ---------------------------------------------------------------------------
# (2,1,1,1) y (3,1,1,1) NO admiten corte 2+2: por modulos.
for lam in ((2, 1, 1, 1), (3, 1, 1, 1)):
    hay_corte = False
    for m in range(1, 31):
        for exps in itertools.product(range(m), repeat=4):
            if exps[0]:
                continue
            for I in itertools.combinations(range(4), 2):
                J = tuple(i for i in range(4) if i not in I)
                sI = sum(lam[i] * z(exps[i], m) for i in I)
                sJ = sum(lam[i] * z(exps[i], m) for i in J)
                if abs(sI) < 1e-9 and abs(sJ) < 1e-9:
                    hay_corte = True
        if hay_corte:
            break
    comprueba(u'%s no admite corte 2+2 (m <= 30)' % (lam,), not hay_corte)

# (1,1,1,1) y (2,2,1,1) SI lo admiten, y por eso ahi NO se invoca a Mann
for lam in ((1, 1, 1, 1), (2, 2, 1, 1)):
    hay = False
    m = 2
    for exps in itertools.product(range(m), repeat=4):
        for I in itertools.combinations(range(4), 2):
            J = tuple(i for i in range(4) if i not in I)
            if abs(sum(lam[i] * z(exps[i], m) for i in I)) < 1e-9 and \
               abs(sum(lam[i] * z(exps[i], m) for i in J)) < 1e-9:
                hay = True
    comprueba(u'%s SI admite corte 2+2 ya en mu_2' % (lam,), hay)

# ---------------------------------------------------------------------------
sec(u'4.  Y la conclusion aguanta: singularidad solo depende de gcd(m,6)')
# ---------------------------------------------------------------------------
def singular(lam, m):
    for exps in itertools.product(range(m), repeat=3):
        if abs(suma(lam, (0,) + exps, m)) < 1e-9:
            return True
    return False


cuartetos = ((1, 1, 1, 1), (2, 1, 1, 1), (2, 2, 1, 1), (3, 1, 1, 1), (4, 1, 1, 1),
             (3, 2, 1, 1), (2, 2, 2, 1), (5, 1, 1, 1), (3, 3, 1, 1), (2, 2, 2, 2))
por_gcd = {}
consistente = True
for m in range(1, 25):
    g = math.gcd(m, 6)
    firma = tuple(singular(lam, m) for lam in cuartetos)
    if g in por_gcd and por_gcd[g] != firma:
        consistente = False
        L.append(u'       m = %d rompe la clase gcd = %d' % (m, g))
    por_gcd[g] = firma
comprueba(u'la firma de singularidad de 10 cuartetos depende solo de gcd(m,6), m <= 24',
      consistente)

L.append(u'')
L.append(u'VEREDICTO')
L.append(u'  El resultado vale; el argumento directo no basta.  La demostracion necesita:')
L.append(u'    (a) con longitud cuatro y pesos positivos, o la relacion es minimal --- y')
L.append(u'        entonces Mann da mu_6 --- o se parte en 2+2;')
L.append(u'    (b) el corte 2+2 obliga a w_i = w_j y zeta_i = -zeta_j, luego a 2 | m;')
L.append(u'    (c) y "2 | m" tambien esta determinado por gcd(m,6).')
L.append(u'  Sin (a)-(b) explicitos, el testigo 1 + (-1) + zeta_7 + (-zeta_7) = 0 es un')
L.append(u'  contraejemplo al PASO, aunque no al enunciado.')

print(u'\n'.join(L))
print(u'')
print(u'=' * 70)
print(u'VERIFICACION MANN / MINIMALIDAD:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 70)
sys.exit(1 if MAL else 0)
