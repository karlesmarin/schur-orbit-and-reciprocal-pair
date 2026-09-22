# -*- coding: utf-8 -*-
u"""P3c_familia_a11111.py -- la formula de la familia V_(a,1,1,1,1) necesita hipotesis sobre p; aqui se miden cuales.

La Proposicion 11 de P3c (etiqueta prop:familia) da

    #V_(a,1,1,1,1)(F_p) = p^2 + [3 + 3 chi_p(D_a) + chi_p(E_a)] p + 1,
    D_a = -3(a^2-4),   E_a = (a^2-4)(a^2-16),

y sin condicion sobre p la formula falla: en a = p = 5 el recuento real es 76 y la formula
extrapolada da 31.

LO QUE SE MIDE.  Sobre una rejilla de (a, p) se comprueban dos condiciones suficientes y se
busca si alguna deja pasar un fallo:

    (C1)  p no divide  6 a (a^2-4) (a^2-16)
    (C2)  p > a + 4

y ademas se localiza EXACTAMENTE donde falla la formula sin condicion, para no adjudicar a la
hipotesis mas de lo que hace falta: en la rejilla los fallos son exactamente los (a, p) con
p | a, de modo que C1 es suficiente pero no necesaria.

Y UN HECHO DE GRUPOS que usa la demostracion de la proposicion para pasar de las once rectas
que da el haz (L y las diez que la cortan) a las veintisiete: el estabilizador de una recta en
W(E_6) tiene orden 1920 (es W(D_5)) y actua FIELMENTE sobre las diez rectas que la cortan.  Se
construye W(E_6) como grupo de permutaciones de las 27 rectas (modelo del plano explotado en
seis puntos: a_i, b_j, c_ij) y se comprueba.

CONTROLES:
  - control positivo: a = 1 es la Clebsch, cuyo recuento esta PUBLICADO en la nota;
  - se excluyen a = 2 y a = 4 porque el enunciado ya los excluye (ahi D_a o E_a se anulan);
  - y se mide si C2 implica C1 en el rango, o son independientes.

Uso:  python P3c_familia_a11111.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

OK = MAL = 0
L = []


def comprueba(et, cond, det=u''):
    global OK, MAL
    if cond:
        OK += 1
        L.append(u'  ok   %s%s' % (et, (u'   ' + det) if det else u''))
    else:
        MAL += 1
        L.append(u'  MAL  %s%s' % (et, (u'   ' + det) if det else u''))


def chi(x, p):
    x %= p
    return 0 if x == 0 else (1 if pow(x, (p - 1) // 2, p) == 1 else -1)


def contar(a, p):
    u"""#V_(a,1,1,1,1)(F_p) por fuerza bruta, eliminando la ultima variable de peso 1."""
    w = [a % p, 1, 1, 1, 1]
    cub = (np.arange(p, dtype=np.int64) ** 3) % p
    g = np.arange(p, dtype=np.int64)
    X2 = g.reshape(p, 1, 1)
    X3 = g.reshape(1, p, 1)
    X4 = g.reshape(1, 1, p)
    base = (cub[X2] + cub[X3] + cub[X4]) % p
    lin = (X2 + X3 + X4) % p
    ceros = 0
    for x1 in range(p):
        X5 = (-(w[0] * x1 + lin)) % p
        F = (w[0] * cub[x1] + base + cub[X5]) % p
        ceros += int(np.count_nonzero(F == 0))
    return (ceros - 1) // (p - 1)


def formula(a, p):
    Da = -3 * (a * a - 4)
    Ea = (a * a - 4) * (a * a - 16)
    return p * p + (3 + 3 * chi(Da, p) + chi(Ea, p)) * p + 1


AES = [a for a in range(1, 13) if a not in (2, 4)]
PRIMOS = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]


def C1(a, p):
    return (6 * a * (a * a - 4) * (a * a - 16)) % p != 0


def C2(a, p):
    return p > a + 4


# control positivo: a = 1 es la Clebsch publicada
malc = [p for p in PRIMOS if contar(1, p) != p * p + (6 + chi(5, p)) * p + 1]
comprueba(u'0. CONTROL POSITIVO: a = 1 reproduce la Clebsch publicada p^2+(6+chi_p(5))p+1',
      not malc, u'falla en %s' % (malc,))
comprueba(u'   ... y la formula de la familia coincide ahi con la de la Clebsch',
      all(formula(1, p) == p * p + (6 + chi(5, p)) * p + 1 for p in PRIMOS if C1(1, p)))

fallos, buenos = [], []
for a in AES:
    for p in PRIMOS:
        m, f = contar(a, p), formula(a, p)
        (buenos if m == f else fallos).append((a, p, m, f))

comprueba(u'1. SIN condicion, la formula falla en %d de %d casos de la rejilla'
      % (len(fallos), len(fallos) + len(buenos)),
      len(fallos) > 0, u'primeros: %s' % (fallos[:5],))
comprueba(u'   ... e incluye el caso a = p = 5: medido 76, formula 31',
      (5, 5, 76, 31) in fallos)

# C1
malC1 = [(a, p, m, f) for (a, p, m, f) in fallos if C1(a, p)]
comprueba(u'2. la condicion C1  (p no divide 6a(a^2-4)(a^2-16))  NO deja pasar ningun fallo',
      not malC1, u'fallos que C1 no filtra: %s' % (malC1[:4],))
cubiertos1 = [(a, p) for (a, p, _, _) in buenos if C1(a, p)]
comprueba(u'   ... y no es vacua: cubre %d casos buenos de la rejilla' % len(cubiertos1),
      len(cubiertos1) > 40)

# C2
malC2 = [(a, p, m, f) for (a, p, m, f) in fallos if C2(a, p)]
comprueba(u'3. la condicion C2  (p > a+4)  tampoco deja pasar ningun fallo',
      not malC2, u'fallos que C2 no filtra: %s' % (malC2[:4],))
cubiertos2 = [(a, p) for (a, p, _, _) in buenos if C2(a, p)]
comprueba(u'   ... y cubre %d casos buenos' % len(cubiertos2), len(cubiertos2) > 30)

# relacion entre las dos
c2_no_c1 = [(a, p) for a in AES for p in PRIMOS if C2(a, p) and not C1(a, p)]
c1_no_c2 = [(a, p) for a in AES for p in PRIMOS if C1(a, p) and not C2(a, p)]
comprueba(u'4. C2 implica C1 en la rejilla (p > a+4 fuerza p no divide el producto)',
      not c2_no_c1, u'contraejemplos: %s' % (c2_no_c1[:4],))
comprueba(u'   ... pero C1 es ESTRICTAMENTE mas debil: %d casos la cumplen sin cumplir C2'
      % len(c1_no_c2), len(c1_no_c2) > 0, u'p.ej. %s' % (c1_no_c2[:4],))

# donde falla exactamente
comprueba(u'5. todo fallo tiene p dividiendo 6a(a^2-4)(a^2-16) --- la hipotesis es suficiente',
      all(not C1(a, p) for (a, p, _, _) in fallos))
con_p_a = sorted((a, p) for a in AES for p in PRIMOS if a % p == 0)
comprueba(u'6. los fallos son EXACTAMENTE los (a, p) de la rejilla con p | a --- la hipotesis no es necesaria',
      sorted((a, p) for (a, p, _, _) in fallos) == con_p_a, u'p | a en la rejilla: %s' % (con_p_a,))
sobra = [(a, p) for (a, p, _, _) in buenos if not C1(a, p)]
comprueba(u'   ... y hay pares excluidos por C1 donde la formula acierta: %d' % len(sobra),
      len(sobra) > 0, u'p.ej. %s' % (sobra[:4],))

# el paso de once rectas a veintisiete
from itertools import combinations as _comb
from sympy.combinatorics import Permutation, PermutationGroup

RECTAS = ([(u'a', i) for i in range(6)] + [(u'b', i) for i in range(6)]
          + [(u'c', i, j) for i, j in _comb(range(6), 2)])
IDX = {r: k for k, r in enumerate(RECTAS)}


def cortan(u, v):
    if u == v or (u[0] == v[0] and u[0] in u'ab'):
        return False
    if set((u[0], v[0])) == set(u'ab'):
        return u[1] != v[1]
    if u[0] == v[0] == u'c':
        return not (set(u[1:]) & set(v[1:]))
    c, o = (u, v) if u[0] == u'c' else (v, u)
    return o[1] in c[1:]


def perm(f):
    return Permutation([IDX[f(r)] for r in RECTAS])


def de_sigma(s):
    def f(r):
        if r[0] in u'ab':
            return (r[0], s[r[1]])
        i, j = sorted((s[r[1]], s[r[2]]))
        return (u'c', i, j)
    return f


def cremona(r):
    u"""La transformacion cuadratica en los puntos 0, 1, 2."""
    P, Q = {0, 1, 2}, {3, 4, 5}
    if r[0] == u'a' and r[1] in P:
        return (u'c',) + tuple(sorted(P - {r[1]}))
    if r[0] == u'c' and set(r[1:]) <= P:
        return (u'a', (P - set(r[1:])).pop())
    if r[0] == u'b' and r[1] in Q:
        return (u'c',) + tuple(sorted(Q - {r[1]}))
    if r[0] == u'c' and set(r[1:]) <= Q:
        return (u'b', (Q - set(r[1:])).pop())
    return r


GEN = [perm(de_sigma([1, 0, 2, 3, 4, 5])), perm(de_sigma([1, 2, 3, 4, 5, 0])), perm(cremona)]
comprueba(u'7. los tres generadores preservan la incidencia de las 27 rectas',
      all(cortan(RECTAS[i], RECTAS[j]) == cortan(RECTAS[g.array_form[i]], RECTAS[g.array_form[j]])
          for g in GEN for i in range(27) for j in range(27)))
W = PermutationGroup(GEN)
comprueba(u'   ... y generan un grupo de orden 51840 = |W(E_6)|', W.order() == 51840,
      u'orden %d' % W.order())
L0 = IDX[(u'a', 0)]
EST = W.stabilizer(L0)
VEC = [k for k in range(27) if cortan(RECTAS[L0], RECTAS[k])]
comprueba(u'   ... el estabilizador de una recta tiene orden 1920 = |W(D_5)|, y la cortan 10 rectas',
      EST.order() == 1920 and len(VEC) == 10, u'orden %d, %d rectas' % (EST.order(), len(VEC)))
NUCLEO = EST
for k in VEC:
    NUCLEO = NUCLEO.stabilizer(k)
comprueba(u'   ... y actua FIELMENTE sobre esas diez: el nucleo es trivial', NUCLEO.order() == 1,
      u'orden del nucleo %d' % NUCLEO.order())

L.append(u'')
L.append(u'  Los fallos de la formula sin condicion, y por que:')
L.append(u'  a    p    medido   formula   p | 6a(a^2-4)(a^2-16)')
L.append(u'  ' + u'-' * 56)
for (a, p, m, f) in fallos:
    L.append(u'  %-4d %-4d %-8d %-9d %s' % (a, p, m, f, u'si' if not C1(a, p) else u'NO'))

L.append(u'')
L.append(u'VEREDICTO')
L.append(u'  La formula de la Proposicion 11 NECESITA rango.  Las dos condiciones C1 y C2')
L.append(u'  funcionan las dos sobre la rejilla, y C2 (p > a+4) implica C1')
L.append(u'  pero es mas fuerte de lo necesario: C1 es la que se enuncia, con C2 como')
L.append(u'  la version comoda "en el contexto de las particiones".  C1 es suficiente y no')
L.append(u'  necesaria: en la rejilla la formula solo falla cuando p | a.')
L.append(u'  El estabilizador W(D_5) de una recta actua fielmente sobre las diez que la cortan,')
L.append(u'  asi que el cuerpo de esas diez es el de las veintisiete.')

print(u'\n'.join(L))
print(u'')
print(u'=' * 72)
print(u'VERIFICACION FAMILIA V_(a,1,1,1,1):  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 72)
sys.exit(1 if MAL else 0)
