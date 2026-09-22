# -*- coding: utf-8 -*-
u"""P3d_columna_n3.py -- la columna n = 3: la identidad de thm:columna contra fuerza bruta.

Con n = 3 y e_1 = 0 las identidades de Newton dan s_k como polinomio en e_2, e_3:

  s_1 = 0,  s_2 = -2e_2,  s_3 = 3e_3,   y   s_k = -e_2 s_{k-2} + e_3 s_{k-3}.

El escalado T -> uT manda (e_2, e_3) -> (u^2 e_2, u^3 e_3), asi que s_k es cuasi-homogeneo y
todo factor de s_k distinto de e_2 y e_3 es un polinomio Q_k en u = e_2^3/e_3^2 (forma
primitiva entera).  El enunciado que se comprueba (thm:columna): para todo primo
p que no divide a 6 k D_k, con D_k = lc(Q_k) disc(Q_k) Res(Q_k, u(4u+27)) (eq:delta; D_k = 1 si
Q_k es constante),

  6 C_k(3,p)/(p-1) = [e_2 | s_k](1 + chi_p(-3)) + 3 [e_3 | s_k]
                     + 6 #{u_0 en F_p : Q_k(u_0) = 0, f_{u_0} escinde del todo},

con f_u(Y) = Y^3 + uY - u, y C_k(3,p) el numero de 3-subconjuntos de F_p con suma 0 y suma
de potencias k-esimas 0.

Se comprueba:
  1. la recursion de Newton frente a s_4, s_5, s_6, y la factorizacion exacta de s_k para
     3 <= k <= 30: s_k = c e_2^a e_3^b Q_k(e_2^3/e_3^2) e_3^(2 deg Q_k), con deg Q_k =
     floor(k/6) - [k = 1 mod 6], D_k distinto de 0, y Q_k constante exactamente en k = 3,4,5,7;
  2. la identidad de thm:columna frente al recuento directo, para 3 <= k <= 14 y todo primo
     5 <= p <= 257 con p que no divide a 6 k D_k; se imprimen, para cada k, los primos
     excluidos y cuantos quedan;
  3. la forma cerrada C_7(3,p) = (p-1)(4 + chi_p(-3))/6 de tab:formas;
  4. rem:ventana, k = 9: en que primos aparece el segundo valor;
  5. controles negativos: la misma formula con el cubico equivocado Y^3 - uY - u, y sin la
     condicion de escision, tiene que FALLAR frente al recuento directo.

Uso:  python P3d_columna_n3.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

import sympy as sp

OK = 0
MAL = 0
L = []

KFACT = 30          # factorizaciones exactas de s_k
KMIN, KMAX = 3, 14  # identidad contra fuerza bruta
PMIN, PMAX = 5, 257


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


def primos(lo, hi):
    return [n for n in range(max(2, lo), hi + 1)
            if all(n % d for d in range(2, int(n ** 0.5) + 1))]


def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def cuentas(p, ks):
    u"""{k: #{T de F_p, |T| = 3, suma 0, suma de k-esimas 0}} para todas las k a la vez."""
    pot = {k: [pow(t, k, p) for t in range(p)] for k in ks}
    c = dict((k, 0) for k in ks)
    for a in range(p):
        for b in range(a + 1, p):
            x = (-a - b) % p
            if x <= b:
                continue
            for k in ks:
                pk = pot[k]
                if (pk[a] + pk[b] + pk[x]) % p == 0:
                    c[k] += 1
    return c


# ---------------------------------------------------------------------------
e2, e3, u = sp.symbols('e2 e3 u')
S = {1: sp.Integer(0), 2: -2 * e2, 3: 3 * e3}
for k in range(4, KFACT + 1):
    S[k] = sp.expand(-e2 * S[k - 2] + e3 * S[k - 3])


def descompone(k):
    u"""(a, b, Q_k): e_2^a e_3^b es la parte monomial de s_k y Q_k(u) el resto, primitivo."""
    facs = sp.factor_list(S[k])[1]
    a = b = 0
    resto = sp.Integer(1)
    for base, m in facs:
        bb = sp.expand(base)
        if bb == e2:
            a += m
        elif bb == e3:
            b += m
        else:
            resto *= bb ** m
    if resto == 1:
        return a, b, sp.Poly(sp.Integer(1), u)
    P = sp.Poly(sp.expand(resto), e2, e3)
    expr = sp.Integer(0)
    for (i, j), co in zip(P.monoms(), P.coeffs()):
        if i % 3 or j % 2:
            return a, b, None          # no seria polinomio en u: lo detecta la comprobacion
        expr += co * u ** (i // 3)
    return a, b, sp.Poly(sp.primitive(sp.expand(expr))[1], u)


def D(Q):
    if Q.degree() == 0:
        return sp.Integer(1)
    return sp.Integer(Q.LC() * sp.discriminant(Q) *
                      sp.resultant(Q.as_expr(), u * (4 * u + 27), u))


# ---------------------------------------------------------------------------
sec(u'1.  Newton con e_1 = 0 y la factorizacion exacta de s_k, 3 <= k <= %d' % KFACT)
# ---------------------------------------------------------------------------
comprueba(u's_4 = 2 e_2^2', sp.expand(S[4] - 2 * e2 ** 2) == 0)
comprueba(u's_5 = -5 e_2 e_3', sp.expand(S[5] + 5 * e2 * e3) == 0)
comprueba(u's_6 = -2 e_2^3 + 3 e_3^2', sp.expand(S[6] - (-2 * e2 ** 3 + 3 * e3 ** 2)) == 0)

DESC = {}
malos_poly, malos_rec, malos_grado, malos_D = [], [], [], []
for k in range(3, KFACT + 1):
    a, b, Q = descompone(k)
    DESC[k] = (a, b, Q)
    if Q is None:
        malos_poly.append(k)
        continue
    d = Q.degree()
    # reconstruccion: s_k / (e_2^a e_3^b) es proporcional a e_3^(2d) Q_k(e_2^3/e_3^2)
    rec = sp.expand(e3 ** (2 * d) * Q.as_expr().subs(u, e2 ** 3 / e3 ** 2))
    cociente = sp.cancel(S[k] / (e2 ** a * e3 ** b * rec))
    if not cociente.is_number or cociente == 0:
        malos_rec.append(k)
    if d != k // 6 - (1 if k % 6 == 1 else 0):
        malos_grado.append((k, d))
    if D(Q) == 0:
        malos_D.append(k)
    L.append(u'   k = %-2d  e2^%d e3^%d   deg Q_k = %d   Q_k = %s'
             % (k, a, b, d, Q.as_expr() if d <= 2 else u'(grado %d)' % d))
comprueba(u'el factor propio de s_k es polinomio en u = e_2^3/e_3^2, 3 <= k <= %d' % KFACT,
          not malos_poly, u'fallos: %s' % malos_poly)
comprueba(u's_k = c e_2^a e_3^b e_3^(2d) Q_k(e_2^3/e_3^2) con c constante, 3 <= k <= %d' % KFACT,
          not malos_rec, u'fallos: %s' % malos_rec)
comprueba(u'deg Q_k = floor(k/6) - [k = 1 mod 6], 3 <= k <= %d' % KFACT,
          not malos_grado, u'fallos: %s' % malos_grado)
comprueba(u'D_k distinto de 0, 3 <= k <= %d' % KFACT, not malos_D, u'fallos: %s' % malos_D)
ceros = [k for k in range(3, KFACT + 1) if DESC[k][2] is not None and DESC[k][2].degree() == 0]
comprueba(u'Q_k constante exactamente en k = 3,4,5,7 (k <= %d)' % KFACT, ceros == [3, 4, 5, 7],
          u'%s' % ceros)


# ---------------------------------------------------------------------------
# la prediccion de thm:columna
# ---------------------------------------------------------------------------
def raices(coefs, p):
    u"""Raices en F_p del polinomio de coeficientes enteros (grado descendente)."""
    r = []
    for x in range(p):
        v = 0
        for c in coefs:
            v = (v * x + c) % p
        if v == 0:
            r.append(x)
    return r


def escinde(a, b, p):
    u"""Y^3 + aY + b tiene tres raices distintas en F_p."""
    return len(raices([1, 0, a % p, b % p], p)) == 3


def predice(k, p, modo=u'bueno'):
    a, b, Q = DESC[k]
    tot = 0
    if a:
        tot += 1 + chi(-3, p)
    if b:
        tot += 3
    if Q.degree() > 0:
        cs = [int(c) for c in Q.all_coeffs()]
        for u0 in raices(cs, p):
            if modo == u'bueno':
                tot += 6 if escinde(u0, -u0, p) else 0          # f_u = Y^3 + uY - u
            elif modo == u'cubico_malo':
                tot += 6 if escinde(-u0, -u0, p) else 0         # Y^3 - uY - u
            elif modo == u'sin_escision':
                tot += 6
    return tot


# ---------------------------------------------------------------------------
sec(u'2.  thm:columna contra el recuento directo, %d <= k <= %d, %d <= p <= %d'
    % (KMIN, KMAX, PMIN, PMAX))
# ---------------------------------------------------------------------------
PS = primos(PMIN, PMAX)
KS = list(range(KMIN, KMAX + 1))
L.append(u'   primos 5 <= p <= %d: %d' % (PMAX, len(PS)))
C = {}
for p in PS:
    C[p] = cuentas(p, KS)

total_casos = 0
fuera_de_hipotesis = []   # (k, p) excluidos donde la identidad falla: solo se informa
for k in KS:
    Dk = D(DESC[k][2])
    excl = [p for p in PS if (6 * k * Dk) % p == 0]
    dentro = [p for p in PS if p not in excl]
    fallos = []
    for p in dentro:
        lhs = sp.Rational(6 * C[p][k], p - 1)
        if lhs != predice(k, p):
            fallos.append((p, lhs, predice(k, p)))
    total_casos += len(dentro)
    for p in excl:
        if sp.Rational(6 * C[p][k], p - 1) != predice(k, p):
            fuera_de_hipotesis.append((k, p))
    comprueba(u'k = %-2d  D_k = %-14s excluidos %-16s quedan %d primos'
              % (k, Dk, excl, len(dentro)), not fallos, u'fallos: %s' % fallos[:3])
L.append(u'')
L.append(u'   casos (k, p) comprobados en total: %d' % total_casos)
L.append(u'   (informativo) primos excluidos donde la identidad falla: %s' % fuera_de_hipotesis)

# ---------------------------------------------------------------------------
sec(u'3.  tab:formas: C_7(3,p) = (p-1)(4 + chi_p(-3))/6')
# ---------------------------------------------------------------------------
P7 = [p for p in PS if 42 % p]
f7 = [p for p in P7 if 6 * C[p][7] != (p - 1) * (4 + chi(-3, p))]
comprueba(u'C_7(3,p) = (p-1)(4+chi_p(-3))/6 en los %d primos 5 <= p <= %d con p no | 42'
          % (len(P7), PMAX), not f7, u'fallos: %s' % f7)
f7b = [p for p in PS if 6 * C[p][7] != (p - 1) * (4 + chi(-3, p))]
comprueba(u'y en los %d primos 5 <= p <= %d, p = 7 incluido (fuera de la hipotesis)'
          % (len(PS), PMAX), not f7b, u'fallos: %s' % f7b)

# ---------------------------------------------------------------------------
sec(u'4.  rem:ventana: k = 9, Q_9 = %s' % DESC[9][2].as_expr())
# ---------------------------------------------------------------------------
D9 = D(DESC[9][2])
P9 = [p for p in PS if (6 * 9 * D9) % p]
vals9 = {}
for p in P9:
    vals9.setdefault(sp.Rational(6 * C[p][9], p - 1), []).append(p)
for v in sorted(vals9):
    L.append(u'   6C/(p-1) = %-3s en %d primos: %s' % (v, len(vals9[v]), vals9[v]))
seg = vals9.get(sp.Integer(9), [])
comprueba(u'k = 9 toma exactamente los valores 3 y 9 en p <= %d, p no | 6*9*D_9 = %d'
          % (PMAX, 6 * 9 * D9), sorted(vals9) == [3, 9])
comprueba(u'el segundo valor aparece exactamente en 71, 83, 97, 131, 193, 257',
          seg == [71, 83, 97, 131, 193, 257], u'%s' % seg)
v53 = [p for p in P9 if p <= 53]
comprueba(u'en los %d primos p <= 53 con p no | %d toma un solo valor' % (len(v53), 6 * 9 * D9),
          set(sp.Rational(6 * C[p][9], p - 1) for p in v53) == {3}, u'%s' % v53)
t53 = [p for p in PS if p <= 53]
comprueba(u'y tambien en los %d primos 5 <= p <= 53 sin exclusion (el recuento directo)'
          % len(t53), set(sp.Rational(6 * C[p][9], p - 1) for p in t53) == {3}, u'%s' % t53)

# ---------------------------------------------------------------------------
sec(u'5.  Controles negativos: la formula con otra condicion tiene que fallar')
# ---------------------------------------------------------------------------
for modo, nombre in ((u'cubico_malo', u'cubico Y^3 - uY - u en lugar de f_u'),
                     (u'sin_escision', u'contar las raices de Q_k sin exigir escision')):
    fallos = []
    for k in KS:
        if DESC[k][2].degree() == 0:
            continue
        Dk = D(DESC[k][2])
        for p in PS:
            if (6 * k * Dk) % p == 0:
                continue
            if sp.Rational(6 * C[p][k], p - 1) != predice(k, p, modo):
                fallos.append((k, p))
    comprueba(u'control negativo: %s falla' % nombre, len(fallos) > 0,
              u'%d casos (k, p) en desacuerdo; primeros %s' % (len(fallos), fallos[:3]))

print(u'\n'.join(L))
print(u'')
print(u'=' * 70)
print(u'VERIFICACION COLUMNA n=3:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 70)
sys.exit(1 if MAL else 0)
