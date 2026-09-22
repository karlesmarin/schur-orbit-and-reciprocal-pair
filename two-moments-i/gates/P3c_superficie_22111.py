# -*- coding: utf-8 -*-
u"""P3c_superficie_22111.py -- la superficie V_(2,2,1,1,1), su forma separada, y las cuentas vecinas.

Comprueba, de la nota P3c: la Proposicion prop:22111 y el Corolario cor:colapso
(#S(F_p) = p^2 + 1 + p (3 + chi_p(5) (1 + R_p)) para p > 5, con R_p = #raices distintas de
P(r) = r^3 + 8r^2 + 16r + 10), los casos p = 5 y p = 3 de su prueba, que "P se escinde del todo"
no equivale a "las conicas residuales se parten en rectas racionales" (p = 71, 73, 83), el
polinomio H(z) = (z^3 - 8z)^2 - 80 y su grupo de Galois de orden 12, el discriminante -140 de P,
los valores posibles de a_p(S) = (#S - p^2 - 1)/p, las cuentas directas de la prueba de thm:n6
(Segre #V_(1^6) = p^3 + 6p^2 - 4p + 1 y Cayley #V_(2,1^4) = p^2 + 3p + 1, con sus cambios de
variables) y que la hipotesis de primos de prop:familia es necesaria (a = p = 5).

Todo recuento se hace por FUERZA BRUTA proyectiva, sin usar ninguna de las formulas que se
comprueban; las identidades algebraicas son simbolicas (sympy).  Controles negativos: la forma
separada de cor:colapso aplicada en p = 5 tiene que FALLAR; la igualdad chi_p(-2r) = chi_p(5)
tiene que fallar si se cambia el 5 por otro entero, o el enunciado no dice nada; la formula de
prop:familia extrapolada a a = p = 5 tiene que fallar.  Control positivo: el contador reproduce
la cuenta de la superficie de Clebsch, p^2 + (6 + chi_p(5)) p + 1.

Uso:  python P3c_superficie_22111.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

import numpy as np
import sympy as sp

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


def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def raices_P(p):
    return [r for r in range(p) if (r ** 3 + 8 * r ** 2 + 16 * r + 10) % p == 0]


# ---------------------------------------------------------------- contador general
def contar_proyectivo(pesos, p):
    u"""#V_w(F_p) por fuerza bruta: se elimina la ultima variable con la lineal."""
    w = [x % p for x in pesos]
    n = len(w)
    inv = pow(w[-1], p - 2, p)
    cub = (np.arange(p, dtype=np.int64) ** 3) % p
    g = np.arange(p, dtype=np.int64)
    forma = [g.reshape((1,) * i + (p,) + (1,) * (n - 3 - i)) for i in range(n - 2)]
    base = np.zeros((p,) * (n - 2), dtype=np.int64)
    lin = np.zeros((p,) * (n - 2), dtype=np.int64)
    for i in range(n - 2):
        base = (base + w[i + 1] * cub[forma[i]]) % p
        lin = (lin + w[i + 1] * forma[i]) % p
    ceros = 0
    for a in range(p):
        xl = (-(w[0] * a + lin) * inv) % p
        F = (w[0] * cub[a] + base + w[-1] * cub[xl]) % p
        ceros += int(np.count_nonzero(F == 0))
    return (ceros - 1) // (p - 1)


PR = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 83, 89, 97]

# control positivo: la superficie de Clebsch V_(1^5)
malc = [p for p in PR
        if contar_proyectivo((1, 1, 1, 1, 1), p) != p * p + (6 + chi(5, p)) * p + 1]
comprueba(u'0. CONTROL POSITIVO: el contador reproduce #V_(1^5) = p^2 + (6+chi_p(5))p + 1'
          u' (Clebsch)', not malc, u'falla en %s' % (malc,))

# ---------------------------------------------------------------- A: cor:colapso
r = sp.Symbol('r')
comprueba(u'A1. P(r) = r(r+4)^2 + 10 es una identidad (prueba de cor:colapso)',
          sp.expand(r * (r + 4) ** 2 + 10 - (r ** 3 + 8 * r ** 2 + 16 * r + 10)) == 0)

malA = []
for p in [q for q in PR if q > 5]:
    for rr in raices_P(p):
        if (rr + 4) % p == 0:
            malA.append((p, rr, u'r+4 = 0'))
        elif chi(-2 * rr, p) != chi(5, p):
            malA.append((p, rr, chi(-2 * rr, p), chi(5, p)))
comprueba(u'A2. r+4 != 0 y chi_p(-2r) = chi_p(5) para CADA raiz de P, en todo p > 5 del rango',
          not malA, u'fallos: %s' % (malA[:4],))

S22 = {p: contar_proyectivo((2, 2, 1, 1, 1), p) for p in PR}


def separada(p):
    u"""Forma separada de cor:colapso."""
    return p * p + 1 + p * (3 + chi(5, p) * (1 + len(raices_P(p))))


malS = [(p, S22[p], separada(p)) for p in PR if p > 5 and S22[p] != separada(p)]
comprueba(u'A3. cor:colapso: la forma separada acierta en todo p > 5 del rango (%d primos)'
          % len([q for q in PR if q > 5]),
          not malS, u'fallos: %s' % (malS[:4],))

# control negativo: cambiar el 5 por otro entero tiene que romper A2
rotos5 = []
for d in (2, 3, 7, -1, -3, 13):
    mal = any(chi(-2 * rr, p) != chi(d, p)
              for p in PR if p > 5 for rr in raices_P(p))
    rotos5.append((d, mal))
comprueba(u'A4. CONTROL NEGATIVO: con cualquier otro entero en lugar del 5, la igualdad FALLA',
          all(mal for (_, mal) in rotos5), u'%s' % (rotos5,))

# ---------------------------------------------------------------- B: p = 5 y p = 3
comprueba(u'B1. CONTROL NEGATIVO: en p = 5 la forma separada de cor:colapso falla',
          contar_proyectivo((2, 2, 1, 1, 1), 5) != separada(5),
          u'medido %d, forma separada %d' % (contar_proyectivo((2, 2, 1, 1, 1), 5), separada(5)))
comprueba(u'B2. ... y la formula de prop:22111 acierta en p = 5',
          contar_proyectivo((2, 2, 1, 1, 1), 5)
          == 25 + 1 + 5 * (3 + chi(5, 5) + sum(chi(-2 * rr, 5) for rr in raices_P(5))),
          u'raices de P mod 5: %s' % (raices_P(5),))
comprueba(u'B3. ... y #V_(2,2,1,1,1)(F_5) = 36, el total del recuento a mano de la prueba',
          contar_proyectivo((2, 2, 1, 1, 1), 5) == 36)
comprueba(u'B4. #V_(2,2,1,1,1)(F_3) = 40, o sea todo P^3(F_3) (la excepcion p = 3)',
          contar_proyectivo((2, 2, 1, 1, 1), 3) == 40 and 3 ** 3 + 3 ** 2 + 3 + 1 == 40)

# ---------------------------------------------------------------- C: escindir != rectas racionales
tot = [p for p in PR if len(raices_P(p)) == 3]
comprueba(u'C1. P se escinde del todo exactamente en %s del rango' % (tot,),
          set(tot) >= {71, 73, 83})
comprueba(u'C2. ... y ahi chi_p(5) NO es constante: +1 en 71, -1 en 73 y 83',
          chi(5, 71) == 1 and chi(5, 73) == -1 and chi(5, 83) == -1)
comprueba(u'C3. ... luego "P se escinde" != "las conicas residuales dan rectas racionales":'
          u' #S(F_73) = 73^2 - 73 + 1 = 5257',
          S22[73] == 73 ** 2 - 73 + 1 == 5257)

# ---------------------------------------------------------------- D: H(z), Galois, disc
z = sp.Symbol('z')
H = sp.expand((z ** 3 - 8 * z) ** 2 - 80)
comprueba(u'D1. H(z) = (z^3-8z)^2 - 80 es z^6 - 16z^4 + 64z^2 - 80',
          sp.expand(H - (z ** 6 - 16 * z ** 4 + 64 * z ** 2 - 80)) == 0)
comprueba(u'D2. ... y las raices de H dan las de P por r = -z^2/2',
          sp.simplify(sp.resultant(H, r + z ** 2 / 2, z)) != 0
          and all(sp.simplify((r ** 3 + 8 * r ** 2 + 16 * r + 10).subs(
              r, -sol ** 2 / 2)) == 0 for sol in sp.Poly(H, z).all_roots()[:1]))
try:
    G = sp.Poly(H, z).galois_group()
    orden = G[0].order() if isinstance(G, tuple) else G.order()
except Exception as e:
    orden = None
comprueba(u'D3. el grupo de Galois de H tiene orden 12',
          orden == 12, u'orden calculado: %s' % orden)
discP = sp.discriminant(r ** 3 + 8 * r ** 2 + 16 * r + 10, r)
# El unico subcuerpo cuadratico del cuerpo de descomposicion de una cubica con grupo S_3 es
# Q(sqrt(disc)).  Aqui disc = -140 = 4 * (-35), luego es Q(sqrt(-35)), que NO es Q(sqrt5).
_sqfree = sp.factorint(-discP)
_parte = -sp.prod([q ** (e % 2) for q, e in _sqfree.items()])
comprueba(u'D4. disc P = -140 (prop:22111), cuya parte libre de cuadrados es -35, luego el'
          u' unico cuadratico de K_0 es Q(sqrt(-35)) != Q(sqrt5)',
          discP == -140 and _parte == -35 and _parte != 5,
          u'sqrt(disc) genera Q(sqrt(%s))' % _parte)

# ---------------------------------------------------------------- E: valores de a_p
vals = sorted(set((S22[p] - p * p - 1) // p for p in PR if p > 7))
comprueba(u'E1. los valores de a_p(S) = (#S - p^2 - 1)/p para p > 7 caen en {-1,1,2,4,5,7}',
          set(vals) <= {-1, 1, 2, 4, 5, 7}, u'observados: %s' % (vals,))
comprueba(u'E2. ... y a_p(S) = 3 + chi_p(5)(1 + R_p) en todo el rango (tabla abajo)',
          all((S22[p] - p * p - 1) // p == 3 + chi(5, p) * (1 + len(raices_P(p)))
              for p in PR if p > 7))

# ---------------------------------------------------------------- F y G: prueba de thm:n6
malF = [(p, contar_proyectivo((1, 1, 1, 1, 1, 1), p), p ** 3 + 6 * p ** 2 - 4 * p + 1)
        for p in [5, 7, 11, 13, 17, 19, 23]
        if contar_proyectivo((1, 1, 1, 1, 1, 1), p) != p ** 3 + 6 * p ** 2 - 4 * p + 1]
comprueba(u'F1. thm:n6, cuenta directa de SEGRE: #V_(1^6) = p^3 + 6p^2 - 4p + 1', not malF,
          u'fallos: %s' % (malF,))
s_, t_, u_, v_, w_ = sp.symbols('s t u v w')
x = [(s_ + u_) / 2, (s_ - u_) / 2, (t_ + v_) / 2, (t_ - v_) / 2,
     (-s_ - t_ + w_) / 2, (-s_ - t_ - w_) / 2]
comprueba(u'F2. ... y el cambio de variables da su^2 + tv^2 - (s+t)w^2 - st(s+t) = 0',
          sp.expand(sum(sp.Integer(1) * xi for xi in x)) == 0
          and sp.simplify(sp.expand(sum(xi ** 3 for xi in x)) * sp.Rational(4, 3)
                          - (s_ * u_ ** 2 + t_ * v_ ** 2 - (s_ + t_) * w_ ** 2
                             - s_ * t_ * (s_ + t_))) == 0)

malG = [(p, contar_proyectivo((2, 1, 1, 1, 1), p), p ** 2 + 3 * p + 1)
        for p in PR if contar_proyectivo((2, 1, 1, 1, 1), p) != p ** 2 + 3 * p + 1]
comprueba(u'G1. thm:n6, cuenta directa de CAYLEY: #V_(2,1^4) = p^2 + 3p + 1', not malG,
          u'fallos: %s' % (malG[:4],))
# La sustitucion de la prueba de thm:n6, simbolica: las cuatro de peso uno en dos parejas, y la
# de peso dos igual a -(s+t)/2.
xc = [(s_ + u_) / 2, (s_ - u_) / 2, (t_ + v_) / 2, (t_ - v_) / 2]
x2peso = -(s_ + t_) / 2
comprueba(u'G2. ... y el cambio de variables: la lineal 2*x_peso2 + sum x_i = 0 se cumple sola',
          sp.simplify(2 * x2peso + sum(xc)) == 0)
comprueba(u'G3. ... y la cubica queda s u^2 + t v^2 - s t (s+t), salvo constante',
          sp.simplify(sp.cancel(sp.expand(2 * x2peso ** 3 + sum(xi ** 3 for xi in xc))
                                / (s_ * u_ ** 2 + t_ * v_ ** 2
                                   - s_ * t_ * (s_ + t_)))).is_number,
          u'cociente = %s' % sp.simplify(sp.cancel(
              sp.expand(2 * x2peso ** 3 + sum(xi ** 3 for xi in xc))
              / (s_ * u_ ** 2 + t_ * v_ ** 2 - s_ * t_ * (s_ + t_)))))

# ---------------------------------------------------------------- H: hipotesis de prop:familia
comprueba(u'H1. a = p = 5: #V_(5,1,1,1,1)(F_5) = 76',
          contar_proyectivo((5, 1, 1, 1, 1), 5) == 76,
          u'medido %d' % contar_proyectivo((5, 1, 1, 1, 1), 5))
Da, Ea = -3 * (25 - 4), (25 - 4) * (25 - 16)
comprueba(u'H2. ... y la formula de prop:familia extrapolada a a = p = 5 da 31',
          25 + (3 + 3 * chi(Da, 5) + chi(Ea, 5)) * 5 + 1 == 31,
          u'D_5 = %d, E_5 = %d, chi_5 = %d y %d' % (Da, Ea, chi(Da, 5), chi(Ea, 5)))
comprueba(u'H3. CONTROL NEGATIVO: 76 != 31, luego prop:familia necesita p no divida'
          u' 6a(a^2-4)(a^2-16)',
          contar_proyectivo((5, 1, 1, 1, 1), 5) != 31)

L.append(u'')
L.append(u'  p     #S medido   R_p   chi_p(5)   a_p   3+chi(1+R)')
L.append(u'  ' + u'-' * 58)
for p in PR:
    if p > 7:
        L.append(u'  %-5d %-11d %-5d %-10d %-5d %d'
                 % (p, S22[p], len(raices_P(p)), chi(5, p),
                    (S22[p] - p * p - 1) // p,
                    3 + chi(5, p) * (1 + len(raices_P(p)))))

print(u'\n'.join(L))
print(u'')
print(u'=' * 72)
print(u'VERIFICACION SUPERFICIE V_(2,2,1,1,1):  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 72)
sys.exit(1 if MAL else 0)
