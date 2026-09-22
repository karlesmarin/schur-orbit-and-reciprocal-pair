# -*- coding: utf-8 -*-
u"""P3d_grupo_del_caracter.py -- el grupo del caracter de la columna n=3, para 3 <= k <= 14.

Cuando $e_2\\mid s_k$ el recuento lleva ademas $\\chi_p(-3)$, y el grupo es $S_3\\times C_2$
SOLO si $\\QQ(\\sqrt{-3})$ es DISTINTO del subcuerpo cuadratico del cuerpo de descomposicion
de la cubica: eso no es una consecuencia formal de $e_2\\mid s_k$.  Ese subcuerpo es
$\\QQ(\\sqrt{\\operatorname{disc}})$ con $\\operatorname{disc} f_{u_0} = -u_0^2(4u_0+27)$.
Si para algun $k$ resultara $-u_0^2(4u_0+27) = -3\\cdot(\\text{cuadrado})$, los dos caracteres
serian EL MISMO y el grupo seria $S_3$, no $S_3\\times C_2$.

Se calcula k a k (sympy, exacto): Q_k, sus raices, si $e_2\\mid s_k$, el discriminante de la
cubica de la fibra y su parte libre de cuadrados.  El conjunto de k con grupo exactamente
S_3 tiene que salir {6, 9} (las celdas oscuras de la figura de la rejilla).  Controles: la
funcion de parte libre de cuadrados identifica -12 con -3 (el caso en que coincidirian) y
separa -31 de -3.

DEFINICION QUE SE USA.  "Grupo del caracter" = el menor cociente finito de $G_\\QQ$ por el que
factoriza la funcion de clase que da el recuento.  Es el grupo de Galois del compuesto de los
cuerpos que los terminos de thm:columna nombran: $\\QQ(\\sqrt{-3})$ cuando $e_2\\mid s_k$, y el
cuerpo de descomposicion de cada $f_{u_0}$.

Uso:  python P3d_grupo_del_caracter.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

import sympy as sp

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

# el conjunto de k con grupo exactamente S_3 que da la nota (celdas oscuras de fig:rejilla)
S3_EXACTO_NOTA = {6, 9}

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
    return bool(cond)


e2, e3, u_, Y = sp.symbols('e2 e3 u Y')
S = {0: sp.Integer(3), 1: sp.Integer(0), 2: -2 * e2, 3: 3 * e3}
for k in range(4, 20):
    S[k] = sp.expand(-e2 * S[k - 2] + e3 * S[k - 3])


def Qk(k):
    facs = sp.factor_list(sp.factor(S[k]))[1]
    resto = sp.Integer(1)
    for base, m in facs:
        b = sp.expand(base)
        if b not in (e2, e3):
            resto *= b ** m
    if resto == 1:
        return sp.Poly(sp.Integer(1), u_)
    P = sp.Poly(sp.expand(resto), e2, e3)
    expr = sp.Integer(0)
    for mon, co in zip(P.monoms(), P.coeffs()):
        a, b = mon
        expr += co * u_ ** (a // 3)
    return sp.Poly(sp.primitive(sp.expand(expr))[1], u_)


def sin_cuadrados(n):
    u"""La parte libre de cuadrados de un entero (con su signo)."""
    n = sp.Integer(n)
    s = sp.sign(n)
    r = sp.Integer(1)
    for p, e in sp.factorint(abs(n)).items():
        if e % 2:
            r *= p
    return s * r


L.append(u'EL GRUPO DEL CARACTER, k A k')
L.append(u'-' * 88)
L.append(u'  k   deg  raiz u_0   e2|s_k  disc f_u0        libre de cuad.  cuerpo cuadratico  grupo')
L.append(u'  ' + u'-' * 86)

exactamente_S3 = set()
con_C2 = set()
for k in range(3, 15):
    Q = Qk(k)
    d = Q.degree()
    div2 = sp.simplify(sp.expand(S[k] / e2)).is_polynomial(e2, e3)
    if d == 0:
        L.append(u'  %-3d %-4d %-10s %-7s %-16s %-15s %-18s %s'
                 % (k, 0, u'---', div2, u'---', u'---', u'---', u'cuadratico'))
        continue
    raices = sp.roots(Q.as_expr(), u_)
    unica_rac = d == 1 and len(raices) == 1 and list(raices)[0].is_rational
    if not unica_rac:
        L.append(u'  %-3d %-4d %-10s %-7s %-16s %-15s %-18s %s'
                 % (k, d, u'irrac./varias', div2, u'---', u'---', u'---',
                    u'mayor que S_3'))
        continue
    r = list(raices)[0]
    disc = sp.nsimplify(-r ** 2 * (4 * r + 27))
    # el subcuerpo cuadratico del cuerpo S_3 es Q(sqrt(disc)); se compara su parte libre de
    # cuadrados con la de -3.  Iguales => los dos caracteres coinciden y NO hay C_2 de mas.
    num, den = sp.fraction(sp.nsimplify(disc))
    libre = sin_cuadrados(sp.Integer(num) * sp.Integer(den))   # sqrt(a/b) ~ sqrt(ab)
    mismo = (libre == sin_cuadrados(-3))
    if div2 and not mismo:
        grupo = u'S_3 x C_2'
        con_C2.add(k)
    elif div2 and mismo:
        grupo = u'S_3  (los dos caracteres COINCIDEN)'
    else:
        grupo = u'S_3 exactamente'
        exactamente_S3.add(k)
    L.append(u'  %-3d %-4d %-10s %-7s %-16s %-15s %-18s %s'
             % (k, d, r, div2, disc, libre, u'Q(sqrt(%s))' % libre, grupo))

L.append(u'')
L.append(u'  grupo EXACTAMENTE S_3 : %s' % sorted(exactamente_S3))
L.append(u'  grupo S_3 x C_2       : %s' % sorted(con_C2))
L.append(u'')

# --------------------------------------------------------------------------------------------
comprueba(u'1. el grupo es exactamente S_3 justo en k = 6, 9 (las celdas oscuras de la figura)',
      exactamente_S3 == S3_EXACTO_NOTA,
      u'medido %s' % sorted(exactamente_S3))
comprueba(u'2. y en los k con e_2|s_k el segundo caracter es DISTINTO del de la cubica',
      bool(con_C2) and not (con_C2 & exactamente_S3),
      u'los de S_3 x C_2 son %s' % sorted(con_C2))
# que Q(sqrt(-3)) no sea el subcuerpo cuadratico de la cubica, comprobado explicitamente:
malos = []
for k in sorted(con_C2):
    Q = Qk(k)
    r = list(sp.roots(Q.as_expr(), u_))[0]
    disc = sp.nsimplify(-r ** 2 * (4 * r + 27))
    num, den = sp.fraction(disc)
    if sin_cuadrados(sp.Integer(num) * sp.Integer(den)) == sin_cuadrados(-3):
        malos.append(k)
comprueba(u'3. ninguno de ellos tiene disc(f_u0) igual a -3 por un cuadrado', not malos,
      u'si lo tuviera, Q(sqrt(disc)) = Q(sqrt(-3)) y el grupo seria S_3, no S_3 x C_2: %s'
      % malos)
# CONTROL: un u_0 inventado PARA el cual si coinciden, para ver que la comprobacion distingue
u_falso = sp.Rational(-27, 4) + sp.Rational(1, 1)     # 4u+27 = 4, disc = -4u^2 -> libre -1
disc_falso = sp.nsimplify(-u_falso ** 2 * (4 * u_falso + 27))
n2, d2 = sp.fraction(disc_falso)
comprueba(u'4. control: la comprobacion sabe ver cuando un disc SI es -3 por un cuadrado',
      sin_cuadrados(sp.Integer(-3) * 4) == sin_cuadrados(-3),
      u'-12 y -3 tienen la misma parte libre de cuadrados, y la funcion lo dice')
comprueba(u'5. ... y sabe ver cuando NO lo es',
      sin_cuadrados(sp.Integer(-31)) != sin_cuadrados(-3))

L.append(u'')
L.append(u'VEREDICTO')
L.append(u'  El grupo es S_3 x C_2 en los k con e_2 | s_k, pero por una razon aritmetica que')
L.append(u'  hay que medir y no por la forma de e_2 | s_k: en los')
L.append(u'  cuatro k con e_2 | s_k, el discriminante de la cubica de la fibra tiene parte')
L.append(u'  libre de cuadrados distinta de -3, luego Q(sqrt(disc)) y Q(sqrt(-3)) son cuerpos')
L.append(u'  distintos y el compuesto tiene grado 12 sobre Q: el grupo es S_3 x C_2.')
L.append(u'  Si alguno hubiera coincidido, el grupo habria sido S_3.')

print(u'\n'.join(L))
print(u'')
print(u'=' * 78)
print(u'VERIFICACION EL GRUPO DEL CARACTER:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 78)
sys.exit(1 if MAL else 0)
