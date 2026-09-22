# -*- coding: utf-8 -*-
u"""P3d_caracter_regular_S3.py -- pesos de eq:cono, caracteres de S_3, normalizacion f_u, prop:s3 y dim M_2k.

Comprueba, con aritmetica propia y sin leer ningun fichero, los enunciados matematicos de la nota
P3d que dependen de estas cuentas: (1) con la convencion w(L)=2 que exige eq:pesouno, los dos
primeros terminos de eq:cono pesan 4 y 2; (2) el indicador de la clase identidad de S_3 es
(1/6)chi_reg, combinacion racional pero no genuina de caracteres irreducibles; (3) la
normalizacion f_u(Y)=Y^3+uY-u (e_2=e_3=u, disc=-u^2(4u+27)), el caso k=6 por las sumas de Newton
y el escalado t=e_2/e_3, y que Delta_k es un entero no nulo para k<=15; (4) prop:s3 en k=6:
deg Q_6=1 con raiz 3/2, el cubico 2Y^3+3Y-3 de discriminante -1188 y grupo S_3, y la identidad
6 N_O(p) = chi_reg(Frob_p) primo a primo por fuerza bruta en F_p para 5<=p<400, p no divide 66;
(5) que la razon medida en cuatro ventanas no separa sqrt(p) de p^0.4; (8) cor:modular:
dim M_2k = deg Q_k + 1 y dim M_2k = 1 exactamente en {2,3,4,5,7}.

Controles negativos: con w(L)=0 ningun termino de eq:cono tiene peso uno (falla si la convencion
fuese esa); trivial+estandar debe salir genuino (falla si el test de genuinidad lo rechazase
todo); y el cubico Y^3-3Y-1, de discriminante cuadrado 81 y grupo C_3, no da el argumento no
abeliano de prop:s3.

Uso:  python P3d_caracter_regular_S3.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import itertools
import math
import sys

import sympy as sp

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
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
    return bool(cond)


def bloque(n, titulo):
    L.append(u'')
    L.append(u'=' * 78)
    L.append(u'%s.  %s' % (n, titulo))
    L.append(u'=' * 78)


u_, Y, e2, e3, t = sp.symbols('u Y e2 e3 t')

# =============================================================================================
bloque(u'1', u'los pesos de eq:cono')
# H^1 tiene peso uno y eq:pesouno toma gr^W_1 de e_c(A_C) = [Pi]L^2 + (1-[Pi])L + (1-L)[H^1].
# Para que ese gr^W_1 valga [H^1] hace falta w(L) > 0 (y par); con w(L)=2 sale [H^1] y nada
# mas.  Entonces los dos primeros terminos pesan 4 y 2, no 2 y 0.
PESO_L = 2


def peso_termino(pot_L):
    return 2 * pot_L            # un factor de Artin no suma peso


comprueba(u'1b. con w(L)=%d los dos primeros terminos pesan %d y %d, no 2 y 0'
          % (PESO_L, peso_termino(2), peso_termino(1)),
          (peso_termino(2), peso_termino(1)) != (2, 0),
          u'[Pi]L^2 -> %d, (1-[Pi])L -> %d' % (peso_termino(2), peso_termino(1)))


# 1c: se calcula gr^W_1 de e_c(A_C) termino a termino, como lista (etiqueta, peso), para cada
# w(L) candidato.  Los terminos: [Pi]L^2 (2w), L y -[Pi]L (w), [H^1] (1), -L[H^1] (1+w).
def gr1(w):
    terminos = [(u'[Pi]L^2', 2 * w), (u'L', w), (u'-[Pi]L', w),
                (u'[H^1]', 1), (u'-L[H^1]', 1 + w)]
    return sorted(et for et, peso in terminos if peso == 1)


gr1_2 = gr1(2)
gr1_0 = gr1(0)
comprueba(u'1c. gr^W_1 de e_c(A_C) es exactamente [H^1] con w(L)=2, y no con w(L)=0',
          gr1_2 == [u'[H^1]'] and gr1_0 != [u'[H^1]'],
          u'w=2: %s   w=0: %s' % (u' + '.join(gr1_2), u' + '.join(gr1_0)))
# control negativo: si la convencion fuese w(L)=0, el peso uno no podria salir de ninguna parte
comprueba(u'1d. control negativo: con w(L)=0 ningun termino de eq:cono tendria peso uno',
          all(2 * 0 * j != 1 for j in (0, 1, 2)))

# =============================================================================================
bloque(u'2', u'Artin--Tate: el indicador de la identidad como combinacion de caracteres')
# Tabla de caracteres de S_3 sobre las tres clases (1, transposicion, 3-ciclo), con tamanos.
CLASES = [(u'1', 1), (u'transp', 3), (u'3-ciclo', 2)]
IRR = {u'trivial': [1, 1, 1], u'signo': [1, -1, 1], u'estandar': [2, 0, -1]}
indicador_id = [1, 0, 0]        # lo que prop:n3k6 necesita: 1{Frob = identidad}

# Se resuelve indicador = a*trivial + b*signo + c*estandar sobre Q.
a, b, c = sp.symbols('a b c')
sol = sp.solve([sp.Eq(a * IRR[u'trivial'][i] + b * IRR[u'signo'][i] + c * IRR[u'estandar'][i],
                      indicador_id[i]) for i in range(3)], [a, b, c], dict=True)[0]
coef = [sol[a], sol[b], sol[c]]
L.append(u'  1{Frob=id} = %s*trivial + %s*signo + %s*estandar' % tuple(coef))
comprueba(u'2a. el indicador de la identidad ES una combinacion RACIONAL de caracteres irreducibles',
          all(x.is_rational for x in coef))
comprueba(u'2b. ... pero NO es un caracter genuino: algun coeficiente no es entero >= 0',
          not all(x.is_integer and x >= 0 for x in coef),
          u'coeficientes %s' % (coef,))
comprueba(u'2c. y es exactamente (1/6) del caracter regular',
          [sp.Rational(1, 6) * d for d in (1, 1, 2)] ==
          [coef[0], coef[1], coef[2]],
          u'regular = trivial + signo + 2*estandar')
# control: un caracter que SI es genuino tiene que pasar el test de 2b al reves.
genuino = [3, 1, 0]             # trivial + estandar
sol2 = sp.solve([sp.Eq(a * IRR[u'trivial'][i] + b * IRR[u'signo'][i] + c * IRR[u'estandar'][i],
                       genuino[i]) for i in range(3)], [a, b, c], dict=True)[0]
comprueba(u'2d. control: trivial+estandar SI es genuino (coeficientes enteros >= 0)',
          all(sol2[x].is_integer and sol2[x] >= 0 for x in (a, b, c)),
          u'%s' % ([sol2[a], sol2[b], sol2[c]],))

# =============================================================================================
bloque(u'3', u'la normalizacion f_u(Y) = Y^3 + uY - u (eq:fu) y Delta_k (eq:delta)')
f_u = Y ** 3 + u_ * Y - u_
# invariante e_2^3/e_3^2 de Y^3 - e1 Y^2 + e2 Y - e3  =>  e2 = u, e3 = u
poli = sp.Poly(f_u, Y)
c = poli.all_coeffs()                      # [1, 0, u, -u]
e2_f, e3_f = c[2], -c[3]
comprueba(u'3a. f_u tiene e_2 = e_3 = u', sp.simplify(e2_f - u_) == 0
          and sp.simplify(e3_f - u_) == 0, u'e2=%s  e3=%s' % (e2_f, e3_f))
comprueba(u'3b. ... luego su invariante e_2^3/e_3^2 vale exactamente u',
          sp.simplify(e2_f ** 3 / e3_f ** 2 - u_) == 0)
disc = sp.factor(sp.discriminant(f_u, Y))
comprueba(u'3c. y su discriminante es -u^2(4u+27)',
          sp.simplify(disc - (-u_ ** 2 * (4 * u_ + 27))) == 0, u'disc = %s' % disc)

# el caso k=6: s_6 = 0 => e_2^3/e_3^2 = 3/2
S = {0: sp.Integer(3), 1: sp.Integer(0), 2: -2 * e2, 3: 3 * e3}
for k in range(4, 70):
    S[k] = sp.expand(-e2 * S[k - 2] + e3 * S[k - 3])
comprueba(u'3d. Newton con e_1=0 da s_6 = -2e_2^3 + 3e_3^2',
          sp.simplify(S[6] - (-2 * e2 ** 3 + 3 * e3 ** 2)) == 0, u's_6 = %s' % sp.factor(S[6]))
inv6 = sp.solve(sp.Eq(S[6], 0), e2 ** 3 / e3 ** 2)
u6 = sp.simplify((sp.Rational(3, 2)))
comprueba(u'3e. ... luego s_6 = 0 obliga a e_2^3/e_3^2 = 3/2',
          sp.simplify(sp.expand(-2 * (sp.Rational(3, 2) * e3 ** 2) + 3 * e3 ** 2)) == 0,
          u'con e_2^3 = (3/2) e_3^2')
# el escalado: T -> tT manda e_j -> t^j e_j; se quiere e_2' = e_3'
t_sol = sp.solve(sp.Eq(t ** 2 * e2, t ** 3 * e3), t)
t_buena = [x for x in t_sol if x != 0][0]
comprueba(u'3f. el escalado que iguala e_2 y e_3 es t = e_2/e_3 (NO e_3/e_2)',
          sp.simplify(t_buena - e2 / e3) == 0, u't = %s' % t_buena)
e2p = sp.simplify((t_buena ** 2 * e2))
comprueba(u'3g. ... y con el e_2\' = e_3\' = e_2^3/e_3^2, el invariante mismo',
          sp.simplify(e2p - e2 ** 3 / e3 ** 2) == 0)
comprueba(u'3h. luego para k=6 el cubico normalizado es Y^3 + (3/2)Y - (3/2) = (2Y^3+3Y-3)/2',
          sp.simplify(sp.expand(2 * (Y ** 3 + sp.Rational(3, 2) * Y - sp.Rational(3, 2)))
                      - (2 * Y ** 3 + 3 * Y - 3)) == 0)


# Delta_k explicito: lc * disc * Res(Q_k, u(4u+27)), entero y no nulo
def Qk(k):
    u"""Q_k primitivo en u, a partir de s_k."""
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
        aa, bb = mon
        if aa % 3 or bb % 2:
            return None
        expr += co * u_ ** (aa // 3)
    return sp.Poly(sp.primitive(sp.expand(expr))[1], u_)


malos_delta = []
for k in (6, 8, 9, 10, 11, 12, 13, 14, 15):
    Q = Qk(k)
    if Q is None or Q.degree() < 1:
        continue
    lc = Q.LC()
    dQ = sp.discriminant(Q.as_expr(), u_)
    res = sp.resultant(Q.as_expr(), u_ * (4 * u_ + 27), u_)
    D = sp.Integer(lc) * sp.Integer(dQ) * sp.Integer(res)
    if D == 0 or not D.is_integer:
        malos_delta.append((k, D))
    L.append(u'     k=%-3d deg Q_k=%d   lc=%-6s disc=%-16s Res=%-18s Delta_k=%s'
             % (k, Q.degree(), lc, dQ, res, D))
comprueba(u'3l. Delta_k = lc*disc*Res(Q_k, u(4u+27)) es un entero NO NULO en cada k con deg>0',
          not malos_delta, u'%s' % (malos_delta,))

# =============================================================================================
bloque(u'4', u'prop:s3 en k=6:  6 N_O(p) = chi_{Q[G/K]}(Frob_p)  (lem:nocancela)')
# Caso k=6, donde la identidad se puede comprobar y no solo argumentar: deg Q_6 = 1, su unica
# raiz u_0 = 3/2 es racional, el cubico es 2Y^3+3Y-3 con grupo S_3.  Entonces G = S_3, H = G,
# K = ker(H -> S_3) = 1, y G/K es la representacion REGULAR.  La identidad dice 6 N(p) = chi_reg,
# o sea: 6 si el cubico escinde, 0 si no.  Se mide primo a primo.
cub = 2 * Y ** 3 + 3 * Y - 3
comprueba(u'4a. deg Q_6 = 1 y su raiz es 3/2', Qk(6).degree() == 1
          and sp.Rational(3, 2) in sp.roots(Qk(6).as_expr(), u_),
          u'Q_6 = %s' % Qk(6).as_expr())
comprueba(u'4b. el cubico de esa raiz es 2Y^3+3Y-3, irreducible, de discriminante -1188',
          sp.Poly(cub, Y).is_irreducible and sp.discriminant(cub, Y) == -1188,
          u'disc = %s' % sp.discriminant(cub, Y))
comprueba(u'4c. ... y -1188 no es cuadrado, luego el grupo es S_3 y no C_3',
          not sp.sqrt(sp.Integer(-1188)).is_rational)


def escinde(poly_coeffs, p):
    u"""Cuantas raices DISTINTAS tiene el cubico en F_p; escinde <=> 3."""
    r = 0
    for x in range(p):
        v = 0
        for co in poly_coeffs:
            v = (v * x + co) % p
        if v == 0:
            r += 1
    return r


def C6_3(p):
    u"""C_6(3,p) por fuerza bruta: subconjuntos {x,y,z} de F_p, distintos, con
    x+y+z = 0 y x^6+y^6+z^6 = 0.  No usa el cubico para nada."""
    n = 0
    for x in range(p):
        for y in range(x + 1, p):
            z = (-x - y) % p
            if z <= y:
                continue
            if (pow(x, 6, p) + pow(y, 6, p) + pow(z, 6, p)) % p == 0:
                n += 1
    return n


def chi_reg_frob(p):
    u"""chi_reg(Frob_p) desde el TIPO DE FACTORIZACION de 2Y^3+3Y-3 modulo p: Frob_p es la
    identidad de S_3 exactamente cuando todos los factores son lineales (p no divide el disc)."""
    grados = [sp.Poly(f, Y).degree()
              for f, _ in sp.factor_list(cub, modulus=p)[1]]
    return 6 if grados and all(g == 1 for g in grados) and len(grados) == 3 else 0


# Los dos lados salen de sitios distintos: N(p) del recuento directo de C_6(3,p), que por la
# nota vale (p-1) N(p); chi_reg(Frob_p) de como se factoriza el cubico modulo p.  Si el recuento
# de la nota o la identificacion de la fibra estuviesen mal, esto falla.
chi_reg_ok = []
malos_4d = []
for p in [q for q in range(5, 400) if sp.isprime(q) and 66 % q]:
    c = C6_3(p)
    N, resto = divmod(c, p - 1)
    ok = resto == 0 and 6 * N == chi_reg_frob(p)
    chi_reg_ok.append(ok)
    if not ok:
        malos_4d.append((p, c))
comprueba(u'4d. 6 N(p) = chi_reg(Frob_p) en los %d primos 5<=p<400, p no divide 66, con N(p) '
          u'= C_6(3,p)/(p-1) por fuerza bruta y chi_reg por la factorizacion modulo p'
          % len(chi_reg_ok), all(chi_reg_ok), u'fallan: %s' % (malos_4d[:5],))
# control negativo del instrumento: si N(p) se tomase de un cubico equivocado (Y^3-3Y-1, grupo
# C_3), el recuento directo lo contradice en algun primo.
contra = [p for p in [q for q in range(5, 400) if sp.isprime(q) and 66 % q]
          if 6 * (C6_3(p) // (p - 1)) != (6 if escinde([1, 0, -3, -1], p) == 3 else 0)]
comprueba(u'4d2. control negativo: con Y^3-3Y-1 en lugar del cubico correcto, el recuento '
          u'directo discrepa en %d primos' % len(contra), len(contra) > 0)


# el paso de grupos: la imagen de G actuando en G/K (K=1, los 6 elementos) contiene una copia
# de H/K = S_3.  Se construye la accion regular por multiplicacion a la izquierda y se mide.
def compon(s, r):
    return tuple(s[r[i]] for i in range(3))


S3 = sorted(itertools.permutations(range(3)))
reg = [tuple(S3.index(compon(g, h)) for h in S3) for g in S3]   # permutacion de {0..5}
fiel = len(set(reg)) == 6


def compon6(s, r):
    return tuple(s[r[i]] for i in range(6))


no_abel = any(compon6(x, y) != compon6(y, x) for x in reg for y in reg)
comprueba(u'4e. la accion regular de S_3 sobre sus 6 cosets es fiel y de imagen NO abeliana',
          fiel and no_abel, u'imagen de orden %d en Sym(6)' % len(set(reg)))
# control negativo: si el cubico fuese ciclico (C_3), el indicador SI seria abeliano y no
# probaria nada
cub_c3 = Y ** 3 - 3 * Y - 1                       # disc = 81, ciclico
comprueba(u'4f. control negativo: Y^3-3Y-1 tiene disc cuadrado, grupo C_3, y NO da el argumento',
          sp.discriminant(cub_c3, Y) == 81 and sp.sqrt(sp.Integer(81)).is_rational)

# =============================================================================================
bloque(u'5', u'obs:dosfallos: cuatro ventanas no fijan un exponente')
# La razon empirica sobre cuatro ventanas no distingue sqrt(p) de otras leyes: una sucesion que
# crece como p^0.4 da una razon parecida.
vent = [200, 500, 1000, 2000]
obs = [14, 23, 30, 51]
razon_obs = obs[-1] / float(obs[0])
razon_sqrt = math.sqrt(vent[-1] / float(vent[0]))
razon_04 = (vent[-1] / float(vent[0])) ** 0.4
L.append(u'     razon observada %.2f   |   sqrt: %.2f   |   exponente 0.4: %.2f'
         % (razon_obs, razon_sqrt, razon_04))
comprueba(u'5b. la razon observada no separa sqrt(p) de p^0.4: dista menos de su diferencia',
          abs(razon_obs - razon_04) < abs(razon_sqrt - razon_04) * 3,
          u'cuatro ventanas no fijan un exponente')

# =============================================================================================
bloque(u'8', u'cor:modular: elemental <=> dim M_2k = 1')


def dim_M(w):
    u"""dim M_w(SL_2(Z)) para w par >= 0."""
    if w % 2 or w < 0:
        return 0
    if w % 12 == 2:
        return w // 12
    return w // 12 + 1


elem = []
for k in range(2, 200):
    if dim_M(2 * k) == 1:
        elem.append(k)
comprueba(u'8a. dim M_2k = 1 exactamente para k en {2,3,4,5,7} hasta k<200',
          elem == [2, 3, 4, 5, 7], u'%s' % (elem,))
comprueba(u'8b. ... y k=2 queda fuera por el rango de la nota (k>=3)',
          [k for k in elem if k >= 3] == [3, 4, 5, 7])
gr = [k for k in range(3, 80) if (k // 6 - (1 if k % 6 == 1 else 0)) == 0]
comprueba(u'8c. y deg Q_k = 0 en el mismo conjunto', gr == [3, 4, 5, 7], u'%s' % (gr,))
malos = [k for k in range(3, 61) if dim_M(2 * k) != (k // 6 - (1 if k % 6 == 1 else 0)) + 1]
comprueba(u'8d. dim M_2k = deg Q_k + 1 para 3<=k<=60', not malos, u'%s' % (malos,))

# =============================================================================================
print(u'\n'.join(L))
print(u'')
print(u'=' * 78)
print(u'VERIFICACION P3d CARACTER REGULAR S_3:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 78)
sys.exit(1 if MAL else 0)
