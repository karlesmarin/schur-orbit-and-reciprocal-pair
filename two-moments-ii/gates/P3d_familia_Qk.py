# -*- coding: utf-8 -*-
u"""P3d_familia_Qk.py -- la familia Q_k: Waring, recursiones, grado modular, cubico y curva D.

Comprueba, en aritmetica simbolica exacta (sympy) y por fuerza bruta sobre F_p, los enunciados
de la nota P3d sobre la familia Q_k: la formula de Newton-Waring eq:waring frente a la
recursion de Newton (2 <= k <= 60); prop:grado (deg Q_k = floor(k/2) - ceil(k/3)); la perdida
de monicidad en k = 35; las dos recursiones de R_k y las funciones simetricas de
beta_i = y_i^6/u^2 calculadas desde s_6 y s_12 sin extraer raices; prop:modular
(dim M_2k = deg Q_k + 1, y dim M_2k = 1 exactamente en k = 3,4,5,7); disc(Y^3+uY-u) =
-u^2(4u+27) y el caso k = 11 (Y^3+Y-1, disc -31, grupo S_3); la cota A < 2^d con sus casos de
igualdad; y la curva de genero cinco: j(E_0) = -235298/13, el discriminante del sextico de D,
#D(F_11) = 11, #D(F_19) = 25 y los polinomios de Weil W_11, W_19 con su ecuacion funcional.

Controles negativos: la cota A < 2^d separa el fallo estricto (k = 9,15,21,33,35) de la
igualdad (k = 6,10), y un control que confunda "<" con "<=" daria MAL; Q_35 debe salir NO
monico y primitivo, rompiendo el patron de los k primos 11,17,23,29; cada W_p debe predecir
el recuento de puntos de D medido por fuerza bruta, que es independiente de W_p.

Uso:  python P3d_familia_Qk.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys
from fractions import Fraction
from math import factorial, gcd

import sympy as sp

e2, e3, u, Y, Z, t, x, v = sp.symbols('e2 e3 u Y Z t x v')

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


def sec(x_):
    L.append(u'')
    L.append(x_)
    L.append(u'-' * len(x_))


KMAX = 60
S = {0: sp.Integer(3), 1: sp.Integer(0), 2: -2 * e2, 3: 3 * e3}
for k in range(4, KMAX + 1):
    S[k] = sp.expand(-e2 * S[k - 2] + e3 * S[k - 3])


def Qk(k):
    resto = sp.Integer(1)
    for base, m in sp.factor_list(sp.factor(S[k]))[1]:
        b = sp.expand(base)
        if b not in (e2, e3):
            resto *= b ** m
    if resto == 1:
        return sp.Poly(sp.Integer(1), u)
    P = sp.Poly(sp.expand(resto), e2, e3)
    expr = sp.Integer(0)
    for (a, b), c in zip(P.monoms(), P.coeffs()):
        assert a % 3 == 0 and b % 2 == 0
        expr += c * u ** (a // 3)
    return sp.Poly(sp.expand(expr), u)


# ---------------------------------------------------------------------------
sec(u'1.  La formula de Newton-Waring para s_k (eq:waring)')
# ---------------------------------------------------------------------------
def waring(k):
    tot = sp.Integer(0)
    for b in range(k // 3 + 1):
        if (k - 3 * b) % 2:
            continue
        a = (k - 3 * b) // 2
        if a + b == 0:
            continue
        c = Fraction((-1) ** a * factorial(a + b - 1), factorial(a) * factorial(b))
        tot += k * sp.Rational(c.numerator, c.denominator) * e2 ** a * e3 ** b
    return sp.expand(tot)


malos = [k for k in range(2, KMAX + 1) if sp.expand(waring(k) - S[k]) != 0]
comprueba(u's_k = k * suma (-1)^a (a+b-1)!/(a!b!) e2^a e3^b, 2 <= k <= %d' % KMAX, not malos,
      u'fallan: %s' % (malos[:5],))
# todos los coeficientes de eq:waring son distintos de cero
sin_cero = True
for k in range(4, 31):
    P = sp.Poly(S[k], e2, e3)
    if any(c == 0 for c in P.coeffs()):
        sin_cero = False
comprueba(u'ningun coeficiente de esa suma se anula', sin_cero)

# ---------------------------------------------------------------------------
sec(u'2.  El grado de Q_k (prop:grado)')
# ---------------------------------------------------------------------------
malos = []
for k in range(3, KMAX + 1):
    d1 = Qk(k).degree()
    d2 = k // 2 - -(-k // 3)                       # floor(k/2) - ceil(k/3)
    d3 = k // 6 - (1 if k % 6 == 1 else 0)
    if not (d1 == d2 == d3):
        malos.append((k, d1, d2, d3))
comprueba(u'deg Q_k = floor(k/2) - ceil(k/3) = floor(k/6) - [k=1 mod 6]', not malos,
      u'fallan: %s' % (malos[:4],))

# ---------------------------------------------------------------------------
sec(u'3.  La monicidad se rompe en k = 35, y por que')
# ---------------------------------------------------------------------------
q35 = Qk(35)
esperado = sp.Poly(35 * u ** 5 - 1225 * u ** 4 + 7007 * u ** 3
                   - 8580 * u ** 2 + 1925 * u - 35, u)
comprueba(u'Q_35 = 35u^5 - 1225u^4 + 7007u^3 - 8580u^2 + 1925u - 35',
      (q35 - esperado).is_zero or sp.expand(q35.as_expr() + esperado.as_expr()) == 0
      or sp.expand(q35.as_expr() - esperado.as_expr()) == 0,
      u'obtenido %s' % q35.as_expr())
comprueba(u'Q_35 es primitivo', gcd(gcd(*[abs(int(c)) for c in q35.all_coeffs()[:3]]),
                                gcd(*[abs(int(c)) for c in q35.all_coeffs()[3:]])) == 1)
comprueba(u'Q_35 NO es monico', abs(int(q35.LC())) != 1, u'lider = %d' % q35.LC())
# 11,17,23,29 primos; 35 no
primos_clase = [k for k in (11, 17, 23, 29, 35)
                if all(k % d for d in range(2, int(k ** 0.5) + 1))]
comprueba(u'11,17,23,29 son primos y 35 no', primos_clase == [11, 17, 23, 29])
monicos = [k for k in (11, 17, 23, 29, 35) if abs(int(Qk(k).LC())) == 1]
comprueba(u'y son exactamente los monicos de la clase k = 5 (mod 6)',
      monicos == [11, 17, 23, 29], u'monicos: %s' % monicos)
# para k primo, k divide a todos los coeficientes de s_k
mal = []
for k in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
    P = sp.Poly(S[k], e2, e3)
    if any(int(c) % k for c in P.coeffs()):
        mal.append(k)
comprueba(u'para k primo, k divide todo coeficiente de s_k', not mal, u'fallan: %s' % mal)

# ---------------------------------------------------------------------------
sec(u'4.  Las dos recurrencias')
# ---------------------------------------------------------------------------
# R_k(u) = u^{-ceil(k/3)} s_k(u,u)
def R(k):
    return sp.simplify(sp.expand(S[k].subs({e2: u, e3: u})) / u ** (-(-k // 3)))


malos = []
for k in range(3, 40):
    izq = sp.simplify(R(k))
    der = sp.simplify(R(k - 3) - (u if k % 3 == 0 else 1) * R(k - 2))
    if sp.simplify(izq - der) != 0:
        malos.append(k)
comprueba(u'R_k = R_{k-3} - u^[3|k] R_{k-2}, 3 <= k < 40', not malos,
      u'fallan: %s' % (malos[:5],))

malos = []
for k in range(0, 22):
    izq = sp.expand(R(k + 18))
    der = sp.expand((3 - 2 * u) * R(k + 12) - (u ** 2 + 6 * u + 3) * R(k + 6) + R(k))
    if sp.simplify(izq - der) != 0:
        malos.append(k)
comprueba(u'R_{k+18} = (3-2u)R_{k+12} - (u^2+6u+3)R_{k+6} + R_k', not malos,
      u'fallan: %s' % (malos[:5],))

# Y SU RAZON: las sumas simetricas de beta_i = y_i^6 / u^2, con y_i las raices de
# Y^3 + uY - u.  No hacen falta las raices --- sympy ni siquiera las ordena sobre ZZ[u]:
# basta con las sumas de potencias, que ya las tenemos.  Con e2 = e3 = u:
#    p_1(beta) = p_6(y)/u^2 ,  p_2(beta) = p_12(y)/u^4 ,  e_3(beta) = (prod y_i)^6/u^6 = 1
# y e_2(beta) = (p_1^2 - p_2)/2 por Newton.
p6 = sp.expand(S[6].subs({e2: u, e3: u}))
p12 = sp.expand(S[12].subs({e2: u, e3: u}))
b1 = sp.simplify(p6 / u ** 2)
b2 = sp.simplify((b1 ** 2 - p12 / u ** 4) / 2)
b3 = sp.simplify(u ** 6 / u ** 6)
comprueba(u'suma de beta_i = 3 - 2u', sp.simplify(b1 - (3 - 2 * u)) == 0, u'%s' % b1)
comprueba(u'suma de beta_i beta_j = u^2 + 6u + 3',
      sp.simplify(b2 - (u ** 2 + 6 * u + 3)) == 0, u'%s' % sp.expand(b2))
comprueba(u'producto de beta_i = 1', sp.simplify(b3 - 1) == 0)

# ---------------------------------------------------------------------------
sec(u'5.  El cubico normalizado y su discriminante')
# ---------------------------------------------------------------------------
fu = Y ** 3 + u * Y - u
comprueba(u'disc(Y^3 + uY - u) = -u^2 (4u + 27)',
      sp.simplify(sp.discriminant(fu, Y) - (-u ** 2 * (4 * u + 27))) == 0,
      u'%s' % sp.factor(sp.discriminant(fu, Y)))
comprueba(u'para k = 11, u = 1 y el cubico es Y^3 + Y - 1',
      Qk(11).as_expr() in (u - 1, -(u - 1)) and sp.simplify(fu.subs(u, 1) - (Y ** 3 + Y - 1)) == 0)
comprueba(u'Y^3 + Y - 1 es irreducible, disc = -31 no es cuadrado, grupo S_3',
      sp.Poly(Y ** 3 + Y - 1, Y).is_irreducible
      and sp.discriminant(Y ** 3 + Y - 1, Y) == -31)

# ---------------------------------------------------------------------------
sec(u'6.  La interpretacion modular: dim M_2k = d_k + 1 (prop:modular)')
# ---------------------------------------------------------------------------
def dim_M(w):
    u"""dim M_w(SL_2(Z)) para w par >= 0."""
    if w % 2 or w < 0:
        return 0
    if w % 12 == 2:
        return w // 12
    return w // 12 + 1


malos = [k for k in range(3, KMAX + 1) if dim_M(2 * k) != Qk(k).degree() + 1]
comprueba(u'dim M_{2k} = deg Q_k + 1, 3 <= k <= %d' % KMAX, not malos,
      u'fallan: %s' % (malos[:5],))
unos = [k for k in range(3, 200) if dim_M(2 * k) == 1]
comprueba(u'dim M_{2k} = 1 exactamente para k en {3,4,5,7}', unos == [3, 4, 5, 7],
      u'%s' % unos)
comprueba(u'los pesos son 6, 8, 10, 14 --- E_6, E_8, E_10, E_14',
      [2 * k for k in unos] == [6, 8, 10, 14])
comprueba(u'y con k = 2, fuera del rango de la nota, dim M_{2k} = 1 exactamente en {2,3,4,5,7}, k<200',
      [k for k in range(2, 200) if dim_M(2 * k) == 1] == [2, 3, 4, 5, 7])

# Segunda ruta para la dimension, sin la formula cerrada: M_* = C[E_4, E_6], asi que dim M_w es
# el numero de monomios E_4^a E_6^b con 4a + 6b = w.
mono = lambda w: sum(1 for a in range(w // 4 + 1) for b in range(w // 6 + 1) if 4 * a + 6 * b == w)
comprueba(u'dim M_w por monomios E_4^a E_6^b = formula cerrada, w par <= 400',
          all(mono(w) == dim_M(w) for w in range(0, 401, 2)))
comprueba(u'dim M_12 = 2 y dim M_14 = 1 (ap:delta), por monomios', mono(12) == 2 and mono(14) == 1)
comprueba(u'S_14 = Delta M_2 = 0: dim M_2 = 0', mono(2) == 0)

# Series de q exactas: E_4 = 1 + 240 sum sigma_3(n) q^n, E_6 = 1 - 504 sum sigma_5(n) q^n.
NQ = 30
sig = lambda n, r: sum(d ** r for d in range(1, n + 1) if n % d == 0)
E4 = [1] + [240 * sig(n, 3) for n in range(1, NQ)]
E6 = [1] + [-504 * sig(n, 5) for n in range(1, NQ)]


def mul(a, b):
    return [sum(a[i] * b[n - i] for i in range(n + 1)) for n in range(NQ)]


E4c, E6c = mul(mul(E4, E4), E4), mul(E6, E6)
disc1728 = [E4c[n] - E6c[n] for n in range(NQ)]
prod = [0] * NQ
prod[1] = 1                                   # q * prod (1 - q^n)^24
for n in range(1, NQ):
    for _ in range(24):
        prod = [prod[i] - (prod[i - n] if i >= n else 0) for i in range(NQ)]
comprueba(u'Delta = q prod(1-q^n)^24 = (E_4^3 - E_6^2)/1728 hasta q^%d' % (NQ - 1),
          all(disc1728[n] == 1728 * prod[n] for n in range(NQ)),
          u'tau(2..5) = %s' % (prod[2:6],))
comprueba(u'tau(2) = -24, tau(3) = 252', prod[2] == -24 and prod[3] == 252)
# s_6 con e_1 = 0 es -2 e_2^3 + 3 e_3^2; en e_2 = -3E_4, e_3 = 2E_6 da 54 E_4^3 + 12 E_6^2.
s6 = sp.expand(-2 * (-3 * x) ** 3 + 3 * (2 * v) ** 2)          # x = E_4, v = E_6
comprueba(u'calF_6/66 = (9E_4^3 + 2E_6^2)/11', sp.simplify(s6 / 66 - (9 * x ** 3 + 2 * v ** 2) / 11) == 0,
          u'calF_6 = %s' % s6)

# ---------------------------------------------------------------------------
sec(u'7.  La cota A < 2^d y los cinco casos que deja')
# ---------------------------------------------------------------------------
# OJO: la afirmacion es A < 2^d, y hay que separar el fallo ESTRICTO de la IGUALDAD.
# Los casos de igualdad son k = 6 y k = 10 y se tratan aparte, porque la desigualdad
# sigue siendo contradictoria salvo que todas las raices valgan 1/2, y u = 3/2, 15/2 no
# lo cumplen.  Un control que mezcle "<" con "<=" declara un fallo donde no lo hay.
estrictos, iguales = [], []
for k in range(3, 60):
    d = Qk(k).degree()
    if d == 0:
        continue
    A = abs(int(Qk(k).LC()))
    if A > 2 ** d:
        estrictos.append(k)
    elif A == 2 ** d:
        iguales.append(k)
comprueba(u'A > 2^d solo en k = 9, 15, 21, 33, 35', estrictos == [9, 15, 21, 33, 35],
      u'%s' % estrictos)
comprueba(u'los casos de igualdad A = 2^d son exactamente k = 6 y k = 10',
      iguales == [6, 10], u'%s' % iguales)
comprueba(u'y sus raices no son 1/2: u = 3/2 y u = 15/2',
      [sp.solve(Qk(k).as_expr(), u) for k in (6, 10)]
      == [ [sp.Rational(3, 2)], [sp.Rational(15, 2)] ])
comprueba(u'para k par, A <= 2', all(abs(int(Qk(k).LC())) <= 2
                                for k in range(4, 60, 2) if Qk(k).degree() > 0))
comprueba(u'para k primo impar >= 5, A = 1',
      all(abs(int(Qk(k).LC())) == 1 for k in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43)))

# ---------------------------------------------------------------------------
sec(u'8.  La curva de genero cinco: E_0 y la curva de genero 2 D')
# ---------------------------------------------------------------------------
E = x ** 3 - 147 * x - 718
c4 = -48 * (-147)
disc_E = -16 * (4 * (-147) ** 3 + 27 * 718 ** 2)
jE = sp.Rational(c4 ** 3, disc_E) if disc_E else None
comprueba(u'j(v^2 = x^3 - 147x - 718) = -235298/13', jE == sp.Rational(-235298, 13),
      u'j = %s' % jE)
sext = t ** 6 + 12 * t ** 4 + 40 * t ** 3 + 36 * t ** 2 + 48 * t + 16
d6 = sp.discriminant(sext, t)
comprueba(u'disc del sextico de D = 2^35 * 3^6 * 13', int(d6) == 2 ** 35 * 3 ** 6 * 13,
      u'disc = %s' % sp.factorint(int(d6)))
comprueba(u'D es lisa de genero 2 (sextico sin raiz doble)', d6 != 0)


def cuenta_D(p):
    u"""#D(F_p) para v^2 = sextico, modelo afin mas los puntos del infinito."""
    n = 0
    for tt in range(p):
        r = int(sext.subs(t, tt)) % p
        if r == 0:
            n += 1
        elif pow(r, (p - 1) // 2, p) == 1:
            n += 2
    # grado 6: dos puntos en el infinito si el coeficiente lider es cuadrado
    n += 2 if pow(1, (p - 1) // 2, p) == 1 else 0
    return n


for p, esp in ((11, 11), (19, 25)):
    comprueba(u'#D(F_%d) = %d' % (p, esp), cuenta_D(p) == esp, u'obtenido %d' % cuenta_D(p))

for p, coef in ((11, [121, -11, -5, -1, 1]), (19, [361, 95, 18, 5, 1])):
    W = sum(c * Z ** i for i, c in enumerate(coef))
    comprueba(u'W_%d es irreducible sobre Q' % p, sp.Poly(W, Z).is_irreducible)
    # p^2 = producto de las raices, y W(1/ (pX)) simetria funcional
    comprueba(u'W_%d tiene termino independiente p^2 = %d' % (p, p * p), coef[0] == p * p)
    # W = X^4 + a1 X^3 + a2 X^2 + p a1 X + p^2, y #C(F_p) = p + 1 + a1.
    # (Con P(X) = prod (X - alpha_i), la suma de raices es -a1 y #C = p+1-sum alpha.)
    a1 = coef[3]
    comprueba(u'W_%d cumple la ecuacion funcional: coef de X es p*a1' % p,
          coef[1] == p * a1, u'%d vs %d' % (coef[1], p * a1))
    N1 = p + 1 + a1
    comprueba(u'W_%d predice #D(F_%d) = %d' % (p, p, N1), N1 == cuenta_D(p),
          u'p+1+a1 = %d, medido %d' % (N1, cuenta_D(p)))

print(u'\n'.join(L))
print(u'')
print(u'=' * 70)
print(u'VERIFICACION FAMILIA Q_k:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 70)
sys.exit(1 if MAL else 0)
