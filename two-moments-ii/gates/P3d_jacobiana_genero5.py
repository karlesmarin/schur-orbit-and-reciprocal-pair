# -*- coding: utf-8 -*-
u"""P3d_jacobiana_genero5.py -- la traza en K_0 (x) E y la jacobiana de genero cinco de P3d.

Comprueba dos grupos de enunciados de la nota P3d.  (1) El lado Artin--Tate es una clase
VIRTUAL en K_0^{ss}(G_Q, Qlbar) (x) E y la traza la determina porque los caracteres
irreducibles son E-linealmente independientes: se mide el rango de las tablas de caracteres de
S_3 y de C_3 (esta sobre E = Q(zeta_3)), y que el indicador de {1} en S_3 tiene coeficientes no
enteros.  (2) La jacobiana de la curva D: v^2 = t^6+12t^4+40t^3+36t^2+48t+16 (apendice de la
jacobiana de genero cinco): por fuerza bruta se cuentan #D(F_p) y #D(F_{p^2}) en p = 11 y 19, de
ahi los polinomios de Frobenius; se certifica su irreducibilidad sobre Q (modulo 2 y modulo 17),
que el polinomio de Weil sobre F_{p^m} es irreducible para m <= 24 (sumas de potencias enteras,
exacto), y que los subcuerpos cuadraticos reales, leidos en los resolventes X^2-X-27 y
X^2+5X-20 de discriminantes 109 y 105, son distintos: eso es lo que cierra End = Z.

Controles negativos, y por que pueden fallar:
  1d. C_3 sobre Q: el polinomio caracteristico de la permutacion ciclica se parte en grados 1 y
      2, dos irreducibles para tres clases.  Si Q escindiera C_3 saldrian tres factores lineales
      y el control fallaria: muestra que el enunciado necesita el cuerpo de escision E.
  2f. y^2 = x^5-1, con multiplicacion compleja por Q(zeta_5): en p = 11 y 31 los dos subcuerpos
      reales COINCIDEN (Q(sqrt(5))), luego el criterio de subcuerpos distintos correctamente no
      certifica End = Z en una curva donde End != Z.  Si el criterio diera subcuerpos distintos
      tambien aqui, no discriminaria nada y el control fallaria.
  2i. modulo 3 el cuartico de p = 11 SI se parte: un certificado de irreducibilidad que valiera
      con cualquier primo no certificaria nada; el testigo es el primo concreto (2 y 17).
  2l-2n. restriccion de Weil T^4-T^2+121 en p = 11: irreducible sobre Q como los dos de D, pero
      su polinomio de Weil se parte ya en m = 2, como (X^2-X+121)^2.  Es simple y NO
      absolutamente simple; si la rutina de 2k no detectara particiones, 2m fallaria.

Uso:  python P3d_jacobiana_genero5.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

import numpy as np
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


def bloque(n, t):
    L.append(u'')
    L.append(u'=' * 78)
    L.append(u'BLOQUE %s.  %s' % (n, t))
    L.append(u'=' * 78)


# =============================================================================================
bloque(u'1', u'el lado Artin--Tate es VIRTUAL: el enunciado vive en K_0 (x) E')
# Los caracteres irreducibles de un grupo son linealmente independientes sobre cualquier cuerpo
# de escision de caracteristica cero: la matriz de caracteres es cuadrada e invertible.  Eso es
# exactamente lo que hace que la traza determine la clase en K_0^{ss} (x) E, aunque el elemento
# no sea la clase de ninguna representacion genuina.
TABLA_S3 = sp.Matrix(3, 3, [1, 1, 1,    # trivial   en clases  1, (12), (123)
                            1, -1, 1,   # signo
                            2, 0, -1])  # estandar
comprueba(u'1a. la tabla de caracteres de S_3 es cuadrada, 3 clases y 3 irreducibles',
          TABLA_S3.shape == (3, 3))
comprueba(u'1b. y tiene rango maximo, luego la traza determina la clase en K_0 (x) E',
          TABLA_S3.rank() == 3, u'det = %s' % TABLA_S3.det())
z = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2
TABLA_C3 = sp.Matrix(3, 3, [1, 1, 1, 1, z, z ** 2, 1, z ** 2, z ** 4])
comprueba(u'1c. y lo mismo para C_3 SOBRE E = Q(zeta_3), que lo escinde',
          sp.simplify(TABLA_C3.det()) != 0, u'tres irreducibles de dimension uno')
# CONTROL NEGATIVO: sobre Q la cosa falla, y por eso el enunciado pide un cuerpo de escision.
M = sp.Matrix(3, 3, [0, 0, 1, 1, 0, 0, 0, 1, 0])
grados_Q = sorted(sp.Poly(f, sp.Symbol('lambda')).degree()
                  for f, _ in sp.factor_list(M.charpoly().as_expr())[1])
comprueba(u'1d. CONTROL NEGATIVO: sobre Q solo hay DOS irreducibles de C_3, de grados 1 y 2',
          grados_Q == [1, 2], u'menos irreducibles que clases: la traza NO basta sobre Q')
# el elemento de K_0 (x) E es virtual de verdad: el indicador de {1} en S_3 tiene coeficientes
# fraccionarios en la base de caracteres irreducibles.
coef = sp.Matrix([sp.Rational(1, 6), sp.Rational(1, 6), sp.Rational(1, 3)])
ind = (TABLA_S3.T * coef).T
comprueba(u'1e. el indicador de {1} en S_3 es (1/6, 1/6, 1/3) en la base de irreducibles',
          list(ind) == [1, 0, 0], u'valores en las tres clases: %s' % list(ind))
comprueba(u'1f. ... con coeficientes NO enteros, luego no es ninguna representacion genuina',
          any(c.q != 1 for c in coef), u'1/6, 1/6, 1/3')

# =============================================================================================
bloque(u'2', u'la jacobiana de genero cinco: que fuerzan de verdad los dos primos')
COEF = [1, 0, 12, 40, 36, 48, 16]          # D: v^2 = t^6+12t^4+40t^3+36t^2+48t+16


def _pol(coef, t, p):
    y = 0
    for c in coef:
        y = (y * t + c) % p
    return y


def puntos1(coef, p, inf):
    Q = set((a * a) % p for a in range(1, p))
    n = inf
    for t in range(p):
        y = _pol(coef, t, p)
        if y == 0:
            n += 1
        elif y in Q:
            n += 2
    return n


def puntos2(coef, p, inf):
    nr = next(u for u in range(2, p) if sp.legendre_symbol(u, p) == -1)

    def mul(a, b):
        return ((a[0] * b[0] + nr * a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    Q = set()
    for a in range(p):
        for b in range(p):
            if (a, b) != (0, 0):
                Q.add(mul((a, b), (a, b)))
    n = inf
    for a in range(p):
        for b in range(p):
            y = (0, 0)
            for c in coef:
                y = mul(y, (a, b))
                y = ((y[0] + c) % p, y[1])
            if y == (0, 0):
                n += 1
            elif y in Q:
                n += 2
    return n


T = sp.Symbol('T')


def frobenius(coef, p, inf):
    """Devuelve (n1, n2, a, b, P, irreducible, sin_raices_de_unidad, disc_real)."""
    n1, n2 = puntos1(coef, p, inf), puntos2(coef, p, inf)
    s1, s2 = p + 1 - n1, p ** 2 + 1 - n2
    a = -s1
    b = (s1 * s1 - s2) // 2
    P = sp.expand(T ** 4 + a * T ** 3 + b * T ** 2 + p * a * T + p ** 2)
    irr = len(sp.factor_list(sp.Poly(P, T))[1]) == 1 \
        and sp.factor_list(sp.Poly(P, T))[1][0][1] == 1 \
        and sp.Poly(sp.factor_list(sp.Poly(P, T))[1][0][0], T).degree() == 4
    rs = np.roots([1, a, b, p * a, p * p])
    uni = False
    for i in range(4):
        for j in range(4):
            if i != j and any(abs((rs[i] / rs[j]) ** m - 1) < 1e-8 for m in range(1, 25)):
                uni = True
    # alpha + p/alpha satisface  X^2 + aX + (b - 2p):  el subcuerpo cuadratico REAL
    disc = sp.Integer(a * a - 4 * (b - 2 * p))
    return n1, n2, a, b, P, irr, (not uni), sp.factorint(disc)


L.append(u'      D:  v^2 = t^6+12t^4+40t^3+36t^2+48t+16   (lc cuadrado: 2 puntos en el infinito)')
datos = {}
for p in (11, 19):
    n1, n2, a, b, P, irr, sinuni, df = frobenius(COEF, p, 2)
    # parte libre de cuadrados del discriminante del subcuerpo real
    libre = 1
    for q, e in df.items():
        if e % 2:
            libre *= int(q)
    datos[p] = (n1, P, irr, sinuni, libre)
    L.append(u'      p=%-3d  #D(F_p)=%-4d #D(F_p^2)=%-5d  P(T) = %s' % (p, n1, n2, P))
    L.append(u'            irreducible: %s   sin cocientes raiz de la unidad: %s   '
             u'subcuerpo real: Q(sqrt(%d))' % (irr, sinuni, libre))

comprueba(u'2a. #D(F_11) = 11', datos[11][0] == 11)
comprueba(u'2b. #D(F_19) = 25', datos[19][0] == 25)
comprueba(u'2c. los dos polinomios de Frobenius son irreducibles sobre Q',
          datos[11][2] and datos[19][2])
comprueba(u'2d. y ningun cociente de raices es raiz de la unidad',
          datos[11][3] and datos[19][3], u'luego J(D) mod p es absolutamente simple en los dos')
comprueba(u'2e. los subcuerpos cuadraticos reales son DISTINTOS',
          datos[11][4] != datos[19][4],
          u'Q(sqrt(%d)) y Q(sqrt(%d)) --- eso es lo que cierra End = Z, y no la irreducibilidad'
          % (datos[11][4], datos[19][4]))
# CONTROL NEGATIVO: y^2 = x^5 - 1 tiene multiplicacion compleja por Q(zeta_5).  Todos sus
# cuerpos de Frobenius viven en Q(zeta_5), asi que los subcuerpos reales COINCIDEN y el criterio
# no puede certificar End = Z --- que es lo correcto, porque End no es Z.
CM = [1, 0, 0, 0, 0, -1]                  # x^5 - 1, grado impar: 1 punto en el infinito
pares = []
for p in (11, 31):                        # p = 1 mod 5, buena reduccion
    _, _, _, _, _, irr, sinuni, df = frobenius(CM, p, 1)
    libre = 1
    for q, e in df.items():
        if e % 2:
            libre *= int(q)
    pares.append(libre)
L.append(u'      CONTROL NEGATIVO y^2 = x^5-1 (CM por Q(zeta_5)): subcuerpos reales '
         u'Q(sqrt(%d)), Q(sqrt(%d))' % tuple(pares))
comprueba(u'2f. CONTROL NEGATIVO: con multiplicacion compleja los dos subcuerpos reales '
          u'COINCIDEN', pares[0] == pares[1],
          u'el criterio correctamente NO certifica End = Z en una curva donde End != Z')
# --- los certificados EXACTOS de los dos cuarticos, al mismo liston que la Tabla 2 ---------
TT = sp.Symbol('T')
CUART = {11: (-1, -5), 19: (5, 18)}


def cuartico(p):
    a, b = CUART[p]
    return TT ** 4 + a * TT ** 3 + b * TT ** 2 + p * a * TT + p ** 2


def irreducible_mod(expr, q, grado):
    fac = sp.factor_list(expr, modulus=q)[1]
    return len(fac) == 1 and fac[0][1] == 1 and sp.Poly(fac[0][0], TT).degree() == grado


comprueba(u'2g. el cuartico de p=11 es ya irreducible modulo 2',
          irreducible_mod(cuartico(11), 2, 4), u'el certificado, verificable en una linea')
comprueba(u'2h. el de p=19 lo es modulo 17', irreducible_mod(cuartico(19), 17, 4))
# CONTROL NEGATIVO: un primo cualquiera NO sirve; el certificado tiene que ser el primo dado.
comprueba(u'2i. CONTROL NEGATIVO: modulo 3 el de p=11 SI se parte, luego el testigo no es '
          u'cualquiera', not irreducible_mod(cuartico(11), 3, 4),
          u'un certificado que valiera con cualquier primo no certificaria nada')
# El polinomio de Weil sobre F_{p^m} es prod (X - alpha_i^m) --- NO es P(X^m), que es otro
# polinomio, de grado 4m.  Se calcula por sumas de potencias enteras: exacto y sin raices.
XX = sp.Symbol('X')


def _potencias(c, N):
    u"""p_k = suma alpha_i^k, con c = [a1,a2,a3,a4] de T^4 + a1 T^3 + a2 T^2 + a3 T + a4."""
    e = [1, -c[0], c[1], -c[2], c[3]]
    p = [4]
    for k in range(1, N + 1):
        v = 0
        if k <= 4:
            for i in range(1, k):
                v += (-1) ** (i + 1) * e[i] * p[k - i]
            v += (-1) ** (k + 1) * k * e[k]
        else:
            for i in range(1, 5):
                v += (-1) ** (i + 1) * e[i] * p[k - i]
        p.append(v)
    return p


def weil_sobre_extension(c, m):
    u"""Coeficientes de prod (X - alpha_i^m), grado 4, por Newton."""
    from fractions import Fraction
    p = _potencias(c, 4 * m)
    q = [p[j * m] for j in range(5)]
    E = [Fraction(1)]
    for k in range(1, 5):
        v = Fraction(0)
        for i in range(1, k + 1):
            v += (-1) ** (i - 1) * E[k - i] * q[i]
        E.append(v / k)
    assert all(x.denominator == 1 for x in E), E
    return [1, -int(E[1]), int(E[2]), -int(E[3]), int(E[4])]


def irr4(coefs):
    pol = sum(co * XX ** (4 - i) for i, co in enumerate(coefs))
    f = sp.factor_list(sp.expand(pol))[1]
    return len(f) == 1 and f[0][1] == 1 and sp.Poly(f[0][0], XX).degree() == 4


CUARTL = {11: [-1, -5, -11, 121], 19: [5, 18, 95, 361]}
comprueba(u'2j. la rutina reproduce P_1 = P en los dos primos',
          all(weil_sobre_extension(CUARTL[p], 1) == [1] + CUARTL[p] for p in (11, 19)),
          u'si P_1 no saliera, todo lo que sigue no mediria nada')
malos = {p: [m for m in range(1, 25) if not irr4(weil_sobre_extension(CUARTL[p], m))]
         for p in (11, 19)}
comprueba(u'2k. el polinomio de Weil sobre F_{p^m} es irreducible para todo m<=24, en los dos',
          not malos[11] and not malos[19],
          u'criterio de Tate, exacto por sumas de potencias enteras --- no por raices en coma '
          u'flotante, y NO es P(X^m)')
# CONTROL NEGATIVO: la restriccion de Weil a F_11 de una eliptica sobre F_{11^2}.  Su
# polinomio es Q(T^2), luego las raices vienen en pares +-alpha y colapsan al elevar: la
# superficie es simple sobre F_11 y NO absolutamente simple.  El criterio la caza en m=2.
NEG = [0, -1, 0, 121]                      # T^4 - T^2 + 121
comprueba(u'2l. CONTROL NEGATIVO: T^4-T^2+121 es irreducible sobre Q, como los dos de D',
          irr4([1] + NEG))
rompe = [m for m in range(1, 13) if not irr4(weil_sobre_extension(NEG, m))]
comprueba(u'2m. ... pero su polinomio de Weil SI se parte, ya en m=2',
          rompe[:1] == [2],
          u'se parte en m = %s --- el criterio distingue lo simple de lo absolutamente simple'
          % rompe)
comprueba(u'2n. ... y ahi es exactamente el cuadrado (X^2-X+121)^2',
          weil_sobre_extension(NEG, 2) == [1, -2, 243, -242, 14641]
          and sp.factor(XX ** 4 - 2 * XX ** 3 + 243 * XX ** 2 - 242 * XX + 14641)
          == (XX ** 2 - XX + 121) ** 2,
          u'la restriccion de Weil de una eliptica, que es el ejemplo clasico')
comprueba(u'2o. los resolventes X^2+a_1X+(a_2-2p) son X^2-X-27 y X^2+5X-20',
          (CUART[11][0], CUART[11][1] - 22) == (-1, -27)
          and (CUART[19][0], CUART[19][1] - 38) == (5, -20))
comprueba(u'2p. sus discriminantes 109 y 105 son libres de cuadrados',
          all(e == 1 for e in sp.factorint(109).values())
          and all(e == 1 for e in sp.factorint(105).values()),
          u'109 primo; 105 = 3*5*7')

print(u'\n'.join(L))
print(u'')
print(u'=' * 78)
print(u'VERIFICACION P3D JACOBIANA DE GENERO CINCO:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 78)
sys.exit(1 if MAL else 0)
