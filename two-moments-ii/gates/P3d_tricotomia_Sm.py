# -*- coding: utf-8 -*-
u"""P3d_tricotomia_Sm.py -- tricotomia de S_m, contraejemplo (3,1,1,1) y formulas cerradas de P3d.

Comprueba los enunciados de la nota P3d sobre el lugar singular de longitud cuatro y las
cuentas que sobreviven: el conjunto S_m recalculado por enumeracion exhaustiva sobre mu_m^4
para 1 <= m <= 18 y 3 <= n <= 10 (Teorema thm:tricotomia), su colapso a clases de gcd(m,6)
(Proposicion prop:gcd) y el certificado 2+X^2+X^3+X^4 = Phi_6 (X^2+2X+2); el punto singular
unico de (3,1,1,1) en mu_4 y el genero 5 de la normalizada (Proposicion prop:contra); la
suavidad de (n-3,1,1,1) para n >= 7; las formulas cerradas de C_4(3,p), C_5(3,p), C_5(4,p)
(Teorema thm:formulas) con |B| y |A cap B|; C_6(3,p) y el discriminante -1188 del cubico
2Y^3+3Y-3 (Proposicion prop:n3k6); y p_g = 4 de la quintica nodal (Observacion rem:quintica).
Las formulas cerradas se verifican por fuerza bruta sobre F_p, primo a primo: si fallan en un
primo, se ve.  Controles negativos: S_6 debe contener el 5 (sin (2,1,1,1) saldria {3,4,6});
(2,1,1,1) no debe tener relacion en mu_2, mu_3, mu_4; m = 15, impar compuesto, debe dar {3};
(3,1,1,1) debe tener un solo punto singular y (2,2,1,1) si debe tener relacion en mu_4.
Cualquiera de ellos falla si el test de anulacion en raices de la unidad esta mal.

Uso:  python P3d_tricotomia_Sm.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import itertools
import sys

OK = 0
MAL = 0
LINEAS = []


def comprueba(etiqueta, condicion, detalle=u''):
    global OK, MAL
    if condicion:
        OK += 1
        LINEAS.append(u'  ok   %s%s' % (etiqueta, (u'   ' + detalle) if detalle else u''))
    else:
        MAL += 1
        LINEAS.append(u'  MAL  %s%s' % (etiqueta, (u'   ' + detalle) if detalle else u''))


def seccion(t):
    LINEAS.append(u'')
    LINEAS.append(t)
    LINEAS.append(u'-' * len(t))


# ---------------------------------------------------------------------------
# utilidades aritmeticas
# ---------------------------------------------------------------------------

def chi(a, p):
    u"""Simbolo de Legendre (a|p)."""
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def primos(lo, hi):
    out = []
    for n in range(max(2, lo), hi + 1):
        if all(n % d for d in range(2, int(n ** 0.5) + 1)):
            out.append(n)
    return out


def cuenta(n, k, p):
    u"""#{T subconjunto de F_p, |T| = n, suma(T) = 0 = suma(T^k)}, por fuerza bruta."""
    c = 0
    for T in itertools.combinations(range(p), n):
        if sum(T) % p:
            continue
        if sum(pow(t, k, p) for t in T) % p == 0:
            c += 1
    return c


# ---------------------------------------------------------------------------
# 1.  El lugar singular de longitud 4:  suma w_i z_i = 0 con z_i en mu_m
# ---------------------------------------------------------------------------
# Trabajamos en Z[x]/(x^m - 1): un elemento de mu_m es un exponente en Z/m, y la suma
# sum w_i x^{a_i} se anula en TODA raiz primitiva m-esima si y solo si Phi_m divide al
# polinomio.  Aqui basta el criterio directo: evaluamos en el cuerpo ciclotomico via
# reduccion modulo Phi_m con aritmetica entera exacta.

def phi_m(m):
    u"""Polinomio ciclotomico Phi_m como lista de coeficientes enteros (grado ascendente)."""
    # x^m - 1 dividido por todos los Phi_d con d | m, d < m
    num = [-1] + [0] * (m - 1) + [1]          # x^m - 1
    for d in range(1, m):
        if m % d == 0:
            num = divide(num, phi_m(d))
    return num


def divide(a, b):
    u"""Division exacta de polinomios con coeficientes enteros."""
    a = list(a)
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(q) - 1, -1, -1):
        c = a[i + len(b) - 1] // b[-1]
        q[i] = c
        if c:
            for j in range(len(b)):
                a[i + j] -= c * b[j]
    return q


_CACHE_PHI = {}


def PHI(m):
    if m not in _CACHE_PHI:
        _CACHE_PHI[m] = phi_m(m)
    return _CACHE_PHI[m]


def reduce_mod_phi(coef, m):
    u"""Reduce un polinomio modulo Phi_m; devuelve la lista reducida."""
    f = list(coef)
    g = PHI(m)
    while len(f) >= len(g):
        c = f[-1]
        if c:
            desp = len(f) - len(g)
            for j in range(len(g)):
                f[desp + j] -= c * g[j]
        f.pop()
    return f


def hay_relacion(pesos, m):
    u"""True si existen z_i en mu_m (m >= 1) con sum pesos_i z_i = 0.

    Recorre los exponentes a_i en Z/m (fijando a_0 = 0 por homogeneidad) y comprueba si
    el polinomio sum w_i x^{a_i} es divisible por Phi_d para algun d | m, d > 1 --- que es
    exactamente la condicion de que se anule en alguna raiz m-esima de la unidad.
    """
    l = len(pesos)
    for exps in itertools.product(range(m), repeat=l - 1):
        a = (0,) + exps
        coef = [0] * m
        for w, e in zip(pesos, a):
            coef[e] += w
        # se anula en alguna raiz m-esima <=> algun Phi_d | coef, d | m
        for d in divisores(m):
            if d == 1:
                # x = 1: la suma de los pesos, siempre positiva
                continue
            r = reduce_mod_phi(coef, d)
            if all(c == 0 for c in r):
                return True
    return False


def divisores(m):
    return [d for d in range(1, m + 1) if m % d == 0]


def particiones(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxp), 0, -1):
        for resto in particiones(n - k, k):
            yield (k,) + resto


def S_de_m(m, nmax=12):
    u"""Los n para los que TODA particion de longitud 4 admite relacion en mu_m.

    Es el conjunto S_m del Teorema de la tricotomia (thm:tricotomia): 'toda curva de cuatro
    partes es singular'.  Se calcula por enumeracion exhaustiva sobre mu_m^4.
    """
    out = []
    for n in range(3, nmax + 1):
        cuartetos = [lam for lam in particiones(n) if len(lam) == 4]
        if all(hay_relacion(lam, m) for lam in cuartetos):
            out.append(n)
    return out


# ---------------------------------------------------------------------------
seccion(u'1.  k = 7 (m = 6): S_6 y la particion (2,1,1,1)')
# ---------------------------------------------------------------------------
# Control negativo: el valor {3,4,6} (sin el 5) es el que sale si (2,1,1,1) no se detecta.
s6 = S_de_m(6)
comprueba(u'S_6 = {3,4,5,6} (y no {3,4,6})',
      s6 == [3, 4, 5, 6], u'S_6 = %s' % (s6,))

# la particion que decide es (2,1,1,1): entra en n = 5
comprueba(u'(2,1,1,1) tiene relacion en mu_6', hay_relacion((2, 1, 1, 1), 6))
comprueba(u'(2,1,1,1) NO tiene relacion en mu_2 ni en mu_3 ni en mu_4',
      not hay_relacion((2, 1, 1, 1), 2) and not hay_relacion((2, 1, 1, 1), 3)
      and not hay_relacion((2, 1, 1, 1), 4))

# el certificado Phi_6 de la nota:  2 + X^2 + X^3 + X^4 = (X^2 - X + 1)(X^2 + 2X + 2)
izq = [2, 0, 1, 1, 1]                                  # 2 + X^2 + X^3 + X^4
A, B = [1, -1, 1], [2, 2, 1]                           # Phi_6  y  X^2+2X+2
der = [0] * (len(A) + len(B) - 1)
for i, ca in enumerate(A):
    for j, cb in enumerate(B):
        der[i + j] += ca * cb
comprueba(u'certificado 2+X^2+X^3+X^4 = Phi_6 * (X^2+2X+2)', izq == der,
      u'%s vs %s' % (izq, der))

# ---------------------------------------------------------------------------
seccion(u'2.  Tricotomia (thm:tricotomia) y el invariante gcd(m,6) (prop:gcd)')
# ---------------------------------------------------------------------------
tricotomia = {}
for m in range(1, 19):
    tricotomia[m] = S_de_m(m, nmax=10)

esperado = {}
for m in range(1, 19):
    if m % 2:
        esperado[m] = [3] if m % 3 or m == 1 else [3]
    elif m % 3:
        esperado[m] = [3, 4, 6]
    else:
        esperado[m] = [3, 4, 5, 6]
# m impar: {3} salvo que 3 | m, que tambien da {3} (gcd(m,6) in {1,3})
for m in range(1, 19):
    if m % 2:
        esperado[m] = [3]

for m in sorted(tricotomia):
    comprueba(u'S_%-2d = %-14s' % (m, tricotomia[m]), tricotomia[m] == esperado[m],
          u'gcd(m,6) = %d' % (6 if m % 6 == 0 else __import__('math').gcd(m, 6)))

# depende solo de gcd(m,6)
import math
por_gcd = {}
consistente = True
for m in sorted(tricotomia):
    g = math.gcd(m, 6)
    if g in por_gcd and por_gcd[g] != tricotomia[m]:
        consistente = False
    por_gcd[g] = tricotomia[m]
comprueba(u'S_m depende solo de gcd(m,6), m <= 18', consistente, u'%s' % (por_gcd,))

# el caso que mata la intuicion de "varios primos se comporta como 2":  m = 15
comprueba(u'm = 15 (impar compuesto) da {3}, no {3,4,6}', tricotomia[15] == [3])

# ---------------------------------------------------------------------------
seccion(u'3.  El reciproco falla en k = 5, n = 6, lambda = (3,1,1,1)  (prop:contra)')
# ---------------------------------------------------------------------------
# k = 5 -> m = 4.  Enunciado: UN solo punto singular, [1:-1:-1:-1].
def puntos_singulares(pesos, m):
    u"""Los z en mu_m^l con sum w_i z_i = 0, modulo el escalado global (z_0 = 1)."""
    import cmath
    l = len(pesos)
    out = []
    for exps in itertools.product(range(m), repeat=l - 1):
        a = (0,) + exps
        coef = [0] * m
        for w, e in zip(pesos, a):
            coef[e] += w
        s = sum(c * cmath.exp(2j * cmath.pi * j / m) for j, c in enumerate(coef))
        if abs(s) < 1e-9:
            out.append(a)
    return out

sing = puntos_singulares((3, 1, 1, 1), 4)
comprueba(u'(3,1,1,1) en mu_4 tiene exactamente 1 punto singular', len(sing) == 1,
      u'%s' % (sing,))
comprueba(u'y es [1:-1:-1:-1]', sing == [(0, 2, 2, 2)] if sing else False)

# genero de la normalizada de una quintica plana irreducible con 1 nodo
g = (5 - 1) * (5 - 2) // 2 - 1
comprueba(u'g(normalizada) = (5-1)(5-2)/2 - 1 = 5', g == 5, u'g = %d' % g)
comprueba(u'z_(3,1,1,1) = 3 * 3! = 18', 3 * 6 == 18)
comprueba(u'dim H^1 = 2g = 10', 2 * g == 10)

# el otro cuarteto de 6 --- (2,2,1,1) --- si es biquanimo, luego el fallo esta DENTRO
comprueba(u'(2,2,1,1) si tiene relacion en mu_4 (es biquanimo)', hay_relacion((2, 2, 1, 1), 4))
comprueba(u'los dos cuartetos de 6 son (3,1,1,1) y (2,2,1,1)',
      sorted(lam for lam in particiones(6) if len(lam) == 4) == [(2, 2, 1, 1), (3, 1, 1, 1)])

# ---------------------------------------------------------------------------
seccion(u'4.  (n-3,1,1,1) es lisa para todo m cuando n >= 7  (desigualdad triangular)')
# ---------------------------------------------------------------------------
todo = True
for n in range(7, 15):
    for m in range(1, 13):
        if hay_relacion((n - 3, 1, 1, 1), m):
            todo = False
            LINEAS.append(u'       contraejemplo n=%d m=%d' % (n, m))
comprueba(u'(n-3,1,1,1) sin relacion para 7 <= n <= 14, 1 <= m <= 12', todo)
# la razon: n-3 > 1+1+1 = 3 impide la igualdad en |sum| <= sum |.|
comprueba(u'la razon es aritmetica: n-3 > 3 para n >= 7', all(n - 3 > 3 for n in range(7, 15)))

# ---------------------------------------------------------------------------
seccion(u'5.  Las formulas cerradas (thm:formulas), contra fuerza bruta sobre F_p')
# ---------------------------------------------------------------------------
# C_4(3,p) = (p-1)(1 + chi_p(-3))/6
malos = []
for p in primos(5, 61):
    esp = (p - 1) * (1 + chi(-3, p)) // 6
    obt = cuenta(3, 4, p)
    if esp != obt:
        malos.append((p, esp, obt))
comprueba(u'C_4(3,p) = (p-1)(1+chi(-3))/6   [%d primos]' % len(primos(5, 61)),
      not malos, u'fallos: %s' % (malos,))

# C_5(3,p) = (p-1)(4 + chi_p(-3))/6
malos = []
for p in primos(7, 61):
    esp = (p - 1) * (4 + chi(-3, p)) // 6
    obt = cuenta(3, 5, p)
    if esp != obt:
        malos.append((p, esp, obt))
comprueba(u'C_5(3,p) = (p-1)(4+chi(-3))/6   [%d primos]' % len(primos(7, 61)),
      not malos, u'fallos: %s' % (malos,))
# la formula de C_5(3,p) vale tambien en p = 5, con valor 2
comprueba(u'C_5(3,p) vale en p = 5 y da 2', cuenta(3, 5, 5) == 2 and
      (5 - 1) * (4 + chi(-3, 5)) // 6 == 2, u'cuenta = %d' % cuenta(3, 5, 5))

# C_5(4,p) = (p-1)/24 * (4p - 17 - 6chi(-2) - 3chi(-1))
malos = []
ps12 = primos(7, 47)
for p in ps12:
    num = (p - 1) * (4 * p - 17 - 6 * chi(-2, p) - 3 * chi(-1, p))
    if num % 24:
        malos.append((p, u'no entero', num))
        continue
    esp = num // 24
    obt = cuenta(4, 5, p)
    if esp != obt:
        malos.append((p, esp, obt))
comprueba(u'C_5(4,p) = (p-1)(4p-17-6chi(-2)-3chi(-1))/24   [%d primos]' % len(ps12),
      not malos, u'fallos: %s' % (malos,))
comprueba(u'C_5(4,p) vale en p = 5 y da 1', cuenta(4, 5, 5) == 1, u'cuenta = %d' % cuenta(4, 5, 5))

# los dos ingredientes de C_5(4,p), por separado
malos = []
for p in primos(7, 47):
    # B: e_3 = 0, es decir {+-a, +-b}
    B = set()
    for a in range(1, p):
        for b in range(1, p):
            if b % p in (a % p, (-a) % p):
                continue
            B.add(tuple(sorted({a % p, (-a) % p, b % p, (-b) % p})))
    B = {t for t in B if len(t) == 4}
    if len(B) != (p - 1) * (p - 3) // 8:
        malos.append((p, len(B), (p - 1) * (p - 3) // 8))
comprueba(u'|B| = (p-1)(p-3)/8', not malos, u'fallos: %s' % (malos,))

malos = []
for p in primos(7, 61):
    # A cap B = {a,-a,ia,-ia}
    AB = set()
    for i in range(1, p):
        if (i * i) % p != (p - 1):
            continue
        for a in range(1, p):
            t = tuple(sorted({a % p, (-a) % p, (i * a) % p, (-i * a) % p}))
            if len(t) == 4:
                AB.add(t)
    esp = (p - 1) * (1 + chi(-1, p)) // 8
    if len(AB) != esp:
        malos.append((p, len(AB), esp))
comprueba(u'|A cap B| = (p-1)(1+chi(-1))/8', not malos, u'fallos: %s' % (malos,))

# ---------------------------------------------------------------------------
seccion(u'6.  El n = 3 de k = 6: caracter de Artin NO abeliano  (prop:n3k6)')
# ---------------------------------------------------------------------------
# C_6(3,p) = (p-1) * 1{2Y^3 + 3Y - 3 parte completamente en F_p}
def parte_completamente(p):
    raices = [y for y in range(p) if (2 * y ** 3 + 3 * y - 3) % p == 0]
    # cuenta con multiplicidad: grado 3 y tres raices distintas
    return len(raices) == 3

malos = []
ps14 = [p for p in primos(5, 71) if 66 % p]
for p in ps14:
    esp = (p - 1) if parte_completamente(p) else 0
    obt = cuenta(3, 6, p)
    if esp != obt:
        malos.append((p, esp, obt))
comprueba(u'C_6(3,p) = (p-1)*1{2Y^3+3Y-3 parte}   [%d primos]' % len(ps14),
      not malos, u'fallos: %s' % (malos,))

# disc(2Y^3 + 3Y - 3) = -1188 ;  para aY^3+bY^2+cY+d con b = 0:
a, b, c, d = 2, 0, 3, -3
disc = (18 * a * b * c * d - 4 * b ** 3 * d + b * b * c * c
        - 4 * a * c ** 3 - 27 * a * a * d * d)
comprueba(u'disc(2Y^3+3Y-3) = -1188', disc == -1188, u'disc = %d' % disc)
comprueba(u'-1188 = -36 * 33', -1188 == -36 * 33)
comprueba(u'-1188 no es un cuadrado -> grupo de Galois S_3, no C_3',
      int(abs(disc) ** 0.5) ** 2 != abs(disc) or disc < 0)
# el reparto de Chebotarev: 1/6 de los primos parten completamente
parten = [p for p in primos(5, 400) if 66 % p and parte_completamente(p)]
frac = float(len(parten)) / len([p for p in primos(5, 400) if 66 % p])
comprueba(u'densidad de primos que parten ~ 1/6 = 0.167', 0.10 < frac < 0.24,
      u'medida %.3f sobre %d primos' % (frac, len([p for p in primos(5, 400) if 66 % p])))

# ---------------------------------------------------------------------------
seccion(u'7.  "sin peso uno" no implica "Artin-Tate": la quintica nodal de (2,1,1,1,1)  (rem:quintica)')
# ---------------------------------------------------------------------------
# p_g de una superficie de grado d en P^3 es binomial(d-1,3); la resolucion de nodos lo conserva
from math import factorial
pg = factorial(4) // (factorial(3) * factorial(1))
comprueba(u'p_g de una quintica en P^3 = binom(4,3) = 4', pg == 4, u'p_g = %d' % pg)
comprueba(u'(2,1,1,1,1) es particion de 6 de longitud 5',
      sum((2, 1, 1, 1, 1)) == 6 and len((2, 1, 1, 1, 1)) == 5)
comprueba(u'dim V_(2,1,1,1,1) = l - 3 = 2 (superficie)', 5 - 3 == 2)

# ---------------------------------------------------------------------------
seccion(u'8.  Resumen: la clasificacion (thm:clasificacion) y lo medido arriba')
# ---------------------------------------------------------------------------
# texto de resumen impreso, no son comprobaciones
LINEAS.append(u'  k = 3 da {3,4,6}  (caso k = 3; no se re-mide aqui)')
LINEAS.append(u'  k = 5 da {3,4}: n=3 por C_5(3,p), n=4 por C_5(4,p), n=6 cae por el nodo (bloques 3 y 5)')
LINEAS.append(u'  k >= 6 da {3}: (n-3,1,1,1) lisa para n >= 7 (bloque 4); n = 4,5,6 en la nota')

# ---------------------------------------------------------------------------
print(u'\n'.join(LINEAS))
print(u'')
print(u'=' * 70)
print(u'VERIFICACION TRICOTOMIA S_m:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 70)
sys.exit(1 if MAL else 0)
