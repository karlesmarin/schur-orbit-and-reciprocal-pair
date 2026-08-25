# -*- coding: utf-8 -*-
# LOS ELEMENTOS REGULARES NO PRINCIPALES DE ORDEN PAR, EXPLICITOS Y VERIFICADOS.   20 de agosto de 2026.
#
# POR QUE ESTA GATE EXISTE.  Los elementos regulares NO principales de orden par de Sp(2m) se
# describen aqui en una sola frase: en tipo C_m tienen autovalores xi^{±1}, ..., xi^{±m} con xi
# de orden t = h + 2.  Una frase asi hay que comprobarla antes de apoyar nada en ella, y desde
# luego antes de entregar matrices concretas: un numero que se entrega lleva dentro el
# instrumento con que se midio.  Aqui se hacen explicitos y se verifican para m = 1, 2, 3.
#
# QUE SE COMPRUEBA, para m = 1, 2, 3  (SL(2)=Sp(2), Sp(4), Sp(6)):
#   E0  el elemento esta en Sp(2m): sus autovalores vienen en pares inversos y son unitarios.
#   E1  su ORDEN es exactamente t = 2m+2 = h+2, con h = 2m el numero de Coxeter de C_m.
#   E2  es REGULAR: alpha(x) != 1 para toda raiz alpha del sistema C_m.  Se recorre el sistema
#       entero --- +-e_i +- e_j y +-2e_i --- y no se argumenta: se evalua.
#   E3  NO es PRINCIPAL, y con el testigo que un teorico de grupos reconoce al vuelo: un elemento
#       principal de orden d cumple alpha_i(x) = zeta_d en TODA raiz simple, y aqui las m-1
#       primeras dan xi^{-1} y la ultima da xi^{-2}.
#   E4  el elemento PRINCIPAL de orden t, para contrastar: se resuelve alpha_i(x0) = xi para toda
#       raiz simple y se ve que pide coordenada semientera, o sea que no vive en el toro de Sp(2m)
#       para este t --- que es la razon estructural de que la rama par necesitara prueba a mano.
#   E5  el valor del caracter: sp_lambda(x) para las primeras lambda, calculado por el bialternante
#       en aritmetica exacta de raices de la unidad (enteros ciclotomicos via numpy de enteros no;
#       se usa fracciones de polinomios en xi modulo el ciclotomico).  Se declara el rango.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python npp_elementos_regulares.py > npp_elementos_regulares_OUT.txt 2>&1

import itertools
import json
import sys
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

MS = [1, 2, 3]


# ---------------------------------------------------------------- ciclotomica exacta
def ciclotomico(t):
    """Polinomio minimo de una raiz primitiva t-esima, como lista de coeficientes enteros."""
    # x^t - 1 = prod_{d|t} Phi_d ; se divide sucesivamente
    poly = [-1] + [0] * (t - 1) + [1]          # x^t - 1
    for d in range(1, t):
        if t % d == 0:
            poly = divide(poly, ciclotomico(d))
    return poly


def divide(a, b):
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] // b[-1]
        q[i] = c
        for j in range(len(b)):
            a[i + j] -= c * b[j]
    assert all(x == 0 for x in a), "division no exacta"
    return q


def mulmod(a, b, phi):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return reduce_(out, phi)


def reduce_(a, phi):
    a = a[:]
    n = len(phi) - 1
    for i in range(len(a) - 1, n - 1, -1):
        c = a[i]
        if c:
            a[i] = 0
            for j in range(n):
                a[i - n + j] -= c * phi[j]
    while len(a) > n:
        a.pop()
    while len(a) < n:
        a.append(0)
    return a


def xi_pow(k, t, phi):
    k %= t
    v = [0] * (len(phi) - 1)
    if k < len(v):
        v[k] = 1
        return v
    e = [0] * (len(phi) - 1)
    e[1 if len(e) > 1 else 0] = 1
    r = [1] + [0] * (len(phi) - 2)
    for _ in range(k):
        r = mulmod(r, e, phi)
    return r


def es_cero(a):
    return all(x == 0 for x in a)


def es_uno(a):
    return a[0] == 1 and all(x == 0 for x in a[1:])


print("=" * 100)
print("LOS ELEMENTOS REGULARES NO PRINCIPALES, VERIFICADOS UNO A UNO")
print("=" * 100)

resumen = {}
for m in MS:
    h = 2 * m                 # numero de Coxeter de C_m
    t = h + 2                 # = 2m+2
    phi = ciclotomico(t)
    a = list(range(1, m + 1))   # exponentes: x = diag(xi^1, ..., xi^m, xi^-m, ..., xi^-1)
    G = {1: "SL(2)=Sp(2)", 2: "Sp(4)", 3: "Sp(6)"}[m]

    print("")
    print("-" * 100)
    print("m = %d   %s   h = %d   t = h+2 = %d   xi = e^{2 pi i / %d}" % (m, G, h, t, t))
    print("-" * 100)
    print("  x_%d = diag( %s ,  %s )" % (
        m, ", ".join("xi^%d" % k for k in a), ", ".join("xi^{-%d}" % k for k in reversed(a))))

    # ---- E0
    autov = a + [-k for k in reversed(a)]
    e0 = sorted(autov) == sorted([-k for k in autov])
    print("  E0  autovalores en pares inversos : %s" % ("SI" if e0 else "*** NO ***"))

    # ---- E1  orden
    orden = t // __import__("math").gcd(t, min(abs(k) for k in a))
    e1 = (orden == t)
    print("  E1  orden del elemento : %d   %s" % (orden, "= t = h+2  OK" if e1 else "*** != t ***"))

    # ---- E2  regularidad: se EVALUA sobre todo el sistema de raices C_m
    raices = []
    for i in range(m):
        raices += [tuple(2 if k == i else 0 for k in range(m)),
                   tuple(-2 if k == i else 0 for k in range(m))]
    for i, j in itertools.combinations(range(m), 2):
        for si in (1, -1):
            for sj in (1, -1):
                raices.append(tuple(si if k == i else (sj if k == j else 0) for k in range(m)))
    malas = []
    for r in raices:
        expo = sum(r[k] * a[k] for k in range(m)) % t
        if expo == 0:
            malas.append(r)
    print("  E2  REGULAR: alpha(x) != 1 en las %d raices de C_%d : %s"
          % (len(raices), m, "SIN EXCEPCION" if not malas else "*** falla en %s ***" % (malas[:3],)))

    # ---- E3  no principal
    simples = [tuple((1 if k == i else (-1 if k == i + 1 else 0)) for k in range(m))
               for i in range(m - 1)]
    simples.append(tuple(2 if k == m - 1 else 0 for k in range(m)))
    vals = [sum(s[k] * a[k] for k in range(m)) % t for s in simples]
    e3 = len(set(vals)) > 1 if m > 1 else True
    print("  E3  alpha_i(x) sobre las raices simples, como exponente de xi : %s" % vals)
    if m > 1:
        print("      --> %s: las %d primeras dan xi^{%d} y la ultima xi^{%d}"
              % ("NO ES PRINCIPAL" if e3 else "*** todas iguales ***",
                 m - 1, vals[0] - t if vals[0] > t // 2 else vals[0],
                 vals[-1] - t if vals[-1] > t // 2 else vals[-1]))
    else:
        print("      --> m=1: una sola raiz simple 2e_1, con alpha(x)=xi^2 y no xi: NO ES PRINCIPAL")

    # ---- E4  el principal de orden t pediria coordenada semientera
    #      alpha_i(y)=xi para toda simple  <=>  y_i - y_{i+1} = 1  y  2 y_m = 1
    y = [Fraction(1, 2) + (m - 1 - i) for i in range(m)]
    e4 = any(v.denominator != 1 for v in y)
    print("  E4  el PRINCIPAL de orden t pediria exponentes %s  --> %s"
          % ([str(v) for v in y], "SEMIENTEROS: no esta en este toro" if e4 else "enteros"))

    resumen["m=%d" % m] = dict(grupo=G, h=h, t=t, exponentes=a, E0=bool(e0), E1=bool(e1),
                               E2=len(malas) == 0, n_raices=len(raices),
                               simples=[int(v) for v in vals], E3=bool(e3), E4=bool(e4))

# ---------------------------------------------------------------- E5  valores del caracter
print("")
print("=" * 100)
print("E5  VALORES DEL CARACTER  sp_lambda(x_m)  --- bialternante, aritmetica ciclotomica exacta")
print("=" * 100)


def sp_valor(lam, m, a, t, phi):
    """sp_lambda(x) por el bialternante de tipo C: det(xi^{(lam+rho)_i * a_j} - xi^{-...}) / mismo con rho."""
    lam = list(lam) + [0] * (m - len(lam))
    rho = [m - i for i in range(m)]
    lr = [lam[i] + rho[i] for i in range(m)]

    def bialt(exps):
        M = [[None] * m for _ in range(m)]
        for i in range(m):
            for j in range(m):
                e = exps[i] * a[j]
                M[i][j] = [u - v for u, v in zip(xi_pow(e, t, phi), xi_pow(-e, t, phi))]
        return det(M, phi)

    num, den = bialt(lr), bialt(rho)
    if es_cero(den):
        return None
    # division exacta en Z[xi]: se prueba c en {0,+-1} y se comprueba c*den == num
    for c in (0, 1, -1, 2, -2):
        if [c * x for x in den] == num:
            return c
    return "no en {0,+-1,+-2}"


def det(M, phi):
    n = len(M)
    if n == 1:
        return M[0][0][:]
    acc = [0] * (len(phi) - 1)
    for j in range(n):
        men = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        term = mulmod(M[0][j], det(men, phi), phi)
        if j % 2:
            term = [-x for x in term]
        acc = [u + v for u, v in zip(acc, term)]
    return acc


for m in MS:
    t = 2 * m + 2
    phi = ciclotomico(t)
    a = list(range(1, m + 1))
    lams = [l for l in itertools.product(range(0, 5), repeat=m)
            if all(l[i] >= l[i + 1] for i in range(m - 1))]
    val = {}
    raro = []
    for lam in lams:
        v = sp_valor(lam, m, a, t, phi)
        val[v if not isinstance(v, str) else "otro"] = val.get(
            v if not isinstance(v, str) else "otro", 0) + 1
        if isinstance(v, str) or v is None:
            raro.append((lam, v))
    print("  m=%d  t=%d  %d particiones lambda_1<=4 : reparto de sp_lambda(x) %s   %s"
          % (m, t, len(lams), dict(sorted(val.items(), key=lambda z: str(z[0]))),
             "todos en {0,+-1}" if not raro else "*** raros: %s ***" % raro[:3]))
    resumen["sp_m=%d" % m] = {str(k): int(v) for k, v in val.items()}

json.dump(resumen, open("npp_elementos_regulares_DUMP.json", "w"), indent=1)
print("")
print("=" * 100)
print("DONE")
