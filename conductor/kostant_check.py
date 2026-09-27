"""Compuerta de la atribucion a Kostant (27-sep-2026).

La nota dice (Prop. AZ y la seccion de historia) que para q = 2m+2 el elemento g_{1/q} de Sp(2m) es el de
Kostant [Kostant76, Thm 3.1], en el que todo caracter irreducible vale 0 o +-1.  El enunciado de Kostant, tal
como lo reproducen Cellini-Moseneder Frajria-Papi (math/0507610, Thm 1.1 con la normalizacion (1.1)), usa
a = exp(2 pi i * 2 rho), con rho llevado a la subalgebra de Cartan por la forma normalizada (theta,theta) = 1/h^v.

Este guion comprueba, para Sp(2m):
  (1) con esa normalizacion, a = exp(2 pi i rho/(2m+2)), rho = (m, ..., 1): sus autovalores en la
      representacion estandar son xi^{+-j}, j = 1..m, xi = exp(2 pi i/(2m+2)); es decir, a = g_{1/(2m+2)};
      y exp(2 pi i rho^v/(2m+2)), rho^v = (m-1/2, ..., 1/2), NO es ese elemento;
  (2) chi_lambda(g_{1/(2m+2)}) esta en {-1, 0, 1} para todo peso dominante de la caja, por la formula de
      Weyl de tipo C (cociente de determinantes), en aritmetica compleja con tolerancia 1e-9.
"""
import cmath
import itertools
from fractions import Fraction


def det(M):
    n = len(M)
    M = [row[:] for row in M]
    d = 1
    for c in range(n):
        piv = max(range(c, n), key=lambda i: abs(M[i][c]))
        if abs(M[piv][c]) < 1e-14:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            for j in range(c, n):
                M[i][j] -= f * M[c][j]
    return d


fallos = 0
for m in range(1, 6):
    # (1) normalizacion: forma (u,v) = k * u.v en R^m, raiz maxima theta = 2 e_1, h^v = m + 1
    k = Fraction(1, 4 * (m + 1))          # (theta,theta) = 4k = 1/h^v
    rho = [m - i for i in range(m)]        # (m, ..., 1)
    x = [2 * k * r for r in rho]           # lambda(x) = (lambda, 2 rho): x = 2k rho
    esperado = [Fraction(m - i, 2 * m + 2) for i in range(m)]
    rho_v = [Fraction(2 * (m - i) - 1, 2) for i in range(m)]
    otro = [r / (2 * m + 2) for r in rho_v]
    ok1 = (x == esperado) and (otro != esperado)
    fallos += not ok1
    print("m=%d  a = exp(2 pi i x), x = %s  == rho/(2m+2): %s;  rho^v/(2m+2) = %s distinto: %s"
          % (m, [str(t) for t in x], x == esperado, [str(t) for t in otro], otro != esperado))

    # (2) caracteres en g_{1/(2m+2)}
    q = 2 * m + 2
    z = [cmath.exp(2j * cmath.pi * j / q) for j in range(1, m + 1)]
    den = det([[z[c] ** (m - i) - z[c] ** (-(m - i)) for c in range(m)] for i in range(m)])
    tope = 3 * q
    valores, casos = set(), 0
    for lam in itertools.product(range(tope), repeat=m):
        if any(lam[i] < lam[i + 1] for i in range(m - 1)):
            continue
        l = [lam[i] + m - i for i in range(m)]
        num = det([[z[c] ** l[i] - z[c] ** (-l[i]) for c in range(m)] for i in range(m)])
        v = num / den
        r = round(v.real)
        if abs(v - r) > 1e-9 or r not in (-1, 0, 1):
            fallos += 1
            print("   FALLO m=%d lambda=%s valor=%s" % (m, lam, v))
        valores.add(r)
        casos += 1
    print("      pesos dominantes con lambda_1 < %d: %d; valores de chi en g_{1/%d}: %s"
          % (tope, casos, q, sorted(valores)))
print("fallos = %d" % fallos)
