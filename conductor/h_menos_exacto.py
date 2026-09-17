# Copyright (c) 2026 Carles Marin. Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# h^-(Q(zeta_n)) EXACTO para n primo impar, contra el valor en coma flotante de bernoulli_W1.h_menos_numerico
# (el que usa la Figura 4 de la nota del conductor).
#
# Metodo exacto (Washington, Thm 4.17): h^- = Q w prod_{chi impar} (-1/2 B_{1,chi}), con Q = 1, w = 2n para n primo.
# Se agrupan los caracteres impares en orbitas de Galois: chi(g) = zeta_{n-1}^k con k impar; la orbita de k es
# {k u : u unidad mod n-1}.  Para cada orbita el producto de -B_{1,chi}/2 es la norma de un elemento de
# Q(zeta_d), d = (n-1)/gcd(k, n-1), que se calcula como resultante con el polinomio ciclotomico (sympy, enteros
# y racionales exactos).  El producto total es un racional; se exige que sea entero.
#
# Uso: python h_menos_exacto.py [NMAX]
import os
import sys
from fractions import Fraction
from math import gcd

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bernoulli_W1 import h_menos_numerico, es_primo, primitiva  # noqa: E402

x = sp.symbols("x")


def h_menos_exacto(n):
    g = primitiva(n)
    log = {}
    y = 1
    for i in range(n - 1):
        log[y] = i
        y = y * g % n
    N = n - 1
    total = Fraction(1)
    vistos = set()
    for k in range(1, N, 2):  # caracteres impares: chi(g) = zeta_N^k con k impar
        if k in vistos:
            continue
        d = N // gcd(k, N)
        orbita = {(k * u) % N for u in range(1, N) if gcd(u, N) == 1}
        vistos |= orbita
        # chi(a) = zeta_d^{(k/gcd) * log a}; tomando el representante primitivo, B_{1,chi} = (1/n) sum a zeta_d^{j log a}
        j = k // gcd(k, N)
        # polinomio P(x) = sum_a a x^{(j log a) mod d}; -B/2 = -P(zeta_d)/(2n)
        coef = {}
        for a in range(1, n):
            e = (j * log[a]) % d
            coef[e] = coef.get(e, 0) + a
        P = sum(c * x**e for e, c in coef.items())
        Phi = sp.cyclotomic_poly(d, x)
        norma = sp.resultant(Phi, P, x)  # = prod_{zeta primitivo} P(zeta), salvo signo (+1: Phi monico)
        grado = sp.totient(d)
        total *= Fraction(int(norma)) * Fraction(-1, 2 * n) ** int(grado)
    h = Fraction(2 * n) * total
    assert h.denominator == 1, (n, h)
    return int(h)


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 90
    malos = 0
    for n in range(5, NMAX):
        if not es_primo(n):
            continue
        he = h_menos_exacto(n)
        hn = h_menos_numerico(n)
        ok = he == hn
        malos += not ok
        print("q'=%3d  h^- exacto %s  flotante %s  %s" % (n, he, hn, "igual" if ok else "*** DISTINTO ***"))
    print("primos 5 <= q' < %d: discrepancias %d" % (NMAX, malos))


if __name__ == "__main__":
    main()
