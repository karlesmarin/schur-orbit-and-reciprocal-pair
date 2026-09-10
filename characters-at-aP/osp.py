# -*- coding: utf-8 -*-
"""LA PISTA DE osp(1|2m), verificada por dimensiones.

Afirma:  h = osp(1|2m)^{Ad(-s)}  =  (+)_c gl_{n_c}  (+)  sp_{2n_0}  (+)  osp(1|2n_{q/2}),

con s = diag(xi^{l_1},...,xi^{l_m}, xi^{-l_1},...,xi^{-l_m}).

Razon: la parte PAR de osp(1|2m) es sp_{2m} y -s actua alli igual que s, asi que su fijo es el
centralizador habitual.  La parte IMPAR es la natural de dimension 2m, y ser fijo por -s es estar
en el autoespacio de -1 de s: solo la clase q/2 sobrevive.  Y sp_{2n} + (2n impares) = osp(1|2n).

Se comprueba: dim de cada parte, calculada a mano sobre los autoespacios, contra la suma de las
dimensiones de los bloques anunciados.
"""
import sys
from collections import Counter
from math import gcd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def clases(ell, q):
    """ocupaciones n_c de las clases plegadas, y cuales son fijas."""
    oc = Counter()
    for x in ell:
        r = x % q
        oc[min(r, q - r)] += 1
    fij = {0} | ({q // 2} if q % 2 == 0 else set())
    return oc, fij


def dim_fijo_par(ell, q):
    """dim del centralizador de s en sp_2m, contado sobre los autoespacios de la natural.
    autovalores de la natural: xi^{+-l_i}.  El centralizador de un elemento semisimple de Sp
    con autoespacios V_a (a != +-1, emparejados con V_{a^-1}) y V_{+1}, V_{-1} es
       (+)_{pares} gl(dim V_a)  (+)  sp(dim V_1)  (+)  sp(dim V_{-1}).
    """
    mult = Counter()
    for x in ell:
        mult[x % q] += 1
        mult[(-x) % q] += 1
    d = 0
    vistos = set()
    for a, n in mult.items():
        if a in vistos:
            continue
        b = (-a) % q
        if a == b:                      # a = +1 (a=0) o a = -1 (a=q/2)
            d += n * (n + 1) // 2       # sp_n con n par: dim = n(n+1)/2
            vistos.add(a)
        else:
            d += n * n                  # gl_n, con n = mult(a) = mult(b)
            vistos.add(a); vistos.add(b)
    return d


def dim_fijo_impar(ell, q):
    """la parte impar es la natural; fijo por Ad(-s) = autoespacio de -1 de s."""
    if q % 2:
        return 0
    mult = Counter()
    for x in ell:
        mult[x % q] += 1
        mult[(-x) % q] += 1
    return mult.get(q // 2, 0)


def dim_anunciada(ell, q):
    oc, fij = clases(ell, q)
    par = 0
    impar = 0
    for c, n in oc.items():
        if c == 0:
            par += n * (2 * n + 1)                 # sp_{2n}
        elif c in fij:
            par += n * (2 * n + 1)                 # sp_{2n} dentro de osp(1|2n)
            impar += 2 * n                          # la parte impar de osp(1|2n)
        else:
            par += n * n                            # gl_n
    return par, impar


print("=" * 92)
print("h = osp(1|2m)^{Ad(-s)}  contra  (+) gl_{n_c} (+) sp_{2n_0} (+) osp(1|2n_{q/2})")
print("=" * 92)
mal = tot = 0
for m in range(2, 8):
    rho = [m - i for i in range(m)]
    for q in range(2, 13):
        for lam0 in range(0, 6):
            for lam1 in range(0, lam0 + 1):
                lam = [lam0, lam1] + [0] * (m - 2)
                ell = [lam[i] + rho[i] for i in range(m)]
                p1, i1 = dim_fijo_par(ell, q), dim_fijo_impar(ell, q)
                p2, i2 = dim_anunciada(ell, q)
                tot += 1
                if (p1, i1) != (p2, i2):
                    mal += 1
                    if mal <= 5:
                        print("   FALLA m=%d q=%d ell=%s  medido=(%d,%d) anunciado=(%d,%d)"
                              % (m, q, ell, p1, i1, p2, i2))
print("   casos: %d   discrepancias: %d" % (tot, mal))
print()
print("   y el ejemplo de la §13.1:  m=6, q=2")
m, q = 6, 2
rho = [m - i for i in range(m)]
lam = [1, 1, 0, 0, 0, 0]
ell = [lam[i] + rho[i] for i in range(m)]
oc, fij = clases(ell, q)
print("      ell=%s   ocupaciones plegadas %s" % (ell, dict(oc)))
print("      par medido %d, impar medido %d" % (dim_fijo_par(ell, q), dim_fijo_impar(ell, q)))
print("      anunciado: sp_%d (+) osp(1|%d)  ->  (%d, %d)"
      % (2 * oc[0], 2 * oc[1], *dim_anunciada(ell, q)))
