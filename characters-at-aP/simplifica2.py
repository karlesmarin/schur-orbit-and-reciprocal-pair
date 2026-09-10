# -*- coding: utf-8 -*-
"""Igual que simplifica.py pero con x_i ENTEROS y comparando multiconjuntos de exponentes:
sympy no puede decidir con exponentes simbolicos, y su 'False' era del instrumento."""
import sys, random
from collections import Counter
from itertools import product
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def raices(tipo, m):
    """C: e_i+-e_j y 2e_i.   B: e_i+-e_j y e_i.   D: solo e_i+-e_j (sin raices cortas)."""
    R = []
    for i in range(m):
        for j in range(i + 1, m):
            for s in (1, -1):
                for t in (1, -1):
                    v = [0] * m; v[i] = s; v[j] = t; R.append(tuple(v))
    if tipo in ("B", "C"):
        for i in range(m):
            for s in (1, -1):
                v = [0] * m
                v[i] = 2 * s if tipo == "C" else s
                R.append(tuple(v))
    return R


def lado_izq(x, R):
    """sum_{alpha in Phi} X^{(x,alpha)}  como Counter de exponentes."""
    c = Counter()
    for a in R:
        c[sum(a[i] * x[i] for i in range(len(x)))] += 1
    return c


def sym2(x):
    """char Sym^2 del estandar de Sp(2m): pesos +-x_i, pares CON repeticion."""
    pesos = [+xi for xi in x] + [-xi for xi in x]
    c = Counter()
    for i in range(len(pesos)):
        for j in range(i, len(pesos)):
            c[pesos[i] + pesos[j]] += 1
    return c


def lam2(x):
    """char Lambda^2 del estandar de SO(2m+1): pesos +-x_i Y EL CERO."""
    return _lam2([+xi for xi in x] + [-xi for xi in x] + [0])


def lam2_D(x):
    """char Lambda^2 del estandar de SO(2m): pesos +-x_i, SIN peso cero.
    ⚠ Esta era la que faltaba.  El articulo decia «B y D» con el +1 de B, y en D_4 daba 32
    raices donde hay 24.  Un tipo que el guion no probaba."""
    return _lam2([+xi for xi in x] + [-xi for xi in x])


def _lam2(pesos):
    c = Counter()
    for i in range(len(pesos)):
        for j in range(i + 1, len(pesos)):
            c[pesos[i] + pesos[j]] += 1
    return c


print("=" * 92)
print("sum_{alpha in Phi} X^{(x,alpha)}  ==  char(adjunto)(x) - rango ?")
print("   tipo C : adjunto = Sym^2(estandar)      tipo B : adjunto = Lambda^2(estandar)")
print("=" * 92)
random.seed(11)
for tipo, adj in (("C", sym2), ("B", lam2), ("D", lam2_D)):
    for m in range(2, 8):
        R = raices(tipo, m)
        malos = 0
        for _ in range(400):
            x = [random.randint(-40, 40) for _ in range(m)]
            izq = lado_izq(x, R)
            der = adj(x)
            der[0] -= m                       # quitar los 'rango' ceros
            der = Counter({k: v for k, v in der.items() if v})
            if izq != der:
                malos += 1
        print("  %s_%d :  400 vectores aleatorios,  discrepancias: %d %s"
              % (tipo, m, malos, "" if malos else " <-- identidad"))
print()
print("=" * 92)
print("Y por tanto el 'lema de residuos' del articulo es, palabra por palabra:")
print("     R_x(X) = char(adjunto)(x) - rango,")
print("y su segunda mitad (R_ell = R_rho) es:  el CARACTER ADJUNTO de ell y de rho")
print("coinciden en Z[X]/(X^q-1).  Eso no menciona el tipo.")
