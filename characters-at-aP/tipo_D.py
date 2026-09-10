# -*- coding: utf-8 -*-
"""ESCRIBI «en tipos B y D el adjunto es Lambda^2 ... con A_x = P_x + 1».

El +1 es el peso CERO de la natural de dimension 2m+1, que es de tipo B.  En tipo D la natural
tiene dimension 2m y NO tiene peso cero.  Mi guion de verificacion (simplifica2.py) probo B y C
y NUNCA probo D --- y yo escribi «B y D» sobre una medida que no cubria D.

Aqui se mide D de verdad, con las dos formulas.
"""
import random, sys
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def raices_D(m):
    R = []
    for i in range(m):
        for j in range(i + 1, m):
            for s in (1, -1):
                for t in (1, -1):
                    v = [0] * m; v[i] = s; v[j] = t; R.append(tuple(v))
    return R


def lado_izq(x, R):
    c = Counter()
    for a in R:
        c[sum(a[i] * x[i] for i in range(len(x)))] += 1
    return c


def lam2(pesos):
    c = Counter()
    for i in range(len(pesos)):
        for j in range(i + 1, len(pesos)):
            c[pesos[i] + pesos[j]] += 1
    return c


print("=" * 92)
print("TIPO D:  suma sobre raices  ==  char Lambda^2(natural) - m ?")
print("   natural de D_m: pesos +-x_i, dimension 2m, SIN peso cero  ->  A_x = P_x")
print("   la que escribi:  A_x = P_x + 1   (que es la de B)")
print("=" * 92)
random.seed(5)
for m in range(2, 8):
    R = raices_D(m)
    mal_D = mal_B = 0
    for _ in range(300):
        x = [random.randint(-30, 30) for _ in range(m)]
        izq = lado_izq(x, R)
        # correcta: Lambda^2 de {+-x_i}
        der = lam2([+v for v in x] + [-v for v in x])
        der[0] -= m
        der = Counter({k: v for k, v in der.items() if v})
        if izq != der:
            mal_D += 1
        # la que escribi: Lambda^2 de {+-x_i, 0}
        derB = lam2([+v for v in x] + [-v for v in x] + [0])
        derB[0] -= m
        derB = Counter({k: v for k, v in derB.items() if v})
        if izq != derB:
            mal_B += 1
    print("  D_%d : #raices=%-3d | con A=P_x  discrepancias %-4d | con A=P_x+1  discrepancias %-4d"
          % (m, len(R), mal_D, mal_B))

print()
print("=" * 92)
print("EL CONTROL DE UNA LINEA, evaluando en X=1  (cada termino cuenta multiplicidades)")
for m in (3, 4, 5):
    n_raices = len(raices_D(m))
    P1 = 2 * m
    A1 = 2 * m + 1
    print("  D_%d : raices reales = %-3d |  (P^2-P)/2 - m = %-4d  |  (A^2-A)/2 - m = %-4d  %s"
          % (m, n_raices, (P1 * P1 - P1) // 2 - m, (A1 * A1 - A1) // 2 - m,
             "<- la que escribi, y falla" if (A1 * A1 - A1) // 2 - m != n_raices else ""))
