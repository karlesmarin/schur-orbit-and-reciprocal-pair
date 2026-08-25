# -*- coding: utf-8 -*-
# SONDA: .por que univ_sp devuelve None?   19 de agosto de 2026.
# No decide nada; enseña las piezas.  Plano, sin if/else a nivel superior (REPL).

import sys

T, NLET, LMAX = 3, 5, 12
kmax = 3 * LMAX + 12
K = CyclotomicField(6, 'zt')
zt = K.gen()
w = zt ** 2
LzK = LaurentPolynomialRing(K, 'z')
zK = LzK.gen()
orbita = [LzK(w ** j) for j in range(T)] + [zK, zK ** -1]

print("kmax =", kmax)
P = orbita[0].parent()
print("parent de las letras:", P)
S = PowerSeriesRing(P, 'u', default_prec=kmax + 2)
u = S.gen()
f = S(1)
for a in orbita:
    f *= sum((a * u) ** j for j in range(kmax + 2))
print("precision de f :", f.prec())
hs = [f[k] for k in range(kmax + 1)]
print("len(hs) =", len(hs))
print("hs[0] =", hs[0])
print("hs[1] =", hs[1])
print("hs[13] =", hs[13])
print("hs[48] =", hs[48])
nones = [k for k in range(len(hs)) if hs[k] is None]
print("posiciones None en hs:", nones[:10], "total", len(nones))
sys.stdout.flush()

# el determinante a mano para lambda = (1,0)
parts = [1]
L = 1
print("")
print("lambda=(1,0): indices que pide univ_sp:")
for i in range(L):
    for j in range(L):
        print("   x = hh(%d),  y = hh(%d)" % (parts[i] - (i + 1) + (j + 1), parts[i] - (i + 1) - (j + 1) + 2))
print("")
print("DONE")
