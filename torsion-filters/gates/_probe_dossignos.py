# -*- coding: utf-8 -*-
# .PUEDEN LOS DOS SIGNOS REALIZAR UNA MISMA RANURA LIBRE?
#
# OBJECION.  "No entry +-2 occurs, so the two signs never both realise a slot" NO se deduce
# de la definicion.  Es cierto, y el error era nuestro: la entrada libre es la SUMA de los eps validos, asi
# que si valen los dos la entrada es (+1)+(-1) = 0, no +-2.  Una entrada 0 esconde exactamente esa
# cancelacion local, y el test que escribi no la ve.
#
# Aritmetica: si los dos signos valen para el mismo (i,j), entonces
#     V_i - t k_1 = X_j = -V_i - t k_2   con k_1, k_2 impares
#  => 2 V_i = t (k_1 - k_2),  y k_1 - k_2 es PAR,  luego V_i es multiplo de t.
# O sea: solo puede pasar en filas cuya clase plegada es 0 -- que son justo las que ningun
# transversal usa, pero que SI son filas de M.  Asi que hay que contarlo, no argumentarlo.
#
# CONTROLES
#   S0  cuantos pares (i,j) admiten los DOS signos, sobre todas las matrices vivas.
#   S1  y de esos, cuantos estan en filas con V_i multiplo de t (deberian ser todos, por la
#       aritmetica de arriba: si sale otro, mi razonamiento esta mal).
#   S2  cuantas entradas valen 0 por esa cancelacion, frente a las que valen 0 por no haber ningun
#       eps.  Es la distincion que la frase del paper daba por hecha.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_dossignos.py

import itertools
import json
import sys
from collections import Counter

from divided_differences import nu_de, plegar
from unimodularidad_barrido import CASOS, columna_nula, matriz

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

S0 = S1 = S2ambos = S2vacio = 0
total_libres = 0
testigos = []

for (t, r, cota) in CASOS:
    mp = (t - 1) // 2
    Rp = mp + r
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu = nu_de(list(Lam), t, r)
        if not nu:
            continue
        V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
        Mx = max(max(abs(v) for v in k) for k in nu)
        permitidos = set()
        for i in range(Rp):
            for e in (1, -1):
                k = 1
                while True:
                    x = e * V[i] - t * k
                    if x < -Mx:
                        break
                    if x <= Mx:
                        permitidos.add(x)
                    k += 2
        for X in itertools.product(sorted(permitidos), repeat=r):
            if any(X[j] <= X[j + 1] for j in range(r - 2)):
                continue
            if not (X[r - 2] > abs(X[r - 1])):
                continue
            if columna_nula(matriz(V, list(X), t, mp, r)):
                continue
            for i in range(Rp):
                for j in range(r):
                    validos = []
                    for e in (1, -1):
                        num = e * V[i] - X[j]
                        if num > 0 and num % t == 0 and (num // t) % 2 == 1:
                            validos.append(e)
                    total_libres += 1
                    if len(validos) == 2:
                        S0 += 1
                        if V[i] % t == 0:
                            S1 += 1
                        S2ambos += 1
                        if len(testigos) < 3:
                            testigos.append({"t": t, "r": r, "Lambda": list(Lam), "X": list(X),
                                             "i": i, "V_i": V[i], "j": j,
                                             "clase": plegar(V[i], t)[0]})
                    elif not validos:
                        S2vacio += 1

print("=" * 96)
print(".PUEDEN LOS DOS SIGNOS REALIZAR LA MISMA RANURA LIBRE?")
print("=" * 96)
print("")
print("  entradas libres examinadas                       : %d" % total_libres)
print("  S0  pares (i,j) con LOS DOS signos validos        : %d" % S0)
print("  S1  de esos, con V_i multiplo de t                : %d" % S1)
print("  S2  entradas 0 por cancelacion de los dos signos  : %d" % S2ambos)
print("      entradas 0 por no haber ningun signo valido   : %d" % S2vacio)
if testigos:
    print("")
    print("  testigos:")
    for w in testigos:
        print("    " + json.dumps(w))
print("")
if S0 == 0:
    print("  LECTURA: no ocurre nunca en la poblacion barrida.  La frase del paper es CIERTA, pero")
    print("  no por la razon que daba: hay que decir que se conto, no deducirlo de la ausencia de 2.")
else:
    print("  LECTURA: SI ocurre.  La frase del paper es falsa y hay que retirarla.")

json.dump({"total_libres": total_libres, "S0": S0, "S1": S1,
           "ceros_por_cancelacion": S2ambos, "ceros_por_vacio": S2vacio,
           "testigos": testigos},
          open("_probe_dossignos_DUMP.json", "w"), indent=1)
print("")
print("DONE")
