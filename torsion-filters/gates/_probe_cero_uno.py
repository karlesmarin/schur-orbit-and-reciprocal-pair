# -*- coding: utf-8 -*-
# .ES M EQUIVALENTE A UNA MATRIZ 0/1?
#
# Las 16 clases canonicas de t=3,r=2 salieron todas con entradas en {0,-1}: ninguna +1.  Si eso vale
# en general, entonces salvo permutaciones con signo de filas y columnas M es una matriz 0/1 -- y la
# unimodularidad total de matrices 0/1 es el terreno clasico (matrices de intervalos, de red,
# Ghouila-Houri), no un problema nuevo.  Eso cambiaria el enunciado de la pregunta abierta.
#
# Pero puede ser un artefacto del canonizador: la forma canonica es el MINIMO lexicografico, que
# prefiere -1 sobre 0 sobre +1.  Que el minimo tenga muchos -1 no implica que no haya ningun +1.
# Asi que se comprueba de verdad, y sobre TODAS las configuraciones del barrido, no solo t=3,r=2.
#
# CONTROLES
#   Z1  toda clase canonica tiene entradas en {0,-1}.
#   Z2  y por tanto -M (o M con las filas negadas) es 0/1: se exhibe el conteo.
#   Z3  senuelo: matrices 0/+-1 ALEATORIAS del mismo tamano y densidad, .cuantas son equivalentes a
#       una 0/1?  Si fuera la mayoria, Z1 no diria nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_cero_uno.py

import itertools
import json
import random
from collections import Counter

from divided_differences import plegar, nu_de, eps_t
from unimodularidad_barrido import (CASOS, canon_signada, matriz, columna_nula)

Z1 = Z1n = 0
Z3 = Z3n = 0
por_caso = []
malos = []
rnd = random.Random(20260816)

for (t, r, cota) in CASOS:
    mp = (t - 1) // 2
    Rp = mp + r
    if Rp > 5:
        por_caso.append({"t": t, "r": r, "cota": cota, "saltado": "R' > 5, canonizar es caro"})
        continue
    vistas = set()
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
            M = matriz(V, X, t, mp, r)
            if columna_nula(M):
                continue
            vistas.add(canon_signada(M))
    ok = 0
    for cl in vistas:
        Z1n += 1
        vals = {v for fila in cl for v in fila}
        if vals <= {0, -1}:
            Z1 += 1
            ok += 1
        elif len(malos) < 3:
            malos.append({"t": t, "r": r, "clase": [list(f) for f in cl]})
    print("  t=%2d r=%d cota %2d :  %4d clases,  %4d con entradas en {0,-1}" % (t, r, cota, len(vistas), ok))
    por_caso.append({"t": t, "r": r, "cota": cota, "clases": len(vistas), "cero_uno": ok})

# Z3  senuelo: matrices 0/+-1 aleatorias con la misma densidad media
dens = 0.45
for n in (3, 4, 5):
    for _ in range(60):
        M = [[rnd.choice([0, 1, -1]) if rnd.random() < dens else 0 for _ in range(n)]
             for _ in range(n)]
        if any(all(M[i][j] == 0 for i in range(n)) for j in range(n)):
            continue
        Z3n += 1
        cl = canon_signada(M)
        if {v for fila in cl for v in fila} <= {0, -1}:
            Z3 += 1

print("")
print("  Z1  toda clase tiene entradas en {0,-1}  : %d de %d" % (Z1, Z1n))
print("  Z3  SENUELO: matrices 0/+-1 aleatorias   : %d de %d" % (Z3, Z3n))
if malos:
    print("")
    print("  !! clases con un +1 irreducible:")
    for m in malos:
        print("    t=%d r=%d" % (m["t"], m["r"]))
        for fila in m["clase"]:
            print("       %s" % fila)
print("")
print("  LECTURA: si Z1 sale entero y Z3 se hunde, M es -- salvo permutaciones con signo -- una")
print("  matriz 0/1, y la conjetura de unimodularidad cae en la teoria clasica en vez de ser nueva.")

json.dump({"Z1": [Z1, Z1n], "Z3_senuelo": [Z3, Z3n], "por_caso": por_caso, "malos": malos},
          open("_probe_cero_uno_DUMP.json", "w"), indent=1)
print("")
print("DONE")
