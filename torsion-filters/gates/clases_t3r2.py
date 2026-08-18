# -*- coding: utf-8 -*-
# LAS CLASES DE MATRICES A t=3, r=2, IMPRESAS.
#
# unimodularidad_barrido.py encontro que a t=3, r=2 los 1191 X vivos dan 220 matrices distintas y
# solo 16 CLASES salvo permutaciones con signo de filas y columnas -- la equivalencia exacta que la
# unimodularidad total no distingue.  Aqui se imprimen.  Si la conjetura es demostrable por
# inspeccion de casos, es sobre estas 16 sobre las que se demuestra.
#
# R' = m' + r = 1 + 2 = 3, asi que son matrices 3x3:
#   columna 0  la ranura CONGELADA, que pide la clase plegada 1;   entrada = signo del plegado
#   columnas 1,2  las ranuras LIBRES, que piden X_1 y X_2;         entrada = suma de los eps validos
#
# Se imprime, por clase: el representante canonico, |det|, si es totalmente unimodular, cuantas
# matrices distintas caen en ella y cuantos pares (Lambda, X) vivos, y un (Lambda, X) testigo.
#
# CONTROLES
#   K1  toda clase es totalmente unimodular (si alguna no lo fuera, la conjetura estaria muerta y
#       este listado seria el contraejemplo).
#   K2  |det| en {0,1} en todas.
#   K3  la suma de los pares por clase reproduce el total de X vivos del barrido (1191 a beta<=6).
#   K4  se repite con la caja grande (beta<=12) para ver que clase NUEVA aparece -- el barrido dio
#       16 clases a <=6 y 17 a <=12, y esa diferencia hay que exhibirla, no solo contarla.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python clases_t3r2.py

import itertools
import json
from collections import Counter, defaultdict

from divided_differences import plegar, nu_de, nu_extendida, eps_t
from unimodularidad_barrido import (det_entero, es_TU, canon_signada, matriz, columna_nula,
                                    c_cerrada)

T, R = 3, 2
MP = (T - 1) // 2
RP = MP + R


def recoger(cota):
    """clase canonica -> (matrices distintas, pares vivos, testigo, det, TU)."""
    clases = defaultdict(lambda: {"distintas": set(), "pares": 0, "testigo": None})
    for Lam in itertools.product(range(cota + 1), repeat=RP):
        if any(Lam[i] < Lam[i + 1] for i in range(RP - 1)):
            continue
        nu = nu_de(list(Lam), T, R)
        if not nu:
            continue
        V = [2 * Lam[i] + 2 * (RP - i) - 1 for i in range(RP)]
        Mx = max(max(abs(v) for v in k) for k in nu)
        permitidos = set()
        for i in range(RP):
            for e in (1, -1):
                k = 1
                while True:
                    x = e * V[i] - T * k
                    if x < -Mx:
                        break
                    if x <= Mx:
                        permitidos.add(x)
                    k += 2
        for X in itertools.product(sorted(permitidos), repeat=R):
            if not (X[R - 2] > abs(X[R - 1])):
                continue
            M = matriz(V, X, T, MP, R)
            if columna_nula(M):
                continue
            cl = canon_signada(M)
            d = clases[cl]
            d["distintas"].add(tuple(tuple(f) for f in M))
            d["pares"] += 1
            if d["testigo"] is None:
                d["testigo"] = (list(Lam), list(X), c_cerrada(nu, X, T, R, Mx))
    return clases


def imprimir(clases, titulo):
    print("")
    print("=" * 96)
    print(titulo)
    print("=" * 96)
    orden = sorted(clases, key=lambda k: (abs(det_entero([list(f) for f in k])), k))
    K1 = K2 = 0
    total_pares = 0
    for n, cl in enumerate(orden, 1):
        M = [list(f) for f in cl]
        d = det_entero(M)
        tu = es_TU(M)
        K1 += 1 if tu else 0
        K2 += 1 if abs(d) <= 1 else 0
        info = clases[cl]
        total_pares += info["pares"]
        Lam, X, cval = info["testigo"]
        print("")
        print("  clase %2d :  det = %+d   TU %s   |  %3d matrices distintas, %4d pares (Lambda,X)"
              % (n, d, "SI" if tu else "NO !!", len(info["distintas"]), info["pares"]))
        print("              testigo  Lambda=%s  X=%s   ->  c = %+d" % (Lam, X, cval))
        for fila in M:
            print("                 [ %s ]" % "  ".join("%+d" % v for v in fila))
    print("")
    print("  K1  toda clase totalmente unimodular : %d de %d" % (K1, len(orden)))
    print("  K2  |det| en {0,1}                   : %d de %d" % (K2, len(orden)))
    print("  K3  pares (Lambda,X) vivos sumados   : %d" % total_pares)
    return orden


c6 = recoger(6)
o6 = imprimir(c6, "t = 3, r = 2, beta_i <= 6   --  LAS 16 CLASES")

c12 = recoger(12)
o12 = imprimir(c12, "t = 3, r = 2, beta_i <= 12  --  la caja grande")

print("")
print("=" * 96)
print("K4  .QUE CLASE APARECE AL AGRANDAR LA CAJA?")
print("=" * 96)
nuevas = [k for k in o12 if k not in c6]
idas = [k for k in o6 if k not in c12]
print("  clases a beta<=6 : %d      clases a beta<=12 : %d" % (len(o6), len(o12)))
print("  nuevas al agrandar : %d     desaparecidas : %d" % (len(nuevas), len(idas)))
for k in nuevas:
    M = [list(f) for f in k]
    Lam, X, cval = c12[k]["testigo"]
    print("")
    print("  NUEVA :  det = %+d   TU %s   testigo Lambda=%s X=%s  ->  c = %+d"
          % (det_entero(M), "SI" if es_TU(M) else "NO !!", Lam, X, cval))
    for fila in M:
        print("             [ %s ]" % "  ".join("%+d" % v for v in fila))

json.dump({"clases_beta6": [[list(f) for f in k] for k in o6],
           "clases_beta12": [[list(f) for f in k] for k in o12],
           "nuevas_al_agrandar": [[list(f) for f in k] for k in nuevas]},
          open("clases_t3r2_DUMP.json", "w"), indent=1)
print("")
print("=" * 96)
print("DONE")
