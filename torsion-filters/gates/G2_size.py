# -*- coding: utf-8 -*-
# EL TAMANO DE G2, que es la PRECONDICION de transplantar la seccion 5 un piso mas abajo.
#
# La seccion 5 de note_t2 descansa entera en |G| <= 2: Purbhoo-van Willigenburg Thm 2.5 es una
# rigidez ENTRE DOS productos de Schur.  Con tres o mas terminos no dice nada.  Antes de reescribir
# nada con G2 hay que saber si |G2| <= 2.
#
# H13  |G2| <= 2.   *** REFUTADA ***  -- ver el _OUT: llega a 10.
#
# LO QUE SE MIDE
#   Q1  distribucion de |G2| sobre TODAS las formas.
#   Q2  distribucion de |G2| sobre la POBLACION OBJETIVO ([Phi]_top = 0 y (i) falsa), que es la
#       unica que importa: son las formas que el estrato de arriba no resuelve.
#   Q3  en cuantas clases difiere cada elemento de G2 del maximizador mas cercano.
#   Q4  no vacuidad: n de la poblacion objetivo se imprime SIEMPRE.
#
# CONSECUENCIA, escrita para que no se me olvide: el corolario del maximizador unico SI se
# transplanta (|G2| = 1 => [Phi]_{D2} es un solo producto no nulo => Phi_t != 0), y cubre 104 de las
# 132 formas objetivo.  |G2| = 2 usa la maquinaria de la seccion 5 tal cual: 8 mas.  |G2| >= 3 son
# 20 formas y quedan ABIERTAS.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python G2_size.py

import itertools
from collections import Counter

from second_stratum import setup, all_transversals, inv_of

CONFIGS = [(4, 2, 15), (4, 3, 15), (6, 2, 15), (6, 3, 17), (8, 2, 17)]


def main():
    dist = Counter()
    dist_t = Counter()
    ncl = Counter()
    tot = 0
    targ = 0

    for (t, r, M) in CONFIGS:
        N = t + 2 * r
        for comb in itertools.combinations(range(M + 1), N):
            beta = tuple(sorted(comb, reverse=True))
            st = setup(beta, t)
            if st is None:
                continue
            cl, E, Cd = st
            if not E:
                continue
            tr = all_transversals(beta, cl, r, t)
            D = max(x[3] for x in tr)
            G = [x for x in tr if x[3] == D]
            rest = [x for x in tr if x[3] < D]
            if not rest:
                continue
            tot += 1
            D2 = max(x[3] for x in rest)
            G2 = [x for x in rest if x[3] == D2]
            dist[len(G2)] += 1
            if len(G) == 2:
                a, b = G
                if inv_of(a[1], r) == inv_of(b[1], r) and a[2] == -b[2]:
                    S = sorted({v for k in E for v in Cd[k]})
                    C = S[0] + S[-1]
                    if set(C - v for v in S) != set(S):
                        targ += 1
                        dist_t[len(G2)] += 1
            for h in G2:
                d = min((sorted(k for k in sorted(cl) if h[0][k] != m[0][k]) for m in G), key=len)
                ncl[len(d)] += 1

    print("formas: %d    Q4 poblacion objetivo: %d" % (tot, targ))
    print("")
    print("Q1  |G2| en TODAS las formas       : %s" % dict(sorted(dist.items())))
    print("Q2  |G2| en la POBLACION OBJETIVO  : %s" % dict(sorted(dist_t.items())))
    print("Q3  clases en que difiere del maximizador mas cercano : %s" % dict(sorted(ncl.items())))
    print("")
    print("H13 |G2| <= 2 : *** REFUTADA ***  (max observado = %d)" % max(dist))
    n1 = dist_t.get(1, 0)
    n2 = dist_t.get(2, 0)
    print("")
    print("ALCANCE de la seccion 5 transplantada, sobre las %d formas objetivo:" % targ)
    print("   |G2| = 1  -> corolario del maximizador unico, Phi_t != 0 : %d" % n1)
    print("   |G2| = 2  -> maquinaria de la seccion 5 (PvW)            : %d" % n2)
    print("   |G2| >= 3 -> ABIERTO, PvW no aplica                      : %d" % (targ - n1 - n2))


if __name__ == "__main__":
    main()
