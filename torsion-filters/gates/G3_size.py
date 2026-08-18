# -*- coding: utf-8 -*-
# EL TAMANO DEL TERCER ESTRATO.  Se puede controlar como se controlo |G| <= 2?
#
# POR QUE G3 NO ES ANALOGO A G2.  El indice correcto de un estrato NO son las transversales sino los
# PARES (g, A) de la expansion de Laplace por columnas pares:
#
#     [Phi_t]_D = sum_{(g,A) : 2 sum A - sum T_g = D}  w(g) eps(A) a_A(z) a_{T_g\A}(1/z)
#
# En D1 solo contribuyen los maximizadores con su mitad superior, asi que |pares| = |G| <= 2.
# En D2 tampoco entran los maximizadores: el teorema gamma < g_int/2 lo prohibe, luego |pares| =
# |G2|.  Pero en D3 = D1 - 4 (cuando gamma = 2) hay DOS fuentes:
#     (a) transversales de grado D1-4 con su mitad superior;
#     (b) las de G2 con su SEGUNDO subconjunto, si su hueco interno vale 2.
# Los maximizadores siguen fuera: gamma < g_int/2 con gamma = 2 fuerza g_int >= 6, o sea su primera
# caida interna esta en D1-6.  Esa asimetria es la novedad del tercer piso y hay que medirla, no
# suponerla.
#
# LO QUE SE MIDE
#   X0  CONTROL: en D1 el numero de pares tiene que ser exactamente |G|, y en D2 exactamente |G2|.
#       Si no, el indice esta mal montado y lo demas sobra.
#   X1  distribucion de |G3| (pares en el tercer grado distinto).
#   X2  DESGLOSE por fuente: cuantos pares de D3 vienen de un maximizador, de G2, o de una
#       transversal de grado D1-4.  Es la pregunta estructural.
#   X3  el perfil de las formas de PROFUNDIDAD 4 -- aquellas donde el tercer estrato es el que
#       prueba Phi_t != 0.  Ahi es donde G3 tiene que ser manejable.
#   X4  no vacuidad: n de cada casilla, y cuantas formas tienen menos de tres grados.
#   X5  SENUELO: |G3| <= 2, el analogo directo del teorema del estrato de arriba.  Se espera que
#       FALLE -- |G2| ya llegaba a 10 -- y si no fallara habria que desconfiar del indice.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python G3_size.py

import itertools
from collections import Counter, defaultdict

from second_stratum import setup, all_transversals
from depth import deltas, dims

CONFIGS = [(4, 2, 15), (4, 3, 15), (6, 2, 16), (6, 3, 17), (8, 2, 17), (6, 4, 18)]


def strata_pairs(beta, t, r):
    """{grado: [(indice de transversal, A, B, w, eps)]} y la lista de transversales."""
    st = setup(beta, t)
    if st is None:
        return None, None
    cl, E, Cd = st
    if not E:
        return None, None
    tr = all_transversals(beta, cl, r, t)
    idx = range(2 * r)
    by = defaultdict(list)
    for gi, (_, T, w, dg) in enumerate(tr):
        sT = sum(T)
        for R in itertools.combinations(idx, r):
            A = tuple(T[a] for a in R)
            B = tuple(T[a] for a in idx if a not in R)
            D = 2 * sum(A) - sT
            eps = -1 if (sum(R) % 2) else 1
            by[D].append((gi, A, B, w, eps))
    return by, tr


def main():
    x0_bad = x0_n = 0
    x1 = Counter()
    x2 = Counter()
    x5_bad = x5_n = 0
    short = 0
    tot = 0
    depth4 = []

    for (t, r, M) in CONFIGS:
        N = t + 2 * r
        for comb in itertools.combinations(range(M + 1), N):
            beta = tuple(sorted(comb, reverse=True))
            by, tr = strata_pairs(beta, t, r)
            if by is None:
                continue
            tot += 1
            degs = sorted(by, reverse=True)
            if len(degs) < 3:
                short += 1
                continue
            D1, D2, D3 = degs[0], degs[1], degs[2]
            dmax_t = max(x[3] for x in tr)
            G = [i for i, x in enumerate(tr) if x[3] == dmax_t]
            rest = [x[3] for x in tr if x[3] < dmax_t]
            G2 = [i for i, x in enumerate(tr) if rest and x[3] == max(rest)]

            # X0
            x0_n += 1
            if len(by[D1]) != len(G) or (rest and len(by[D2]) != len(G2)):
                x0_bad += 1

            n3 = len(by[D3])
            x1[n3] += 1
            x5_n += 1
            if n3 > 2:
                x5_bad += 1

            Gs, G2s = set(G), set(G2)
            for (gi, A, B, w, eps) in by[D3]:
                if gi in Gs:
                    x2['de un maximizador'] += 1
                elif gi in G2s:
                    x2['de G2 (segundo subconjunto)'] += 1
                else:
                    x2['de una transversal mas baja'] += 1

            # profundidad, reusando los pares ya construidos en vez de reexpandir
            d = None
            for k in range(0, 13, 2):
                D = dmax_t - k
                if D in by and sum(w * eps * dims(A, B, r)
                                   for (_, A, B, w, eps) in by[D]) != 0:
                    d = k
                    break
            if d == 4:
                depth4.append((t, r, beta, n3, len(G2)))
        print("   hecho t=%d r=%d M=%d (acumulado %d)" % (t, r, M, tot), flush=True)

    print("")
    print("X0 CONTROL |pares en D1| == |G|  y  |pares en D2| == |G2| : %d fallos de %d"
          % (x0_bad, x0_n))
    print("X4 formas con menos de tres grados distintos (sin D3) : %d de %d" % (short, tot))
    print("")
    print("X1 distribucion de |G3| (pares en el tercer grado):")
    for k in sorted(x1):
        print("      %3d : %d" % (k, x1[k]))
    print("   maximo observado: %d" % (max(x1) if x1 else 0))
    print("")
    print("X5 SENUELO |G3| <= 2 : %d fallos de %d" % (x5_bad, x5_n))
    print("")
    print("X2 DESGLOSE por fuente de los pares de D3:")
    for k in sorted(x2):
        print("      %-32s : %d" % (k, x2[k]))
    print("")
    print("X3 formas de PROFUNDIDAD 4 (el tercer estrato es el que decide) : %d" % len(depth4))
    if depth4:
        c = Counter(x[3] for x in depth4)
        print("      |G3| en ellas : %s" % dict(sorted(c.items())))
        c2 = Counter(x[4] for x in depth4)
        print("      |G2| en ellas : %s" % dict(sorted(c2.items())))
        for x in depth4[:6]:
            print("      t=%d r=%d |G3|=%d |G2|=%d beta=%s" % (x[0], x[1], x[3], x[4], list(x[2])))


if __name__ == "__main__":
    main()
