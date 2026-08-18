# -*- coding: utf-8 -*-
# EL CERTIFICADO DE DIMENSION, y con el un barrido que llega donde el polinomio no llegaba.
#
# LA IDEA.  Las relaciones son LINEALES, asi que se pueden mirar con grupos algebraicos.
# Dividiendo por Vdm^2, el diccionario manda P(T) = a_H(z) a_L(1/z) a  (prod z)^d s_atil s_astar,
# que vive en el ANILLO DE REPRESENTACIONES de GL(r).  Alli
#
#     sum_{g in G2} w(g) P(T_g) = 0     <=>     SUM_{w=+1} s_a s_b  =  SUM_{w=-1} s_a s_b
#
# o sea: dos SUMAS DIRECTAS de productos tensoriales de GL(r) son isomorfas.  Y dos representaciones
# son isomorfas si y solo si coinciden sus CARACTERES, luego no hace falta expandir polinomios:
# basta evaluar en puntos del toro.  El mas barato es z = 1, que da la dimension de Weyl
#
#     dim(lambda) = prod_{i<j} (lambda_i - lambda_j + j - i) / (j - i)
#
# y el invariante entero
#
#     Delta(beta)  :=  sum_{g in G2} w(g) * dim(atil(T_g)) * dim(astar(T_g))
#
#     Delta != 0   ==>   [Phi]_{D2} != 0   ==>   Phi_t != 0.
#
# Es SANO por construccion (un caracter de una representacion nula es nulo).  Que ademas sea COMPLETO
# es una medida, no un teorema: en las 132 formas de la poblacion objetivo del rango anterior,
# Delta = 0 exactamente cuando el polinomio es 0 (128 / 128 y 4 / 4).
#
# POR QUE IMPORTA.  Convierte 24 formas que estaban MEDIDAS en 24 CERTIFICADOS -- un entero no nulo
# es una prueba, no una observacion -- y quita la expansion polinomica, que era el cuello de botella.
# Con eso el barrido llega a configuraciones que antes no se podian tocar.
#
# LA PREGUNTA QUE ESTE BARRIDO CONTESTA.  En el rango anterior los 4 supervivientes estaban TODOS en
# t=6, r=3 -- la configuracion mas grande -- y t=8, r=2 no producia ni una sola forma objetivo.  Eso
# es la firma de un rango corto, no de un conjunto excepcional finito.  Aqui se anaden t=8 r=3,
# t=6 r=4 y t=10 r=2, que no se habian tocado nunca, y se sube M en las viejas.
#
#   Z1  poblacion objetivo y supervivientes POR CONFIGURACION.  Si los supervivientes siguen
#       apareciendo solo en la config mas grande, el rango sigue siendo el que manda.
#   Z2  |G2| de cada superviviente.  En el rango anterior todos tenian |G2| = 2, o sea caian en la
#       Proposicion del segundo estrato.  Si aparece uno con |G2| >= 3, ESO si seria nuevo.
#   Z3  no vacuidad: se imprime SIEMPRE el n de cada config, incluidas las que dan 0.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python dim_certificate.py

import itertools
from collections import Counter
from fractions import Fraction

from second_stratum import setup, all_transversals, inv_of

CONFIGS = [(4, 2, 15), (4, 3, 15), (6, 2, 17), (6, 3, 18),
           (8, 2, 18), (8, 3, 19), (6, 4, 19), (10, 2, 19)]


def dim_gl(lam):
    r = len(lam)
    d = Fraction(1)
    for i in range(r):
        for j in range(i + 1, r):
            d *= Fraction(lam[i] - lam[j] + j - i, j - i)
    return int(d)


def halves(T, r):
    """(atil, astar) del diccionario a PvW."""
    H, L = T[:r], T[r:]
    alpha = [H[i] - (r - 1 - i) for i in range(r)]
    atil = tuple(a - alpha[-1] for a in alpha)
    Ls = [L[0] - L[r - 1 - i] for i in range(r)]
    astar = tuple(Ls[i] - (r - 1 - i) for i in range(r))
    return atil, astar


def main():
    shapes = Counter()
    targ = Counter()
    surv = Counter()
    surv_list = []
    g2_of_surv = Counter()

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
            shapes[(t, r)] += 1
            tr = all_transversals(beta, cl, r, t)
            D = max(x[3] for x in tr)
            G = [x for x in tr if x[3] == D]
            if len(G) != 2:
                continue
            a, b = G
            if not (inv_of(a[1], r) == inv_of(b[1], r) and a[2] == -b[2]):
                continue                              # [Phi]_top != 0
            S = sorted({v for k in E for v in Cd[k]})
            C = S[0] + S[-1]
            if set(C - v for v in S) == set(S):
                continue                              # (i) cierta
            rest = [x for x in tr if x[3] < D]
            if not rest:
                continue
            targ[(t, r)] += 1
            D2 = max(x[3] for x in rest)
            G2 = [x for x in rest if x[3] == D2]
            delta = 0
            for (_, T, w, _) in G2:
                at, ast = halves(T, r)
                delta += w * dim_gl(at) * dim_gl(ast)
            if delta == 0:
                surv[(t, r)] += 1
                g2_of_surv[len(G2)] += 1
                if len(surv_list) < 30:
                    surv_list.append((t, r, beta, len(G2)))
        print("   hecho t=%d r=%d M=%d : %d formas, %d objetivo, %d supervivientes"
              % (t, r, M, shapes[(t, r)], targ.get((t, r), 0), surv.get((t, r), 0)), flush=True)

    print("")
    print("Z1  %-12s %10s %10s %14s" % ("config", "formas", "objetivo", "supervivientes"))
    for k in sorted(shapes):
        print("    t=%-2d r=%-2d   %10d %10d %14d"
              % (k[0], k[1], shapes[k], targ.get(k, 0), surv.get(k, 0)))
    print("")
    print("Z3 totales: objetivo %d, supervivientes %d" % (sum(targ.values()), sum(surv.values())))
    print("Z2 |G2| de los supervivientes: %s" % dict(sorted(g2_of_surv.items())))
    print("")
    for s in surv_list:
        print("    t=%d r=%d |G2|=%d beta=%s" % (s[0], s[1], s[3], list(s[2])))


if __name__ == "__main__":
    main()
