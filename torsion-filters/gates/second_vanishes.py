# -*- coding: utf-8 -*-
# SE ANULA DE VERDAD EL SEGUNDO ESTRATO?  La pregunta que decide el alcance de la seccion 7.
#
# EL ERROR QUE ESTO CORRIGE.  Al escribir la Observacion 7.5 conte como ABIERTAS las formas con
# |G2| >= 3, porque Purbhoo-van Willigenburg es una rigidez ENTRE DOS productos y no hay teorema
# para tres.  Pero "no hay teorema que las cubra" no es lo mismo que "no se resuelven": si
# [Phi]_{D2} no es cero, la forma esta resuelta y no hace falta ninguna rigidez.  Conte lo que
# APARECE sin teorema en vez de lo que SOBREVIVE a la medida.
#
# Y hay una razon estructural para no buscar ese teorema: la independencia lineal de los productos
# {s_lambda s_mu} es FALSA por recuento de dimensiones -- en grado 6 hay sum_k p(k)p(6-k) = 65
# productos en un espacio de dimension p(6) = 11.  Las relaciones lineales son abundantes.  Lo
# rigido en PvW es la IGUALDAD DE DOS productos, no una relacion lineal cualquiera.  Un analogo
# para tres terminos no existe en general.
#
# QUE SE CALCULA.  Para cada forma de la POBLACION OBJETIVO ([Phi]_top = 0 y (i) falsa) se construye
# explicitamente
#       [Phi]_{D2}  =  sum_{g in G2} w(g) * P(T_g),      P(T) = a_H(z) * a_L(1/z)
# como polinomio de Laurent con coeficientes enteros (expansion de los dos alternantes r x r), y se
# comprueba si es cero.  No se agrupa por INV: eso solo veria cancelaciones DENTRO de una clase, y
# la pregunta es justamente si clases distintas pueden cancelarse entre si.
#
# LO QUE SE MIDE
#   W1  tabla (|G2|, [Phi]_{D2} = 0) con las cuatro esquinas.
#   W2  cuantas formas quedan RESUELTAS (estrato no nulo => Phi_t != 0).
#   W3  las supervivientes, listadas con su beta -- son las que hay que atacar.
#   W4  no vacuidad: n de la poblacion objetivo se imprime SIEMPRE.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python second_vanishes.py

import itertools
from collections import defaultdict, Counter

from second_stratum import setup, all_transversals, inv_of, perm_sign

CONFIGS = [(4, 2, 15), (4, 3, 15), (6, 2, 15), (6, 3, 17), (8, 2, 17)]


def alternant(X, inverse):
    """det(z_j^{X_i}) o det(z_j^{-X_i}) como dict exponente -> entero."""
    r = len(X)
    out = defaultdict(int)
    for s in itertools.permutations(range(r)):
        e = tuple((-X[s[j]] if inverse else X[s[j]]) for j in range(r))
        out[e] += perm_sign(s)
    return {k: v for k, v in out.items() if v}


def P_poly(T, r):
    """P(T) = a_H(z) * a_L(1/z)."""
    A = alternant(T[:r], False)
    B = alternant(T[r:], True)
    out = defaultdict(int)
    for ea, va in A.items():
        for eb, vb in B.items():
            out[tuple(x + y for x, y in zip(ea, eb))] += va * vb
    return {k: v for k, v in out.items() if v}


def main():
    tot = 0
    table = Counter()
    survivors = []

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
            if len(G) != 2:
                continue
            a, b = G
            if not (inv_of(a[1], r) == inv_of(b[1], r) and a[2] == -b[2]):
                continue                      # [Phi]_top != 0
            S = sorted({v for k in E for v in Cd[k]})
            C = S[0] + S[-1]
            if set(C - v for v in S) == set(S):
                continue                      # (i) cierta
            rest = [x for x in tr if x[3] < D]
            if not rest:
                continue
            tot += 1
            D2 = max(x[3] for x in rest)
            G2 = [x for x in rest if x[3] == D2]
            acc = defaultdict(int)
            for (_, T, w, _) in G2:
                for k, v in P_poly(T, r).items():
                    acc[k] += w * v
            zero = not any(acc.values())
            table[(len(G2), zero)] += 1
            if zero:
                survivors.append((t, r, beta, len(G2)))

    print("W4 POBLACION OBJETIVO ([Phi]_top = 0 y (i) falsa) : %d" % tot)
    print("")
    print("W1   |G2|   [Phi]_{D2} = 0 ?   formas")
    for k in sorted(table):
        print("      %3d      %-3s             %d" % (k[0], "SI" if k[1] else "NO", table[k]))
    print("")
    resolved = sum(v for k, v in table.items() if not k[1])
    proved = table.get((1, False), 0)
    print("W2 formas con [Phi]_{D2} != 0, luego Phi_t != 0 : %d de %d" % (resolved, tot))
    print("     de ellas PROBADAS por el Corolario del segundo maximizador unico (|G2| = 1) : %d"
          % proved)
    print("     las otras %d tienen |G2| >= 2: el estrato no se anula, pero eso esta MEDIDO,"
          % (resolved - proved))
    print("     no probado -- no hay rigidez para tres o mas productos, ni puede haberla en general.")
    print("")
    print("W3 SUPERVIVIENTES (el segundo estrato tambien se anula) : %d" % len(survivors))
    for s in survivors:
        print("      t=%d r=%d |G2|=%d beta=%s" % (s[0], s[1], s[3], list(s[2])))
    if survivors and all(s[3] == 2 for s in survivors):
        print("   TODOS tienen |G2| = 2, luego caen en la Proposicion del segundo estrato (PvW).")


if __name__ == "__main__":
    main()
