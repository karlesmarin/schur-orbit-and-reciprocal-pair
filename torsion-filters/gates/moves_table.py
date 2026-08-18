# -*- coding: utf-8 -*-
# LA TABLA DE CUATRO CASOS DEL LEMA DE LOS DOS MOVIMIENTOS, testada fila por fila.
#
# POR QUE.  gap_inequality.py verifico  0 <= delta < g_int  reconstruyendo la transversal y
# recalculando el grado.  Pero NUNCA comprobo que delta valga la FORMULA que la tabla imprime en
# cada caso.  La tabla se le ensena al lector como un hecho y no la testaba nadie.
# Ver [[the-column-nobody-tested]].
#
# LA TABLA, tal como aparece en la nota (u = T ordenado decreciente, H = u_0..u_{r-1}, L = el resto,
# k* = clase de u_{r-1} = min H, k' = clase de u_r = max L, g_k = el elemento de la transversal en
# la clase k):
#
#   caso                              delta
#   subida,  g_{k*}  > u_r            u_{r-1} - g_{k*}
#   subida,  g_{k*}  < u_r            u_{r-1} + g_{k*} - 2 u_r
#   bajada,  g_{k'}  < u_{r-1}        g_{k'} - u_r
#   bajada,  g_{k'}  > u_{r-1}        2 u_{r-1} - g_{k'} - u_r
#
# Los cuatro casos son EXHAUSTIVOS y disjuntos porque g_{k*} y g_{k'} no estan en T (son elementos
# de la transversal) y u_{r-1} > u_r, asi que ninguno puede ser igual a u_{r-1} ni a u_r.
#
# LO QUE SE MIDE
#   T1  delta REAL (recalculado desde la transversal nueva) == delta de la FORMULA de su caso.
#       0 fallos o la tabla esta mal.
#   T2  cobertura: las cuatro filas tienen que darse.  Si alguna sale 0 veces, esa fila esta
#       SIN TESTAR y hay que decirlo, no contarla como verificada.
#   T3  SENUELO: cruzar las formulas -- aplicar a cada caso la formula del caso hermano.  Tiene que
#       FALLAR.  Si acertara, la tabla no estaria diciendo nada.
#   T4  exhaustividad: ninguna forma puede caer fuera de los cuatro casos.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python moves_table.py

import itertools
from collections import Counter

from second_stratum import setup, all_transversals

CONFIGS = [(4, 2, 15), (4, 3, 15), (6, 2, 16), (6, 3, 17), (8, 2, 17), (2, 2, 13), (2, 3, 15)]


def deg(T, r):
    return sum(T[:r]) - sum(T[r:])


def main():
    t1_bad = t1_n = 0
    t3_bad = t3_n = 0
    cover = Counter()
    outside = 0
    bad_ex = []

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
            D1 = max(x[3] for x in tr)
            for (sel, T, w, d) in tr:
                if d != D1:
                    continue
                u1, u2 = T[r - 1], T[r]              # u_{r-1} = min H,  u_r = max L
                kst, kpr = u1 % t, u2 % t
                for (k, direction) in ((kst, 'up'), (kpr, 'down')):
                    if k not in Cd:
                        continue
                    ck = Cd[k]
                    cur = sel[k]
                    i = ck.index(cur)
                    j = i - 1 if direction == 'up' else i + 1
                    if not (0 <= j < len(ck)):
                        continue
                    newsel = dict(sel)
                    newsel[k] = ck[j]
                    chosen = set(newsel.values())
                    Tn = tuple(sorted((v for v in beta if v not in chosen), reverse=True))
                    if len(Tn) != 2 * r:
                        continue
                    delta_real = D1 - deg(Tn, r)

                    g = cur                           # el elemento de la transversal que se mueve
                    if direction == 'up':
                        if g > u2:
                            case, f = 'subida g>u_r', u1 - g
                            decoy = u1 + g - 2 * u2
                        elif g < u2:
                            case, f = 'subida g<u_r', u1 + g - 2 * u2
                            decoy = u1 - g
                        else:
                            outside += 1
                            continue
                    else:
                        if g < u1:
                            case, f = 'bajada g<u_{r-1}', g - u2
                            decoy = 2 * u1 - g - u2
                        elif g > u1:
                            case, f = 'bajada g>u_{r-1}', 2 * u1 - g - u2
                            decoy = g - u2
                        else:
                            outside += 1
                            continue

                    cover[case] += 1
                    t1_n += 1
                    if delta_real != f:
                        t1_bad += 1
                        if len(bad_ex) < 5:
                            bad_ex.append((t, r, beta, case, delta_real, f))
                    t3_n += 1
                    if delta_real == decoy:
                        t3_bad += 1
        print("   hecho t=%d r=%d M=%d" % (t, r, M), flush=True)

    print("")
    print("T1 delta REAL == formula de su fila : %d fallos de %d" % (t1_bad, t1_n))
    for e in bad_ex:
        print("      CONTRAEJEMPLO t=%d r=%d beta=%s  caso=%s  real=%d  formula=%d" % e)
    print("")
    print("T2 COBERTURA de las cuatro filas:")
    for c in ['subida g>u_r', 'subida g<u_r', 'bajada g<u_{r-1}', 'bajada g>u_{r-1}']:
        n = cover.get(c, 0)
        print("      %-20s : %-8d %s" % (c, n, "" if n else "*** SIN TESTAR ***"))
    print("")
    print("T3 SENUELO (formula del caso hermano) : coincide en %d de %d  -- tiene que ser BAJO"
          % (t3_bad, t3_n))
    print("T4 casos fuera de la tabla (g igual a u_{r-1} o u_r, imposible) : %d" % outside)


if __name__ == "__main__":
    main()
