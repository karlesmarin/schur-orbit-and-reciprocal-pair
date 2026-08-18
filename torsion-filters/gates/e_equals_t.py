# -*- coding: utf-8 -*-
# LA CONTINGENCIA e = t, ARCHIVADA.  Es el recibo que faltaba.
#
# POR QUE EXISTE ESTE FICHERO.  La nota cita "0 de 216 con e < t" y "12 de 196 con e = t" como la
# evidencia de que  e = t  es NECESARIO para sobrevivir al segundo estrato.  Esas dos cifras salieron
# de una ejecucion EN LINEA que nunca se archivo.  Al auditar la trazabilidad numero a numero, 216 y
# 196 aparecian en conj_crit_t2_OUT y middle_block_OUT por pura COINCIDENCIA numerica -- ficheros de
# otra medicion completamente distinta.  Eso es peor que no encontrar nada: parece una traza y no lo
# es.  Si esas cifras van a un paper publicado, necesitan recibo reproducible.
# Ver [[save-the-outputs-not-just-the-scripts]].
#
# MISMA POBLACION Y MISMAS CONFIGURACIONES que dim_certificate.py, para que 216 + 196 = 412 cuadre
# con su Z3.  Si no cuadra, una de las dos mediciones esta mal y hay que pararlo todo.
#
#   E1  contingencia completa (e == t) x (superviviente), las cuatro casillas.
#   E2  el total tiene que ser 412 y los supervivientes 12, identicos a dim_certificate_OUT.
#       Es un control cruzado entre dos guiones: si discrepan, algo esta mal.
#   E3  desglose por configuracion de las formas con e = t, para ver donde es POSIBLE y donde no.
#   E4  la barrera aritmetica: e = t exige t <= 2r, porque sum(n_k - 1) = 2r se reparte entre e
#       clases con n_k >= 2.  Se comprueba que ninguna forma con t > 2r tiene e = t.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python e_equals_t.py

import itertools
from collections import Counter

from second_stratum import setup, all_transversals, inv_of
from dim_certificate import dim_gl, halves, CONFIGS


def main():
    cont = Counter()
    bycfg = Counter()
    e_eq_t_when_impossible = 0
    tot = 0
    surv = 0

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
                continue                                   # [N_t]_top != 0
            S = sorted({v for k in E for v in Cd[k]})
            C = S[0] + S[-1]
            if set(C - v for v in S) == set(S):
                continue                                   # (i) cierta
            rest = [x for x in tr if x[3] < D]
            if not rest:
                continue
            tot += 1
            D2 = max(x[3] for x in rest)
            G2 = [x for x in rest if x[3] == D2]
            delta = 0
            for (_, T, w, _) in G2:
                at, ast = halves(T, r)
                delta += w * dim_gl(at) * dim_gl(ast)
            is_surv = (delta == 0)
            surv += is_surv
            et = (len(E) == t)
            cont[(et, is_surv)] += 1
            if et:
                bycfg[(t, r, is_surv)] += 1
                if t > 2 * r:
                    e_eq_t_when_impossible += 1
        print("   hecho t=%d r=%d M=%d (objetivo acumulado %d)" % (t, r, M, tot), flush=True)

    n_et = cont[(True, True)] + cont[(True, False)]
    n_lt = cont[(False, True)] + cont[(False, False)]

    print("")
    print("E1 CONTINGENCIA sobre la poblacion objetivo ([N_t]_top = 0 y (i) falsa):")
    print("      e==t    superviviente   formas")
    for k in sorted(cont):
        print("      %-7s %-15s %d" % (k[0], k[1], cont[k]))
    print("")
    print("   e = t : %d formas, de las cuales sobreviven %d" % (n_et, cont[(True, True)]))
    print("   e < t : %d formas, de las cuales sobreviven %d   <-- NECESIDAD si es 0"
          % (n_lt, cont[(False, True)]))
    print("")
    print("E2 CONTROL CRUZADO contra dim_certificate_OUT.txt (Z3: objetivo 412, supervivientes 12)")
    print("      total objetivo aqui : %d   %s" % (tot, "OK" if tot == 412 else "*** DISCREPA ***"))
    print("      supervivientes aqui : %d   %s" % (surv, "OK" if surv == 12 else "*** DISCREPA ***"))
    print("")
    print("E3 formas con e = t, por configuracion:")
    for k in sorted(bycfg):
        print("      t=%d r=%d superviviente=%-5s : %d" % (k[0], k[1], k[2], bycfg[k]))
    print("")
    print("E4 barrera aritmetica: e = t exige t <= 2r.  Formas con e = t y t > 2r : %d  %s"
          % (e_eq_t_when_impossible,
             "OK (ninguna, como debe)" if e_eq_t_when_impossible == 0 else "*** IMPOSIBLE ***"))


if __name__ == "__main__":
    main()
