# -*- coding: utf-8 -*-
# CIERRE DE LOS DOS CABOS QUE ABRIO second_stratum.py.
#
# R1  El salto de UNA sola clase.  Mi C4 predecia 2 clases y en 65172 casos es 1: bajar al
#     incremento siguiente DENTRO de la misma clase.  Esos ni se testaron en H10.  Con k_out = k_in
#     la congruencia degenera a  g_ext = 0 (mod t).  Se testa, con su senuelo.
#
# R2  H12 en su forma FINA.  El intercambio sale de la pareja empatada {k, k+t/2} exactamente
#     cuando  (g_ext / 2)  no es  0 (mod t/2).  Es lo que explicaria las 12 formas de la poblacion
#     objetivo que NO tocan g_com, y los 3 casos raros de los 20 originales.  Tabla de contingencia
#     completa, con las cuatro casillas -- no solo la diagonal.
#
# R3  La rama de TRASLACION, que es la deuda declarada del paso 5 en t = 2.  Derivado hoy:
#     si T_B = T_A + m entonces u = v = r y  T_A u T_B  son DOS progresiones aritmeticas del mismo
#     paso y longitud r+1.  En t = 2 eso es beta ENTERO.  Aqui se busca un beta asi que ademas sea
#     forma anulante.  Si no existe en el rango, la rama se puede borrar del caso t = 2 y la deuda
#     muere; si existe, hay que dejarla escrita.  n se imprime SIEMPRE.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python second_stratum_refine.py

import itertools
from collections import defaultdict

from second_stratum import setup, deg_of, inv_of, all_transversals

CONFIGS = [(4, 2, 15), (4, 3, 15), (6, 2, 15), (6, 3, 16), (8, 2, 17), (2, 2, 13),
           (2, 3, 15), (2, 4, 17)]


def two_aps(vals):
    """True si el conjunto es union de dos progresiones del MISMO paso y la MISMA longitud."""
    v = sorted(vals)
    n = len(v)
    if n % 2 or n < 4:
        return False
    half = n // 2
    for split in itertools.combinations(range(1, n), half - 1):
        pass
    # exhaustivo por particiones en dos mitades del mismo tamano
    idx = set(range(n))
    for A in itertools.combinations(range(n), half):
        B = sorted(idx - set(A))
        a = [v[i] for i in A]
        b = [v[i] for i in B]
        da = {a[i + 1] - a[i] for i in range(half - 1)}
        db = {b[i + 1] - b[i] for i in range(half - 1)}
        if len(da) == 1 and da == db:
            return True
    return False


def main():
    r1_bad = r1_n = 0
    r1_decoy = 0
    cont = defaultdict(int)          # (fuera_de_pareja, gext/2 != 0 mod t/2) -> n
    ap_hits = []
    ap_n = 0
    tot = 0

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
            Dmax = max(x[3] for x in tr)
            G = [x for x in tr if x[3] == Dmax]
            rest = [x for x in tr if x[3] < Dmax]
            if not rest:
                continue
            tot += 1
            D2 = max(x[3] for x in rest)
            G2 = [x for x in rest if x[3] == D2]
            g_ext = Dmax - D2
            tied = set()
            if len(G) == 2:
                tied = {k for k in sorted(cl) if G[0][0][k] != G[1][0][k]}

            for h in G2:
                diff = min((sorted(k for k in sorted(cl) if h[0][k] != m[0][k])
                            for m in G), key=len)
                if len(diff) == 1:
                    r1_n += 1
                    if g_ext % t != 0:
                        r1_bad += 1
                    if g_ext % t == 0 and (g_ext // 2) % max(1, t // 2) != 0:
                        r1_decoy += 1
                elif len(diff) == 2 and tied:
                    outside = bool(set(diff) - tied)
                    pred = (t > 2) and ((g_ext // 2) % (t // 2) != 0)
                    cont[(outside, pred)] += 1

            # R3: rama de traslacion.  Solo interesa donde HAY empate.
            if len(G) == 2:
                TA, TB = G[0][1], G[1][1]
                union = set(TA) | set(TB)
                if len(union) == 2 * r + 2:
                    ap_n += 1
                    if two_aps(union):
                        m_ok = None
                        for m in range(-M, M + 1):
                            if m and set(x + m for x in TA) == set(TB):
                                m_ok = m
                        ap_hits.append((t, r, beta, m_ok))

    print("FORMAS: %d" % tot)
    print("")
    print("R1  salto de UNA clase.  prediccion g_ext = 0 (mod t)")
    print("      fallos: %d de %d" % (r1_bad, r1_n))
    print("      SENUELO g_ext/2 = 0 (mod t/2) sobre los que pasan: %d violaciones" % r1_decoy)
    print("")
    print("R2  contingencia   (sale de la pareja empatada)  x  (prediccion g_ext/2 != 0 mod t/2)")
    print("      sale=SI  pred=SI : %d      <- acuerdo" % cont[(True, True)])
    print("      sale=NO  pred=NO : %d      <- acuerdo" % cont[(False, False)])
    print("      sale=SI  pred=NO : %d      <- DESACUERDO" % cont[(True, False)])
    print("      sale=NO  pred=SI : %d      <- DESACUERDO" % cont[(False, True)])
    print("")
    print("R3  rama de TRASLACION: union de los dos maximizadores con 2r+2 elementos")
    print("      candidatos examinados : %d" % ap_n)
    print("      son dos progresiones gemelas : %d" % len(ap_hits))
    real = [h for h in ap_hits if h[3] is not None]
    print("      y ademas T_B = T_A + m de verdad : %d" % len(real))
    for h in real[:10]:
        print("          t=%d r=%d beta=%s  m=%s" % (h[0], h[1], list(h[2]), h[3]))


if __name__ == "__main__":
    main()
