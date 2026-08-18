# -*- coding: utf-8 -*-
# EL SEGUNDO ESTRATO.  Arriba esta probado; el de abajo NO es separable (H7 murio en Gbot_anatomy).
# Esta es la tercera puerta: el escalon INMEDIATAMENTE siguiente a Dmax.
#
# POR QUE AHI Y NO ABAJO.  De la bitacora del 12: "en las 20 excepciones el primer grado que
# sobrevive esta 2 por debajo de Dmax en 17 casos y 4 en los otros 3".  La cancelacion NO es
# profunda.  El estrato de abajo basta empiricamente pero es el objeto mas hostil que tenemos
# (permanente, |Gbot| hasta 16, sin separabilidad, sin teorema de rigidez).  El de Dmax-g sigue
# siendo un ALTERNANTE, y si esta limpio, PvW Thm 2.5 vuelve a aplicar tal cual.
#
# LA DERIVACION QUE SE TESTA AQUI (a mano sobre beta=[10,9,7,4,3,2,1,0], t=4, r=2: sale 2, y la
# clase que se mueve es la 2, que es de las NO empatadas -- o sea g_com).
#
#   HUECO EXTERNO   g_ext = Dmax - max{deg(T_g) : g no maximizador}
#                   Es el menor  d_i - d_j > 0  con i <= r < j sobre los incrementos ordenados.
#   HUECO INTERNO   g_int = min sobre maximizadores de  2*(u_{r-1} - u_r)
#                   Segundo grado DENTRO de un mismo A(T): mover u_{r-1} de H a L cuesta eso.
#   LIMPIEZA        g_ext < g_int  =>  el estrato Dmax-g_ext no recibe nada de los maximizadores,
#                   luego es una suma de P(T) = a_H(z)a_L(1/z) puros y el argumento de la seccion 5
#                   de note_t2 se repite un piso mas abajo.
#
# H10  el hueco externo cumple  g_ext = 2*(k_out - k_in)  (mod t),  con k_out la clase del
#      incremento que sale y k_in la del que entra.  Analogo exacto de Delta_k(j) = 2k (mod t).
# H11  g_ext < g_int  (el segundo estrato esta LIMPIO).
# H12  el intercambio del segundo estrato TOCA g_com, o sea mueve alguna clase fuera de la pareja
#      empatada {k, k+t/2}.  Prediccion fina: lo hace exactamente cuando  g_ext/2  no es  0 mod t/2.
#      Con g_ext = 2 y t >= 4 siempre lo hace; con g_ext = 4 y t = 4 nunca.  Eso explicaria los 3
#      casos raros de los 20.
#
# CONTROLES QUE PUEDEN FALLAR
#   C1  senuelo de congruencia:  g_ext = (k_out - k_in) mod t  (sin el 2).  Debe fallar mucho.
#   C2  senuelo de congruencia:  g_ext = 2*(k_out + k_in) mod t.  Debe fallar mucho.
#   C3  no vacuidad: hay que contar cuantas formas tienen |G| = 2 y e > 2 -- el caso DURO, el unico
#       donde g_com no es vacio.  Si sale 0, la medida no dice nada y hay que decirlo.
#   C4  el numero de clases en que difieren maximizador y segundo estrato.  La derivacion predice
#       EXACTAMENTE 2.  Si sale otra cosa, la reformulacion de un-entero-por-clase no captura el
#       segundo estrato y H10 cae con ella.
#   C5  la poblacion objetivo: formas con [Phi]_top = 0 y (i) falsa.  Son las que el estrato de
#       arriba no ve.  Toda la brecha de t >= 4 estaba caracterizada por ~20 de esas; aqui se
#       engordan y se re-mide en ellas la profundidad.  n se imprime SIEMPRE.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python second_stratum.py   (aritmetica entera pura, no hace falta Sage)

import itertools
import sys
from collections import defaultdict


def perm_sign(seq):
    """signo de la permutacion que ordena seq (elementos distintos)."""
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


def setup(beta, t):
    """clases de residuos; None si alguna esta vacia (hipotesis de ocupacion (O))."""
    cl = defaultdict(list)
    for i, b in enumerate(beta):
        cl[b % t].append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    C = {k: sorted((beta[i] for i in cl[k]), reverse=True) for k in E}
    return dict(cl), E, C


def deg_of(T, r):
    return sum(T[:r]) - sum(T[r:])


def inv_of(T, r):
    """INV(T) del diccionario a PvW: (h_r - l_1, {atil, astar} como MULTICONJUNTO).
    Dos T dan el mismo P(T) = a_H(z)a_L(1/z) exactamente cuando coincide INV (control D2)."""
    H = T[:r]
    L = T[r:]
    alpha = [H[i] - (r - 1 - i) for i in range(r)]
    atil = tuple(a - alpha[-1] for a in alpha)
    Lstar = [L[0] - L[r - 1 - i] for i in range(r)]
    astar = tuple(Lstar[i] - (r - 1 - i) for i in range(r))
    return (H[-1] - L[0], tuple(sorted([atil, astar])))


def all_transversals(beta, cl, r, t):
    """[(clave por clase, T, w, deg)] sobre todas las transversales."""
    N = len(beta)
    out = []
    keys = sorted(cl)
    for pick in itertools.product(*[cl[k] for k in keys]):
        P = sorted(pick)
        Ps = set(P)
        T = tuple(beta[i] for i in range(N) if i not in Ps)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        sel = {k: beta[i] for k, i in zip(keys, pick)}
        out.append((sel, T, w, deg_of(T, r)))
    return out


def analyse(beta, t, r):
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(beta, cl, r, t)
    Dmax = max(x[3] for x in tr)
    G = [x for x in tr if x[3] == Dmax]
    rest = [x for x in tr if x[3] < Dmax]
    if not rest:
        return None
    D2 = max(x[3] for x in rest)
    G2 = [x for x in rest if x[3] == D2]
    g_ext = Dmax - D2
    g_int = min(2 * (x[1][r - 1] - x[1][r]) for x in G)

    # [Phi]_top = 0  <=>  |G| = 2, mismo INV, y signos opuestos.
    top_zero = False
    if len(G) == 2:
        a, b = G
        top_zero = (inv_of(a[1], r) == inv_of(b[1], r)) and (a[2] == -b[2])

    # (i): S concentrico.  S = valores de exceso, C = min S + max S.
    S = sorted({v for k in E for v in Cd[k]})
    Cc = S[0] + S[-1]
    cond_i = set(Cc - v for v in S) == set(S)

    # las clases empatadas: en las que difieren los dos maximizadores.
    tied = set()
    if len(G) == 2:
        tied = {k for k in sorted(cl) if G[0][0][k] != G[1][0][k]}

    # el paso al segundo estrato: en que clases difiere, y de que clase sale / entra el incremento.
    swaps = []
    for h in G2:
        for m in G:
            diff = sorted(k for k in sorted(cl) if h[0][k] != m[0][k])
            if len(diff) == 2:
                # el que BAJA en la clase (elige un elemento mas pequeno) gano un incremento en A;
                # orientacion: k_in recibe, k_out cede.
                k1, k2 = diff
                if h[0][k1] < m[0][k1]:
                    k_in, k_out = k1, k2
                else:
                    k_in, k_out = k2, k1
                swaps.append((len(diff), k_out, k_in))
                break
        else:
            best = min(
                (sorted(k for k in sorted(cl) if h[0][k] != m[0][k]) for m in G), key=len)
            swaps.append((len(best), None, None))
    return dict(Dmax=Dmax, g_ext=g_ext, g_int=g_int, nG=len(G), nG2=len(G2),
                e=len(E), top_zero=top_zero, cond_i=cond_i, tied=tied,
                swaps=swaps, cl=cl)


CONFIGS = [(4, 2, 15), (4, 3, 15), (6, 2, 15), (6, 3, 17), (8, 2, 17), (8, 3, 19)]
CAP = 40000


def main():
    tot = 0
    hard = 0                      # C3: |G| = 2 y e > 2
    h10_bad = h10_n = 0
    c1_bad = c1_n = 0
    c2_bad = c2_n = 0
    h11_bad = h11_n = 0
    h12_touch = h12_no = 0
    c4 = defaultdict(int)
    target = []                   # C5: top_zero y (i) falsa
    depth = defaultdict(int)

    for (t, r, M) in CONFIGS:
        N = t + 2 * r
        if M < N - 1:
            continue
        n_cfg = 0
        for comb in itertools.combinations(range(M + 1), N):
            beta = tuple(sorted(comb, reverse=True))
            a = analyse(beta, t, r)
            if a is None:
                continue
            n_cfg += 1
            tot += 1
            if a['nG'] == 2 and a['e'] > 2:
                hard += 1
            for (nd, k_out, k_in) in a['swaps']:
                c4[nd] += 1
                if k_out is None:
                    continue
                h10_n += 1
                if (a['g_ext'] - 2 * (k_out - k_in)) % t != 0:
                    h10_bad += 1
                c1_n += 1
                if (a['g_ext'] - (k_out - k_in)) % t != 0:
                    c1_bad += 1
                c2_n += 1
                if (a['g_ext'] - 2 * (k_out + k_in)) % t != 0:
                    c2_bad += 1
                if a['tied']:
                    if {k_out, k_in} - a['tied']:
                        h12_touch += 1
                    else:
                        h12_no += 1
            h11_n += 1
            if not (a['g_ext'] < a['g_int']):
                h11_bad += 1
            depth[a['g_ext']] += 1
            if a['top_zero'] and not a['cond_i']:
                target.append((t, r, beta, a['g_ext'], a['g_int'],
                               a['e'], a['nG2'],
                               any(k_out is not None and bool({k_out, k_in} - a['tied'])
                                   for (_, k_out, k_in) in a['swaps'])))
            if n_cfg >= CAP:
                break
        print("  cfg t=%d r=%d M=%d : %d formas" % (t, r, M, n_cfg))
        sys.stdout.flush()

    print("")
    print("FORMAS TOTALES (ocupacion (O) y exceso no vacio): %d" % tot)
    print("C3 no vacuidad -- |G| = 2 Y e > 2 (el caso DURO, g_com no vacio): %d" % hard)
    print("")
    print("C4  clases en que difiere el segundo estrato de un maximizador (prediccion: 2)")
    for k in sorted(c4):
        print("      %s clases : %d" % (k, c4[k]))
    print("")
    print("H10 g_ext = 2*(k_out - k_in) mod t : %d fallos de %d" % (h10_bad, h10_n))
    print("C1  SENUELO  g_ext = (k_out - k_in) mod t : %d fallos de %d" % (c1_bad, c1_n))
    print("C2  SENUELO  g_ext = 2*(k_out + k_in) mod t : %d fallos de %d" % (c2_bad, c2_n))
    print("")
    print("H11 g_ext < g_int (segundo estrato LIMPIO) : %d fallos de %d" % (h11_bad, h11_n))
    print("H12 el intercambio sale de la pareja empatada : %d si, %d no" % (h12_touch, h12_no))
    print("")
    print("profundidad del segundo estrato (Dmax - D2):")
    for k in sorted(depth):
        print("      %s : %d" % (k, depth[k]))
    print("")
    print("C5 POBLACION OBJETIVO -- [Phi]_top = 0 y (i) FALSA : n = %d" % len(target))
    if target:
        cl_ = sum(1 for x in target if x[3] < x[4])
        tc = sum(1 for x in target if x[7])
        print("      de esas, segundo estrato limpio (g_ext < g_int) : %d" % cl_)
        print("      de esas, el intercambio toca g_com              : %d" % tc)
        dd = defaultdict(int)
        for x in target:
            dd[x[3]] += 1
        print("      profundidad en la poblacion objetivo:")
        for k in sorted(dd):
            print("          %s : %d" % (k, dd[k]))
        print("      primeras 8:")
        for x in target[:8]:
            print("          t=%d r=%d beta=%s g_ext=%d g_int=%d e=%d |G2|=%d toca_gcom=%s"
                  % (x[0], x[1], list(x[2]), x[3], x[4], x[5], x[6], x[7]))


if __name__ == "__main__":
    main()
