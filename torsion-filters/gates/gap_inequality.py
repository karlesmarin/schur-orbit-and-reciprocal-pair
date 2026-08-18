# -*- coding: utf-8 -*-
# LA DESIGUALDAD  g_ext < g_int,  ahora con prueba y no solo medida.
#
# g_int(T) = 2*(u_{r-1} - u_r)   caida INTERNA de A(T): mover u_{r-1} de H a L.
# g_ext    = Dmax - D2           salto EXTERNO al siguiente grado sobre TODAS las transversales.
# Si g_ext < g_int para todo maximizador, el estrato Dmax - g_ext NO recibe nada de los
# maximizadores y es suma pura de P(T) = a_H(z)a_L(1/z): PvW Thm 2.5 aplica un piso mas abajo.
#
# EL MOVIMIENTO.  Un error mio, y es el que abre la prueba: crei que las competidoras venian de
# INTERCAMBIAR incrementos entre clases.  Falso -- la forma de prefijo caracteriza los OPTIMOS, no
# las demas transversales.  El movimiento bueno es desplazar UN elemento de la transversal un solo
# paso dentro de su clase:
#
#   SUBIDA en k* = clase de u_{r-1} = min H.   g_{k*} = c_{k*,j+1}  ->  c_{k*,j} = u_{r-1}.
#       si g_{k*} > u_r :  H'' = H - u_{r-1} + g_{k*},  L'' = L,  delta = u_{r-1} - g_{k*}
#                          y  delta < u_{r-1} - u_r < g_int.
#       si g_{k*} < u_r :  delta = u_{r-1} + g_{k*} - 2*u_r,  y  delta < g_int  <=>  g_{k*} < u_{r-1}. OK
#   BAJADA en k'* = clase de u_r = max L.      g_{k'*} = c_{k',j'+1}  ->  c_{k',j'+2} = u_r.
#       si g_{k'*} < u_{r-1} : delta = g_{k'*} - u_r  <  u_{r-1} - u_r  <  g_int.
#       si g_{k'*} > u_{r-1} : delta = 2*u_{r-1} - g_{k'*} - u_r,  y  delta < g_int  <=>  g_{k'*} > u_r. OK
#   CASO DEGENERADO k* = k'*  (u_{r-1}, g_k, u_r consecutivos en la MISMA clase): lo cubre la
#       bajada, con delta = g_k - u_r.
#
# TEOREMA, mas fuerte de lo que escribi la primera vez:   g_ext  <  u_{r-1} - u_r  =  g_int / 2.
# En CADA UNO de los cuatro casos la desigualdad  delta < u_{r-1} - u_r  ES la hipotesis de su
# propio caso:  g_{k*} > u_r ,  g_{k*} < u_r ,  g_{k'*} < u_{r-1} ,  g_{k'*} > u_{r-1}.  O sea que
# el estrato Dmax - g_ext esta limpio con FACTOR 2 de margen, no por los pelos.
# (Lo escribi primero como "solo en el caso degenerado" y era falso: vale en los cuatro.  Lo cazo
#  P5, que puse como senuelo esperando que fallara y salio 157386 de 157386 -- ver la correccion en
#  gap_inequality_OUT.txt: NO es un control, es un teorema que no habia visto.)
#
# delta >= 0 por maximalidad.  LO UNICO que puede romper el argumento es delta = 0 en los DOS
# movimientos a la vez: entonces no hay competidor ESTRICTO y la cota no produce nada.  Eso es P2.
#
# LO QUE SE MIDE, cada cosa capaz de fallar
#   P1  0 <= delta < g_int para los dos movimientos, en TODO maximizador.  Es la prueba de arriba;
#       si falla, el analisis de casos esta mal.
#   P2  cuantas veces los DOS movimientos dan delta = 0 -- el unico agujero del argumento.
#   P3  g_ext <= min(delta positivos).  La prueba solo da una COTA; esto comprueba que la cota es
#       correcta y de paso cuanto se pierde.
#   P4  g_ext < g_int, la conclusion, medida directamente.
#   P5  SENUELO que tiene que fallar: "delta < g_int/2", que solo vale en el caso degenerado.
#   P6  no vacuidad: cuantas formas tienen k* = k'* (el caso degenerado) y cuantas |G| = 2.
#   P7  RANGO ANCHO.  El barrido de second_stratum.py llegaba a M = 19 y el caso degenerado que
#       construi a mano necesita beta hasta 24.  Aqui se anaden clases DISPERSAS a proposito
#       (pasos grandes dentro de una clase), que es donde vive el contraejemplo si existe.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python gap_inequality.py

import itertools
import random
from collections import defaultdict

from second_stratum import setup, all_transversals


def moves(beta, t, r):
    """(g_int, [deltas de los dos movimientos], g_ext, k*==k'*, |G|) o None."""
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
    g_ext = Dmax - max(x[3] for x in rest)

    out = []
    for (sel, T, w, d) in G:
        H, L = T[:r], T[r:]
        u1, u2 = H[-1], L[0]
        g_int = 2 * (u1 - u2)
        kstar = u1 % t
        kpstar = u2 % t
        deltas = []
        for (k, target) in ((kstar, 'up'), (kpstar, 'down')):
            if k not in Cd:
                continue
            ck = Cd[k]
            cur = sel[k]
            i = ck.index(cur)
            j = i - 1 if target == 'up' else i + 1
            if not (0 <= j < len(ck)):
                continue
            newsel = dict(sel)
            newsel[k] = ck[j]
            Tn = tuple(sorted((v for v in beta if v not in set(newsel.values())),
                              reverse=True))
            if len(Tn) != 2 * r:
                continue
            deltas.append(Dmax - (sum(Tn[:r]) - sum(Tn[r:])))
        out.append((g_int, deltas, kstar == kpstar))
    return out, g_ext, len(G)


def sweep(shapes, tag, acc):
    for (beta, t, r) in shapes:
        res = moves(beta, t, r)
        if res is None:
            continue
        per, g_ext, nG = res
        acc['n'] += 1
        if nG == 2:
            acc['nG2'] += 1
        best = None
        for (g_int, deltas, degen) in per:
            if degen:
                acc['degen'] += 1
            acc['gmin'] = min(acc.get('gmin', 10 ** 9), g_int)
            for dl in deltas:
                acc['p1_n'] += 1
                if not (0 <= dl < g_int):
                    acc['p1_bad'] += 1
                    acc['p1_ex'].append((t, r, beta, dl, g_int))
                if dl < g_int / 2.0:
                    acc['p5_half'] += 1
                else:
                    acc['p5_not'] += 1
                if dl > 0:
                    best = dl if best is None else min(best, dl)
            if deltas and max(deltas) == 0:
                acc['p2'] += 1
                acc['p2_ex'].append((t, r, beta))
        if best is not None:
            acc['p3_n'] += 1
            if g_ext > best:
                acc['p3_bad'] += 1
                acc['p3_ex'].append((t, r, beta, g_ext, best))
            acc['slack'][g_ext - best] += 1
        acc['p4_n'] += 1
        if not all(g_ext < g_int for (g_int, _, _) in per):
            acc['p4_bad'] += 1
            acc['p4_ex'].append((t, r, beta, g_ext, [x[0] for x in per]))


def main():
    acc = defaultdict(int)
    acc['p1_ex'] = []
    acc['p2_ex'] = []
    acc['p3_ex'] = []
    acc['p4_ex'] = []
    acc['slack'] = defaultdict(int)

    # A) exhaustivo, rango moderado
    shapes = []
    for (t, r, M) in [(2, 2, 12), (4, 2, 14), (4, 3, 14), (6, 2, 14), (8, 2, 16)]:
        N = t + 2 * r
        for comb in itertools.combinations(range(M + 1), N):
            shapes.append((tuple(sorted(comb, reverse=True)), t, r))
    sweep(shapes, 'exhaustivo', acc)
    print("A) exhaustivo: %d formas validas" % acc['n'])

    # B) P7: clases DISPERSAS a proposito -- se construye beta por clases con pasos grandes,
    #    que es el regimen donde vive el caso degenerado (mi testigo a mano llegaba a 24).
    rnd = random.Random(20260812)
    shapes = []
    for _ in range(60000):
        t = rnd.choice([2, 4, 6])
        r = rnd.choice([2, 3])
        N = t + 2 * r
        sizes = [1] * t
        for _ in range(2 * r):
            sizes[rnd.randrange(t)] += 1
        beta = []
        ok = True
        for k in range(t):
            base = rnd.randrange(0, 6)
            vals = []
            cur = k + t * base
            for _ in range(sizes[k]):
                vals.append(cur)
                cur += t * rnd.choice([1, 1, 2, 3, 5, 8])
            beta.extend(vals)
        if len(set(beta)) != N:
            continue
        shapes.append((tuple(sorted(beta, reverse=True)), t, r))
    before = acc['n']
    sweep(shapes, 'disperso', acc)
    print("B) disperso (P7): %d formas validas mas" % (acc['n'] - before))

    print("")
    print("P6 no vacuidad: caso degenerado k* = k'* : %d ;  formas con |G| = 2 : %d"
          % (acc['degen'], acc['nG2']))
    print("")
    print("P1  0 <= delta < g_int en los dos movimientos : %d fallos de %d"
          % (acc['p1_bad'], acc['p1_n']))
    for e in acc['p1_ex'][:5]:
        print("      CONTRAEJEMPLO t=%d r=%d beta=%s delta=%d g_int=%d" % e)
    print("P2  los DOS movimientos dan delta = 0 (el agujero) : %d" % acc['p2'])
    for e in acc['p2_ex'][:5]:
        print("      t=%d r=%d beta=%s" % e)
    print("P3  g_ext <= min(delta positivos) : %d fallos de %d" % (acc['p3_bad'], acc['p3_n']))
    for e in acc['p3_ex'][:5]:
        print("      t=%d r=%d beta=%s g_ext=%d cota=%d" % e)
    print("P4  g_ext < g_int (LA CONCLUSION) : %d fallos de %d" % (acc['p4_bad'], acc['p4_n']))
    for e in acc['p4_ex'][:5]:
        print("      t=%d r=%d beta=%s g_ext=%d g_int=%s" % e)
    print("P5  SENUELO 'delta < g_int/2' : %d cumplen, %d NO (tiene que haber de los dos)"
          % (acc['p5_half'], acc['p5_not']))
    print("")
    print("holgura de la cota  g_ext - min(delta positivos)  (0 = la cota es exacta):")
    for k in sorted(acc['slack'])[:10]:
        print("      %s : %d" % (k, acc['slack'][k]))


if __name__ == "__main__":
    main()
