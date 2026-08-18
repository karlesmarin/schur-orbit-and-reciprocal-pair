# -*- coding: utf-8 -*-
# LA CAPA: que le falta a la condicion necesaria para ser criterio.  14 de agosto de 2026.
#
# DONDE ESTAMOS.  peel_zero.py dejo esto: el locus de ceros es cerrado bajo pelado (705 ceros, 0
# fallos), pero el reciproco NO vale -- en t=4 r=2 hay 708 formas con C=tau y pelado nulo y solo 128
# se anulan.  La pregunta es que separa esas 128.
#
# EL CONCEPTO, ANTES DEL PROGRAMA.  Phi == 0 obliga a que TODA pieza graduada del numerador se anule,
# y en particular la de arriba.  Por (29) esa pieza es
#
#       [det]_{D1} = const * V * sum_{g in G} w(g) P(T_g),
#
# asi que con |G| = 2 se anula exactamente cuando P(T_a) = P(T_b) Y w(g_a) = -w(g_b): hacen falta
# LAS DOS COSAS.  La Corolario 8.31 solo extrae la primera mitad -- la reflexion, via 8.28 -- y de
# ahi saca C = tau.  El SIGNO no aparece en C = tau por ningun lado.  Ese es el sospechoso: mi
# predicado de ayer pedia reflexion y pelado y se dejaba fuera la condicion de signo.
#
# P(T_a) = P(T_b) se lee de inv_of() -- el invariante del diccionario a PvW, (h_r - l_1, {atil,
# astar}) -- sin expandir nada.  Luego los tres candidatos de abajo son TODOS del voraz: coste dos
# ordenaciones, no una expansion polinomica.
#
# LOS CANDIDATOS, ESCRITOS ANTES DE CORRER
#   P1  [Phi]_top = 0      |G|=2  &  inv_of(T_a) = inv_of(T_b)  &  w_a = -w_b     <- con el signo
#   P2  el pelado se anula a rango r-1                                            <- el de ayer
#   P3  P1 y P2                                                                   <- LA APUESTA
#   P4  g_com es tau-simetrico
#   P5  S entero es tau-simetrico   (la condicion de t=2; en t>=4 deberia ser demasiado fuerte)
#   P6  P1 y P4
#   P7  P1 en TODOS los pisos de la bandera hasta rango 1
#
# CONTROL.  Un predicado que no puede separar: "e = 2".  Si saliera exacto, el arnes no mide.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python layer_condition.py

import itertools
import sys

from second_stratum import setup, all_transversals, inv_of
from depth_histogram import measure
from flag import pelar

CFG = [(4, 2, 17), (6, 2, 17), (8, 2, 18), (4, 3, 15), (6, 3, 15), (10, 2, 19)]


def betas(t, r, W):
    N = t + 2 * r
    for mid in itertools.combinations(range(1, W + 1), N - 1):
        yield tuple(sorted(mid, reverse=True)) + (0,)


def anat(b, t, rp):
    """todo lo que hace falta, de una sola pasada del voraz.  None si falta una clase."""
    st = setup(b, t)
    if st is None:
        return None                       # rama (a): Phi == 0, pero no es poblacion de esta pregunta
    cl, E, Cd = st
    if not E:
        return None
    tr = all_transversals(b, cl, rp, t)
    D = max(x[3] for x in tr)
    G = [x for x in tr if x[3] == D]
    S = sorted({v for k in E for v in Cd[k]})
    incr = []
    for k in E:
        c = sorted(Cd[k], reverse=True)
        incr += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    incr.sort(reverse=True)
    tau = incr[rp - 1] if len(incr) >= rp else None
    C = S[0] + S[-1]
    zero = measure([(x[2], x[1]) for x in tr], rp) is None
    top0 = (len(G) == 2
            and inv_of(G[0][1], rp) == inv_of(G[1][1], rp)
            and G[0][2] == -G[1][2])
    gcom = None
    if len(G) == 2:
        sa, sb = G[0][0], G[1][0]
        gcom = sorted(v for k, v in sa.items() if sb.get(k) == v)
    return dict(G=G, S=S, tau=tau, C=C, zero=zero, top0=top0, gcom=gcom, e=len(E))


def flag_top0(b, t, rp):
    """P7: [Phi]_top = 0 en todos los pisos de la bandera hasta rango 1."""
    cur, k = b, rp
    while k >= 1:
        a = anat(cur, t, k)
        if a is None:
            return True                   # una clase vacia: el piso se anula outright
        if not a["top0"]:
            return False
        nx = pelar(cur, t)
        if nx is None or len(nx) != t + 2 * (k - 1) or k == 1:
            return True
        cur, k = nx, k - 1
    return True


def preds(b, t, r, a):
    """los candidatos, todos legibles del voraz salvo P2/P7 que pelan."""
    tau, S = a["tau"], a["S"]
    pb = pelar(b, t)
    if pb is not None and len(pb) == t + 2 * (r - 1) and r >= 2:
        ap = anat(pb, t, r - 1)
        P2 = True if ap is None else ap["zero"]
    else:
        P2 = False
    P1 = a["top0"]
    P4 = a["gcom"] is not None and sorted(tau - v for v in a["gcom"]) == sorted(a["gcom"])
    P5 = sorted(tau - v for v in S) == S
    return {"P1 [Phi]_top=0": P1,
            "P2 pelado=0": P2,
            "P3 P1&P2": P1 and P2,
            "P4 g_com tau-sim": P4,
            "P5 S tau-sim": P5,
            "P6 P1&P4": P1 and P4,
            "P7 top0 en la bandera": flag_top0(b, t, r),
            "CTL e=2": a["e"] == 2}


def run():
    names = None
    print("=" * 100)
    print("POBLACION: ocupacion, |G| = 2 y C = tau  (la condicion necesaria del Corolario 8.32)")
    print("=" * 100)
    tot = {}
    for (t, r, W) in CFG:
        pop = nz = 0
        acc = {}
        for b in betas(t, r, W):
            a = anat(b, t, r)
            if a is None or a["tau"] is None:
                continue
            if len(a["G"]) != 2 or a["C"] != a["tau"]:
                continue
            pop += 1
            z = a["zero"]
            nz += z
            P = preds(b, t, r, a)
            if names is None:
                names = list(P)
            for k, v in P.items():
                d = acc.setdefault(k, [0, 0, 0, 0])       # TP, FP, FN, TN
                d[0 if (v and z) else 1 if (v and not z) else 2 if (not v and z) else 3] += 1
        print()
        print("  t=%d r=%d W=%d :  %d formas en poblacion, %d se anulan" % (t, r, W, pop, nz))
        print("     %-24s %6s %6s %6s %6s   %s" % ("predicado", "TP", "FP", "FN", "TN", "exacto"))
        for k in names:
            TP, FP, FN, TN = acc[k]
            print("     %-24s %6d %6d %6d %6d   %s"
                  % (k, TP, FP, FN, TN, "SI" if FP == 0 and FN == 0 else ""))
            g = tot.setdefault(k, [0, 0, 0, 0])
            for i in range(4):
                g[i] += acc[k][i]
    print()
    print("=" * 100)
    print("GLOBAL")
    print("=" * 100)
    print("     %-24s %6s %6s %6s %6s   %s" % ("predicado", "TP", "FP", "FN", "TN", "exacto"))
    exact = []
    for k in names:
        TP, FP, FN, TN = tot[k]
        ok = FP == 0 and FN == 0
        exact.append(ok)
        print("     %-24s %6d %6d %6d %6d   %s" % (k, TP, FP, FN, TN, "SI" if ok else ""))
    ctl = tot["CTL e=2"]
    if ctl[1] == 0 and ctl[2] == 0:
        print("\n  *** el CONTROL sale exacto: el arnes no esta midiendo nada")
        return 1
    print("\n  el control (e=2) falla, como debe: FP=%d FN=%d" % (ctl[1], ctl[2]))
    return 0 if any(exact[:-1]) else 0


if __name__ == "__main__":
    sys.exit(run())
