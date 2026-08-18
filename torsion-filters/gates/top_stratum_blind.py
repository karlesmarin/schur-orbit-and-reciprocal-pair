# -*- coding: utf-8 -*-
"""EL ESTRATO DE ARRIBA NO PUEDE DECIDIR g_com -- y aqui esta el testigo que lo prueba.

EL PROBLEMA.  Falta una sola ecuacion para cerrar el criterio: tau - g_j = g_{j*}, o sea que g_com
sea tau-simetrico.  El Corolario 8.32 da todo lo demas y no da eso.  La pregunta era si se podia
sacar con mas trabajo del mismo estrato.

LA RESPUESTA ES NO, Y ES UNA PRUEBA, no una impresion.  Este guion agrupa las formas por
(tau, S\g_com, clases empatadas) -- que es TODO lo que el estrato de grado maximo ve, porque
P(T_a) y P(T_b) se calculan sobre T_a, T_b = S \ g_a, S \ g_b, y esos conjuntos omiten g_com por
definicion -- y busca grupos que contengan a la vez una forma que se anula y otra que no.

EXISTEN.  Ejemplo minimo, t=4 r=2, tau=12, empatadas {0,2}, S\g_com = [0,1,2,10,11,12]:

    beta = (12,11,10, 9,3,2,1,0)   g_com = [3,9]   3+9 = 12 = tau   SE ANULA
    beta = (12,11,10, 7,5,2,1,0)   g_com = [5,7]   5+7 = 12 = tau   SE ANULA
    beta = (12,11,10, 5,3,2,1,0)   g_com = [3,5]   3+5 =  8         no
    beta = (12,11,10, 9,7,2,1,0)   g_com = [7,9]   7+9 = 16         no

Cuatro formas con EL MISMO estrato de arriba y dos comportamientos.  Luego ningun argumento que use
solo [det]_{D1} = 0 puede cerrar la ecuacion: la informacion no esta ahi.  Hay que bajar al primer
estrato en el que un elemento de g_com entra en T.

Y de paso el grupo es una confirmacion limpia del criterio en una familia controlada: dentro de cada
grupo se anulan EXACTAMENTE las de g_com tau-simetrico.  En el testigo de nueve miembros, las tres
simetricas y ninguna de las seis restantes.

Authors: Carles Marin, Claude (AI assistant).
"""
import sys, itertools
sys.path.insert(0, r"E:\proyectos\Curiosity\research\orbit-pair\gates")
from fractions import Fraction as F
from criterion_control import betas, lam_of, value

PTS = [[F(3, 2), F(5, 3), F(7, 4)], [F(5, 2), F(7, 3), F(9, 4)],
       [F(4, 3), F(9, 5), F(11, 6)], [F(7, 5), F(11, 4), F(13, 7)]]


def is_zero(b, t, r):
    lam = lam_of(list(b))
    for p in PTS:
        if value(lam, t, r, p[:r]) != 0:
            return False
    return True


def greedy(b, t, r):
    cl = {}
    for v in b:
        cl.setdefault(v % t, []).append(v)
    if len(cl) < t:
        return None
    E = {k: sorted(v, reverse=True) for k, v in cl.items() if len(v) >= 2}
    if not E:
        return None
    inc = []
    for k, c in E.items():
        inc += [(c[i] + c[i + 1], k) for i in range(len(c) - 1)]
    inc.sort(key=lambda x: -x[0])
    if len(inc) < r + 1:
        return None
    tau = inc[r - 1][0]
    if inc[r][0] != tau:
        return None                       # |G| = 1
    kj = {k: sum(1 for (d, kk) in inc if kk == k and d > tau) for k in E}
    tied = sorted({k for (d, k) in inc if d == tau})
    if len(tied) != 2:
        return None
    gcom = sorted(E[k][kj[k]] for k in E if k not in tied)
    S = sorted(x for c in E.values() for x in c)
    if S[0] + S[-1] != tau:
        return None                       # C != tau
    rest = sorted(set(S) - set(gcom))
    if sorted(tau - x for x in rest) != rest:
        return None                       # S \ g_com no es tau-simetrico
    return tau, tuple(rest), tuple(gcom), tuple(tied)


CFG = [(4, 2, 17), (6, 2, 17), (8, 2, 17), (4, 3, 15), (6, 3, 15)]
print("=" * 100)
print("PREGUNTA: dos formas con el MISMO tau, el MISMO S\\g_com y las MISMAS clases empatadas,")
print("          una que se anula y otra que no.  Si existe, el estrato de arriba NO puede decidir")
print("          y la ecuacion que falta tiene que venir de un estrato mas bajo.")
print("=" * 100)
found = 0
for (t, r, W) in CFG:
    groups = {}
    for b in betas(t, r, W):
        G = greedy(b, t, r)
        if G is None:
            continue
        tau, rest, gcom, tied = G
        groups.setdefault((tau, rest, tied), []).append((b, gcom))
    hits = 0
    for key, lst in groups.items():
        if len(lst) < 2:
            continue
        zs = [(b, g, is_zero(b, t, r)) for (b, g) in lst]
        if any(z for _, _, z in zs) and any(not z for _, _, z in zs):
            hits += 1
            if found < 3:
                found += 1
                tau, rest, tied = key
                print("\n  TESTIGO  t=%d r=%d   tau=%d  tied=%s" % (t, r, tau, list(tied)))
                print("           S\\g_com = %s   (tau-simetrico)" % list(rest))
                for b, g, z in zs:
                    print("           beta=%-34s g_com=%-14s  %s"
                          % (list(b), list(g), "SE ANULA" if z else "no se anula"))
    print("  t=%d r=%d: %d grupos con las dos cosas" % (t, r, hits))
print()
print("  Si hay testigos: probado que [det]_{D1} = 0 no determina g_com.")
