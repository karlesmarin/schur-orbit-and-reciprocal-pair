# -*- coding: utf-8 -*-
# LA SATURACION CONTRA g_com -- 15 de agosto de 2026.
#
# DE DONDE VIENE.  LAW_RUNAWAY_EXTREMES.md midio, con familias explicitas hasta anchura ~276, que la
# profundidad se parte en dos sumandos con comportamientos OPUESTOS: los peldaños SIN SOPORTE crecen
# sin freno, y los que de verdad CANCELAN se SATURAN.  Lo que quedo abierto ese dia, escrito alli con
# todas las letras, es CUAL ES EL INVARIANTE: "lo que sobrevive -- y es lo que hay que probar -- es la
# saturacion, NO el numero", porque el valor depende de la semilla (4, 24, 18, 42 en las cuatro
# familias medidas).
#
# LA HIPOTESIS QUE SE MIDE AQUI.  Una consulta externa propone cerrar la Conjetura 8.44 pelando
# g_com por sus extremos -- probar c_1+c_{e-2}=C, luego c_2+c_{e-3}=C, etc. -- empujado por estratos
# soportados sucesivos, y conjetura que la cancelacion genuina esta "controlada por la profundidad de
# reflexion".  Si eso es asi, CANCEL tendria que ser una FUNCION de |g_com| a (t,r) fijos, porque
# |g_com|/2 es exactamente el numero de capas que quedan por reflejar.
#
# ESTO ES LO QUE SE MIDE, y no otra cosa:
#
#     a (t,r) fijos, ¿CANCEL queda determinado por |g_com|?
#
# Si a un mismo |g_com| le corresponden varios CANCEL, la respuesta es NO y el invariante es otro.
# Se imprime la tabla de contingencia entera, no un resumen: un "no" util dice ADEMAS cuanto se
# dispersa.
#
# POBLACION.  Los SUPERVIVIENTES: Phi != 0 con el estrato de arriba cancelado (CANCEL >= 1).  Son los
# unicos donde CANCEL es un numero finito interesante -- en los que se anulan cancela todo.
#
# CONTROL.  Se imprime tambien la tabla contra e = |E|, que es |g_com|+2 y por tanto la MISMA
# informacion.  Si CANCEL saliera determinado por |g_com| habria que comprobar que no es un artefacto
# de que |g_com| y (t,r) ya fijen la forma casi del todo; por eso se imprime cuantas betas distintas
# hay en cada celda.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python cancel_vs_gcom.py

import itertools
from collections import defaultdict

from second_stratum import setup, all_transversals
from depth_histogram import measure


def anatomia(beta, t, r):
    """(e, |g_com|, |G|) o None."""
    cl = {}
    for v in beta:
        cl.setdefault(v % t, []).append(v)
    if len(cl) < t:
        return None
    for i in cl:
        cl[i].sort(reverse=True)
    E = sorted(i for i in cl if len(cl[i]) >= 2)
    if not E:
        return None
    ks = [len(cl[i]) - 1 for i in E]
    best, G = None, []
    for sel in itertools.product(*[range(k + 1) for k in ks]):
        if sum(sel) != r:
            continue
        tot = sum(cl[E[a]][k] + cl[E[a]][k + 1]
                  for a, ki in enumerate(sel) for k in range(ki))
        if best is None or tot > best:
            best, G = tot, [sel]
        elif tot == best:
            G.append(sel)
    if not G:
        return None
    gs = [set(cl[E[a]][ki] for a, ki in enumerate(sel)) for sel in G]
    gcom = set.intersection(*gs)
    return len(E), len(gcom), len(G)


def betas(t, r, W):
    N = t + 2 * r
    for mid in itertools.combinations(range(1, W + 1), N - 1):
        yield tuple(sorted(mid, reverse=True)) + (0,)


CFG = [(4, 2, 16), (6, 2, 16), (8, 2, 17), (4, 3, 15), (6, 3, 15)]

print("=" * 96)
print("¿CANCEL queda determinado por |g_com| a (t,r) fijos?  Poblacion: supervivientes con CANCEL>=1")
print("=" * 96)
gran_det = gran_no = 0
for (t, r, W) in CFG:
    tab = defaultdict(lambda: defaultdict(int))
    nG_seen = defaultdict(set)
    n = 0
    for beta in betas(t, r, W):
        st = setup(beta, t)
        if st is None:
            continue
        cl, E, Cd = st
        if not E:
            continue
        tr = all_transversals(beta, cl, r, t)
        m = measure([(x[2], x[1]) for x in tr], r)
        if m is None:                       # se anula: cancela todo
            continue
        first, cancel, spec, allB = m
        if not cancel:                      # el estrato de arriba NO cancela
            continue
        a = anatomia(beta, t, r)
        if a is None:
            continue
        e, g, nG = a
        tab[g][len(cancel)] += 1
        nG_seen[g].add(nG)
        n += 1
    det = sum(1 for g in tab if len(tab[g]) == 1)
    nodet = sum(1 for g in tab if len(tab[g]) > 1)
    gran_det += det
    gran_no += nodet
    print("\n--- t=%d r=%d W=%d :  %d supervivientes con CANCEL>=1" % (t, r, W, n))
    print("    |g_com| | |G| |  CANCEL -> cuantas betas            | ¿determinado?")
    for g in sorted(tab):
        cells = '  '.join('%d:%d' % (c, tab[g][c]) for c in sorted(tab[g]))
        print("    %7d | %3s | %-34s | %s"
              % (g, ','.join(map(str, sorted(nG_seen[g]))), cells,
                 'SI' if len(tab[g]) == 1 else 'NO (%d valores)' % len(tab[g])))

print()
print("=" * 96)
print("  valores de |g_com| con CANCEL determinado: %d ;  con CANCEL disperso: %d"
      % (gran_det, gran_no))
print("=" * 96)
