# -*- coding: utf-8 -*-
# ¿ESTA LA CANCELACION GENUINA CONTROLADA POR g_com?  El test decisivo -- 15 de agosto de 2026.
#
# POR QUE EL BARRIDO NO SERVIA.  cancel_vs_gcom.py pregunto lo mismo sobre los barridos exhaustivos y
# salio "SI, determinado" en las cuatro configuraciones.  Es una CONFIRMACION VACUA: en esos rangos
# TODOS los supervivientes tienen CANCEL = 1.  Sin variacion en la variable dependiente no hay test,
# y un empate no es un acuerdo (ver a-decoy-that-ties-means-untested).
#
# EL TEST QUE SI TIENE CONTENIDO.  Las familias de saturation.py separan los dos extremos de S en
# pasos de s, y a lo largo de la trayectoria CANCEL sube y luego se congela (2 2 2 ... 3 3 4 4 ... 4).
# Pero esa construccion mueve SOLO el maximo y el minimo de S: no toca ninguna otra clase.  Asi que
# si g_com se queda quieto mientras CANCEL sube,
#
#       CANCEL NO es funcion de |g_com|,
#
# y la conjetura de que la cancelacion genuina esta "controlada por la profundidad de reflexion"
# queda refutada en una sola trayectoria, sin necesidad de barrer nada.
#
# Se imprime la trayectoria entera -- W, prof, CANCEL, sin-soporte, e, |g_com|, |G| -- para que se vea
# QUE se mueve y que no.  Un negativo util dice tambien cual es la variable que si acompaña.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python cancel_traj_gcom.py

import itertools
import os

from second_stratum import setup

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]

CASOS = [(4, 2, 26, (18, 17, 11, 8, 7, 6, 1, 0), 14),
         (4, 2, 28, (18, 17, 11, 8, 7, 6, 1, 0), 14),
         (6, 3, 21, (16, 15, 14, 13, 11, 6, 5, 4, 3, 2, 1, 0), 12),
         (8, 3, 22, (20, 18, 15, 14, 13, 8, 7, 6, 5, 4, 3, 2, 1, 0), 10),
         (6, 4, 21, (20, 19, 17, 16, 15, 13, 8, 7, 6, 5, 4, 3, 1, 0), 10)]


def anatomia(beta, t, r):
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
    gc = set.intersection(*gs)
    return len(E), sorted(gc), len(G)


def extremos_S(beta, t):
    cl, E, Cd = setup(beta, t)
    S = sorted({v for k in E for v in Cd[k]})
    return S[-1], S[0]


print("=" * 100)
print("LA TRAYECTORIA: ¿se mueve g_com cuando CANCEL sube?")
print("=" * 100)
veredicto = []
for (t, r, s, seed, J) in CASOS:
    print("\n--- t=%d r=%d  paso s=%d  semilla %s" % (t, r, s, seed))
    print("      j |    W | prof | CANCEL | sin-sop |  e | |g_com| | g_com                | |G|")
    hi, lo = extremos_S(seed, t)
    cans, gcs = [], []
    for j in range(J + 1):
        b = tuple(sorted([(x + s * j if x == hi else (x - s * j if x == lo else x)) for x in seed],
                         reverse=True))
        rec = probe(b, t, r)
        if rec is None or not rec['surv'] or rec['prof'] is None:
            break
        a = anatomia(b, t, r)
        if a is None:
            break
        e, gc, nG = a
        cans.append(rec['vac_cancelan'])
        gcs.append(tuple(gc))
        print("    %3d | %4d | %4s | %6s | %7s | %2d | %7d | %-20s | %d"
              % (j, b[0] - b[-1], rec['prof'], rec['vac_cancelan'], rec['vac_sin_soporte'],
                 e, len(gc), str(gc)[:20], nG))
    if cans:
        movio_cancel = len(set(cans)) > 1
        movio_g = len(set(len(g) for g in gcs)) > 1
        veredicto.append((t, r, min(cans), max(cans), movio_cancel, movio_g))
        print("    -> CANCEL %d..%d (%s) ;  |g_com| %s (%s)"
              % (min(cans), max(cans), 'SE MUEVE' if movio_cancel else 'quieto',
                 sorted(set(len(g) for g in gcs)),
                 'SE MUEVE' if movio_g else 'QUIETO'))

print()
print("=" * 100)
for (t, r, c0, c1, mc, mg) in veredicto:
    if mc and not mg:
        print("  t=%d r=%d : CANCEL sube de %d a %d con |g_com| QUIETO"
              "  ->  CANCEL NO es funcion de |g_com|" % (t, r, c0, c1))
    elif mc and mg:
        print("  t=%d r=%d : los dos se mueven -- no concluye por si solo" % (t, r))
    else:
        print("  t=%d r=%d : CANCEL no varia en esta trayectoria -- no testa" % (t, r))
print("=" * 100)
