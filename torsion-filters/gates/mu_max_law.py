# -*- coding: utf-8 -*-
# LA LEY DEL PESO SIMPLECTICO MAXIMO.  15 de agosto de 2026.
#
#     mu_max(beta)  =  top( Newt(N_beta) )  -  ( N-1, N-3, ..., N-2r+1 )
#
# donde N_beta es el NUMERADOR (el bialternante sin dividir), top(.) es el representante dominante
# maximo en dominancia de su poligono de Newton REAL -- el del soporte que sobrevive, no el de los
# candidatos --, y el vector que se resta es top(Newt(N_delta)), que solo depende de (t, r).
#
# DE DONDE SALE, Y NO ES UN AJUSTE.  Por Ostrowski, Newt(fg) = Newt(f) (+) Newt(g), luego
#
#     Newt(N_beta) = Newt(Phi) (+) Newt(N_delta)
#
# porque Phi = N_beta / N_delta.  Y para politopos de ORBITA de Weyl con u, v dominantes se tiene
# P(u) (+) P(v) = P(u+v) -- las funciones soporte se suman en la camara dominante y la simetria de
# Weyl lo extiende.  Asi que SI los tres politopos son politopos de orbita, la ley es forzada:
# restar los vertices dominantes.  Lo unico que hay que medir es que de verdad lo sean.
#
# EL CAMINO HASTA AQUI, porque las dos versiones ANTERIORES fallaron y por la misma razon:
#   - adivinando formulas cerradas contra invariantes de g_com: nada pasaba del 35 %;
#   - con  mu_max = v_max(beta) - v_max(delta),  donde v_max es el maximo de los CANDIDATOS v(T)
#     sobre las transversales: 39,2 %.  Falla porque el candidato maximo SE PUEDE CANCELAR -- que es
#     justo el fenomeno que ya habiamos medido en las 16 formas sin certificado.
#   El extremo de N no es el extremo de los candidatos: es el extremo de lo que SOBREVIVE.  Es la
#   tercera vez hoy que ese mismo descuido cambia un resultado.
#
# CONTROLES, y los tres pueden fallar
#   C0  la ley se comprueba contra mu_max calculado EN SAGE por el bialternante y la ramificacion
#       GL->Sp (sp_law.sage), que no comparte ni codigo ni libreria con esto.
#   C1  ¿es UNICO el vertice dominante del soporte real?  Si no lo fuera, Newt(N) no seria un
#       politopo de orbita y la ley no tendria por que valer.  Se mide en cinco configuraciones.
#   C2  SEÑUELO: la misma ley restando (N-1, N-2, ..., N-r), que es el otro vector "natural".
#       Tiene que FALLAR.  Si acertara igual, la ley no distingue nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python mu_max_law.py     (necesita sp_law_DUMP.json, producido por sp_law.sage en Sage)

import json
import os

from collision_graph import atoms, fibras
from peel_zero import betas, occupied, phi_zero


def domina(a, b):
    sa = sb = 0
    for x, y in zip(a, b):
        sa += x
        sb += y
        if sa < sb:
            return False
    return True


def top_real(beta, t, r):
    """el/los maximos en dominancia del soporte REAL de N_beta.  None si N_beta == 0."""
    at = atoms(beta, t, r)
    if at is None:
        return None
    f = fibras(at)
    sop = [m for m, l in f.items() if sum(c for c, _ in l) != 0]
    if not sop:
        return None
    dom = {tuple(sorted((abs(x) for x in m), reverse=True)) for m in sop}
    return [v for v in dom if not any(u != v and domina(u, v) for u in dom)]


def denom_top(t, r):
    """forma cerrada: (N-1, N-3, ..., N-2r+1).  Se comprueba, no se supone."""
    N = t + 2 * r
    return tuple(N - (2 * i + 1) for i in range(r))


print("=" * 104)
print("C1  ¿ES UNICO EL VERTICE DOMINANTE DEL SOPORTE REAL DE N?  -- lo que sostiene la ley")
print("=" * 104)
print("")
print("   t  r   W |  n | top real UNICO | top(N_delta) medido | forma cerrada (N-1,N-3,...) | igual")
for (t, r, W, tope) in [(4, 2, 14, 40), (6, 2, 15, 40), (8, 2, 16, 40), (4, 3, 13, 40), (6, 3, 15, 40)]:
    N = t + 2 * r
    DELTA = tuple(range(N - 1, -1, -1))
    td = top_real(DELTA, t, r)
    cerr = denom_top(t, r)
    n = uni = 0
    for b in betas(t, r, W):
        if not occupied(b, t):
            continue
        z = phi_zero(b, t, r)
        if z is None or z:
            continue
        tb = top_real(b, t, r)
        if tb is None:
            continue
        n += 1
        uni += (len(tb) == 1)
        if n >= tope:
            break
    if not n:
        print("  %3d %2d %3d | POBLACION VACIA -- esta fila no dice nada" % (t, r, W))
        continue
    print("  %3d %2d %3d | %2d | %5d (%5.1f%%) | %-19s | %-27s | %s"
          % (t, r, W, n, uni, 100.0 * uni / n, str(td), str(cerr),
             td == [cerr] if td else "--"))

# ===================================================================== C0 + C2 ===================
print("")
print("=" * 104)
print("C0  LA LEY, contra los mu_max calculados en SAGE   +   C2  el señuelo, que tiene que fallar")
print("=" * 104)
print("")
if not os.path.exists("sp_law_DUMP.json"):
    print("   FALTA sp_law_DUMP.json -- correr antes sp_law.sage en Sage.  Sin el no se dice nada.")
else:
    D = json.load(open("sp_law_DUMP.json"))
    t, r = 4, 2
    N = t + 2 * r
    vd = denom_top(t, r)
    senuelo = tuple(N - 1 - i for i in range(r))
    ok = mal = ok_s = 0
    fallos = []
    for d in D:
        tb = top_real(tuple(d['beta']), t, r)
        if tb is None or len(tb) != 1:
            mal += 1
            continue
        pred = tuple(a - b for a, b in zip(tb[0], vd))
        pre_s = tuple(a - b for a, b in zip(tb[0], senuelo))
        if list(pred) == d['mu']:
            ok += 1
        else:
            fallos.append((d, tb[0], pred))
        ok_s += (list(pre_s) == d['mu'])
    print("   poblacion: %d formas con g_com ASIMETRICO, mu_max calculado en Sage" % len(D))
    print("")
    print("   LEY      mu_max = top(N_beta) - %-12s : %3d de %d   (%5.1f %%)"
          % (str(vd), ok, len(D), 100.0 * ok / len(D)))
    print("   SEÑUELO  mu_max = top(N_beta) - %-12s : %3d de %d   (%5.1f %%)   <- debe fallar"
          % (str(senuelo), ok_s, len(D), 100.0 * ok_s / len(D)))
    print("   sin vertice unico: %d" % mal)
    if fallos:
        print("")
        print("   fallos de la ley (hasta 6):")
        for (d, tb, pred) in fallos[:6]:
            print("     %-27s top=%-9s pred=%-9s real=%s"
                  % (str(tuple(d['beta'])), str(tb), str(pred), str(tuple(d['mu']))))
    print("")
    if ok == len(D) and ok_s < len(D):
        print("   -> la ley PASA y el señuelo FALLA: el vector que hay que restar es el de la")
        print("      forma cerrada, y no otro.")

print("")
print("=" * 104)
print("  LO QUE LA LEY ES Y LO QUE NO ES.  Es exacta y esta forzada por Ostrowski, pero NO es una")
print("  formula cerrada en beta: para evaluarla hay que saber QUE vertice sobrevive en N_beta, que")
print("  es el mismo descenso de siempre.  Lo que hace es trasladarlo al nivel intrinseco y quitar")
print("  de en medio el denominador, que ahora es un vector explicito que solo depende de (t, r).")
print("=" * 104)
