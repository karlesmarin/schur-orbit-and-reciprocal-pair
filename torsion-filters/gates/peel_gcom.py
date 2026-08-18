# -*- coding: utf-8 -*-
# EL PELADO CONTRA g_com -- 15 de agosto de 2026.
#
# LA PREGUNTA, Y POR QUE NO ES LA OBVIA.  peel_zero.py midio que el locus es cerrado bajo pelado:
# Phi_{t,r}(beta) == 0  ==>  Phi_{t,r-1}(pelar beta) == 0, 705 ceros y 0 fallos.  La tentacion es
# preguntar si la Prop. 8.32 (los extremos de S evitan g_com) vale para la forma pelada.  ESO ES
# VACUO: la pelada se anula, y 8.32 es un TEOREMA para toda forma que se anula.  Preguntarlo seria
# medir un teorema.
#
# La pregunta con contenido es otra.  Lo que falta del criterio es la Conjetura 8.44,
#
#       Phi == 0  ==>  C - g_com = g_com,
#
# y lo que el paper YA prueba (eq:Ssymgen) es la reflexion sobre S \ g_com, o sea sobre todo MENOS
# la parte comun.  Si al pelar el centro se conserva, C' = C, entonces la forma pelada aporta la
# reflexion sobre OTRO trozo de S -- el suyo, S' \ g'_com -- respecto del MISMO centro.  Iterando,
# la union de esos trozos podria cubrir g_com, y entonces 8.44 sale por induccion en r.
#
# Asi que se mide, sobre las formas que SE ANULAN:
#
#   M1  C' = C ?                      -- si el centro no se conserva, la idea entera muere aqui
#   M2  |g'_com| < |g_com| ?          -- el pelado encoge la parte comun?
#   M3  g_com \subseteq (S' \ g'_com) union {los dos extremos quitados} ?
#                                     -- lo que hay que cubrir, cubierto en UN paso
#   M4  iterando hasta r=1, la union de los trozos con reflexion probada cubre S entera?
#
# CONTROL.  Sobre formas que NO se anulan pero tienen |G|=2 (o sea g_com existe), M1 tiene que
# fallar a menudo: si C'=C saliera siempre, no estaria midiendo la anulacion sino una identidad.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python peel_gcom.py

import itertools
from collections import Counter

from peel_zero import phi_zero


def clases(beta, t):
    cl = {}
    for v in beta:
        cl.setdefault(v % t, []).append(v)
    for i in cl:
        cl[i].sort(reverse=True)
    return cl


def anatomia(beta, t, r):
    """(S, C, g_com, |G|) o None si falta una clase o no hay exceso."""
    cl = clases(beta, t)
    if len(cl) < t:
        return None
    E = sorted(i for i in cl if len(cl[i]) >= 2)
    if not E:
        return None
    S = sorted(v for i in E for v in cl[i])
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
    gs = [set(cl[E[a]][ki] for a, ki in enumerate(sel)) for sel in G]
    gcom = set.intersection(*gs) if gs else set()
    return S, min(S) + max(S), gcom, len(G)


def pelar(beta, t):
    a = anatomia(beta, t, 1)          # el rango no afecta a S
    if a is None:
        return None
    S = a[0]
    hi, lo = max(S), min(S)
    return tuple(x for x in beta if x != hi and x != lo)


def betas(t, r, W):
    N = t + 2 * r
    for mid in itertools.combinations(range(1, W + 1), N - 1):
        yield tuple(sorted(mid, reverse=True)) + (0,)


CFG = [(4, 2, 15), (6, 2, 15), (8, 2, 16), (4, 3, 14)]

print("=" * 100)
print("EL PELADO CONTRA g_com -- lo que aporta cada capa, medido sobre las formas que SE ANULAN")
print("=" * 100)
print("  t  r  W  |  ceros | M1 C'=C | M2 |g'|<|g| | M3 g_com cubierto | M4 iterando cubre S |"
      " control C'=C en no-nulos")
tot_ceros = tot_m1 = tot_m3 = tot_m4 = 0
for (t, r, W) in CFG:
    ceros = m1 = m2 = m3 = m4 = 0
    ctrl_tot = ctrl_hit = 0
    for beta in betas(t, r, W):
        a = anatomia(beta, t, r)
        if a is None:
            continue
        S, C, gcom, nG = a
        cero = phi_zero(beta, t, r) is True
        if not cero:
            if nG == 2:
                ctrl_tot += 1
                b2 = pelar(beta, t)
                a2 = anatomia(b2, t, r - 1) if b2 else None
                if a2 and a2[1] == C:
                    ctrl_hit += 1
            continue
        ceros += 1
        b2 = pelar(beta, t)
        a2 = anatomia(b2, t, r - 1) if b2 else None
        if a2 is None:
            continue
        S2, C2, gcom2, nG2 = a2
        if C2 == C:
            m1 += 1
        if len(gcom2) < len(gcom):
            m2 += 1
        quitados = {max(S), min(S)}
        if gcom <= (set(S2) - gcom2) | quitados:
            m3 += 1
        # M4: iterar hasta agotar el rango, acumulando los trozos con reflexion probada
        cubierto = (set(S) - gcom) | quitados
        b, rr = b2, r - 1
        while rr >= 1:
            aa = anatomia(b, t, rr)
            if aa is None:
                break
            Sx, Cx, gx, _ = aa
            if Cx != C:
                break
            cubierto |= (set(Sx) - gx) | {max(Sx), min(Sx)}
            nb = pelar(b, t)
            if nb is None or rr == 1:
                break
            b, rr = nb, rr - 1
        if set(S) <= cubierto:
            m4 += 1
    tot_ceros += ceros
    tot_m1 += m1
    tot_m3 += m3
    tot_m4 += m4
    pct = (100.0 * ctrl_hit / ctrl_tot) if ctrl_tot else 0.0
    print("  %2d %2d %2d  | %6d | %7d | %11d | %17d | %19d | %5d/%-5d = %4.1f%%"
          % (t, r, W, ceros, m1, m2, m3, m4, ctrl_hit, ctrl_tot, pct))

print()
print("=" * 100)
print("  ceros: %d   M1 C'=C: %d   M3 g_com cubierto en un paso: %d   M4 iterando cubre S: %d"
      % (tot_ceros, tot_m1, tot_m3, tot_m4))
print("=" * 100)
