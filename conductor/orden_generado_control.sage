# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# CONTROL DE og_orden (orden_generado.sage) CONTRA L.order DE SAGE.
#
#   C1  en casillas BARATAS (prod de grados pequeno), og_orden(L, simetricas) y L.order(simetricas)
#       tienen que dar el MISMO reticulo (HNF igual) y el mismo indice en Z[a].
#   C2  SENUELO: og_orden con UNA simetrica quitada tiene que dar un reticulo DISTINTO en al menos una
#       casilla donde Sage tambien lo da distinto.  Si og_orden devolviera siempre Z[a] (o siempre lo
#       mismo), C1 pasaria en las casillas donde A = Z[a] y no habria medido nada.  Se cuentan aparte
#       las casillas de C1 con indice > 1: esas son las que miden.
#   C3  la casilla que motivo esto, (7,19), con og_orden: tiene que terminar, y se imprime el tiempo.
#
# REPL-SAFE.
#
# Authors: Carles Marin, Claude (AI assistant).
import os
import sys
import time
load("orden_generado.sage")
x = polygen(QQ, 'x')
def psi(q):
    K = CyclotomicField(q)
    z = K.gen()
    return (z + z**(-1)).minpoly(x)
def etas_de(m, L):
    a = L.gen()
    out = [a]
    if m >= 2:
        out.append(a*a - 2)
    for j in range(2, m):
        out.append(a*out[j-1] - out[j-2])
    return out[:m]
def simetricas(m, L):
    R = PolynomialRing(L, 'T')
    T = R.gen()
    P = prod([T - w for w in etas_de(m, L)])
    return [(-1)**i * P[m-i] for i in range(1, m+1)]
def hnf_orden(L, A):
    return og_hnf(L, [L(b) for b in A.basis()])
BARATAS = [(2,20),(3,21),(3,28),(4,20),(4,27),(4,36),(5,33),(6,21),(6,28),(7,21),(7,28),(2,13),(3,15),(4,15),(3,16),(5,16),(4,24),(5,24),(6,13),(6,15),(6,16),(6,17)]
c1med = 0
c1mal = 0
c1idx = 0
c2med = 0
c2dist = 0
print("   m    q   d  prod(grados)   idx_sage  idx_og   HNF igual   decoy: sage distinto / og distinto")
for (m, q) in BARATAS:
    L = NumberField(psi(q), 'a')
    d = L.degree()
    a = L.gen()
    OL = L.order(a)
    ss = simetricas(m, L)
    pg = prod([s.absolute_minpoly().degree() for s in ss])
    if pg > 300000:
        print("  %3d %4d %3d %12d   saltada: cara para Sage" % (m, q, d, pg))
        continue
    # rango incompleto (p.ej. m = numero de clases: todas las simetricas son racionales) es un
    # resultado a COMPARAR, no a descartar: los dos constructores tienen que fallar juntos.
    try:
        As = L.order(ss)
    except ValueError:
        As = None
    try:
        Ao = og_orden(L, ss)
    except ValueError:
        Ao = None
    if As is None or Ao is None:
        c1med += 1
        juntos = (As is None) and (Ao is None)
        if not juntos:
            c1mal += 1
        print("  %3d %4d %3d %12d   rango incompleto: sage %s, og %s   %s" % (m, q, d, pg, "falla" if As is None else "ok", "falla" if Ao is None else "ok", "SI" if juntos else "*** NO ***"))
        continue
    Hs = hnf_orden(L, As)
    Ho = hnf_orden(L, Ao)
    iS = ZZ(As.index_in(OL))
    iO = ZZ(Ao.index_in(OL))
    igual = (Hs == Ho) and (iS == iO)
    c1med += 1
    if not igual:
        c1mal += 1
    if iS > 1:
        c1idx += 1
    ds = "-"
    if len(ss) >= 2:
        sub = ss[:-1]
        try:
            Hs2 = hnf_orden(L, L.order(sub))
        except ValueError:
            Hs2 = "rango"
        try:
            Ho2 = hnf_orden(L, og_orden(L, sub))
        except ValueError:
            Ho2 = "rango"
        if True:
            sd = isinstance(Hs2, str) or (Hs2 != Hs)
            od = isinstance(Ho2, str) or (Ho2 != Ho)
            if isinstance(Hs2, str) != isinstance(Ho2, str):
                c1mal += 1
                print("     *** rango: sage dice %s, og dice %s en (%d,%d) sin la ultima ***" % (Hs2 if isinstance(Hs2, str) else "ok", Ho2 if isinstance(Ho2, str) else "ok", m, q))
            elif not isinstance(Hs2, str) and ((Hs2 == Ho2) != True):
                c1mal += 1
                print("     *** el senuelo da reticulos distintos en sage y og en (%d,%d) ***" % (m, q))
            c2med += 1
            if sd and od:
                c2dist += 1
            ds = "%s / %s" % (sd, od)
    print("  %3d %4d %3d %12d   %8s %7s   %-9s   %s" % (m, q, d, pg, iS, iO, "SI" if igual else "*** NO ***", ds))
    sys.stdout.flush()
print("")
print("   C1  og_orden == L.order :  medidas %d  (con indice > 1: %d)   discrepancias %d   %s" % (c1med, c1idx, c1mal, ("PASA" if c1mal == 0 else "*** NO PASA ***") if c1idx > 0 else "*** sin casillas de indice > 1: no mide ***"))
print("   C2  senuelo (quitar una simetrica cambia el anillo, y og lo ve):  medidas %d   ambos distintos %d   %s" % (c2med, c2dist, "PASA" if c2dist > 0 else "*** el senuelo nunca cambia nada: C1 no distingue ***"))
t0 = time.time()
L = NumberField(psi(19), 'a')
A = og_orden(L, simetricas(7, L))
print("   C3  (7,19) con og_orden:  indice %s   %.1f s" % (A.index_in(L.order(L.gen())), time.time() - t0))
