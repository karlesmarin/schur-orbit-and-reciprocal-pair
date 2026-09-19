#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""medir_atajo.py -- cuanto ahorra el puente, medido y no afirmado.

El atajo calcula el indice del orden NO REAL a partir del real:
        v_p([O:B]) = 2 v_p([O_+ : R]) + P d.
El camino lento lo calcula directamente, en grado phi(q); el atajo, en grado phi(q)/2.  Aqui se
cronometran los dos en las mismas celdas y se comprueba ademas que COINCIDEN, que es lo que da
derecho a usar el rapido.

Uso:  python medir_atajo.py

Authors: Carles Marin, Claude (AI assistant).  Public domain (CC0).
"""
import time

from segmento import euler_phi, exps_GL, indice, puente, valuacion

CELDAS = [(5, 45, 3), (7, 63, 3), (9, 63, 7), (11, 99, 3), (13, 117, 3),
          (7, 35, 5), (11, 55, 5), (15, 45, 3), (23, 207, 3)]

print("%-5s %-6s %-4s %-6s %-9s %-9s %-11s %-11s %-8s %s" % (
    "n", "q", "p", "phi(q)", "slow", "short", "t slow (s)", "t short (s)", "x", "agree"))
tl = ta = 0.0
ok = True
for (n, q, p) in CELDAS:
    r = q // p ** valuacion(q, p)
    t0 = time.perf_counter()
    lento = indice(exps_GL(n), q, p, dim=euler_phi(r))[0]
    t1 = time.perf_counter()
    rapido, P, d, hip = puente(n, q, p)
    t2 = time.perf_counter()
    dl, da = t1 - t0, t2 - t1
    tl += dl
    ta += da
    ig = (lento == rapido)
    ok = ok and ig
    print("%-5d %-6d %-4d %-6d %-9d %-9s %-11.3f %-11.3f %-8s %s" % (
        n, q, p, euler_phi(q), lento, rapido, dl, da,
        ("%.1f" % (dl / da)) if da > 0 else "-", "yes" if ig else "*** NO ***"))

print()
print("total: slow %.3f s, shortcut %.3f s, factor %.1f" % (tl, ta, tl / ta if ta else 0))
print("the two paths give the same integer in every cell: %s" % ("yes" if ok else "NO"))
print()
print("What is measured is how the work is split, not an optimal implementation: both routes use")
print("the same layer code, and the only difference is the DEGREE at which they close -- phi(q) on")
print("the slow path and phi(q)/2 on the shortcut.  The factor grows with phi(q).")
