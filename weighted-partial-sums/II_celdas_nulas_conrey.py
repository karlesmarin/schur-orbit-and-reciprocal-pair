# -*- coding: utf-8 -*-
r"""II_celdas_nulas_conrey.py -- intervalos donde la serie seno de Conrey se anula identicamente.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

f_q(x) = sum_n (n/q) sin(2 pi n x)/n^2 (q = 3 mod 4).  Por la identidad f_q(x) = (2 pi^2/q^{3/2}) F_q(qx),
f_q es identicamente 0 en [r/q, (r+1)/q] exactamente cuando r esta en R_q.
Se comprueba con la SERIE (truncada a NT terminos, error de cola < 1/NT) en puntos interiores:
  q = 103: celdas r = 47, 51, 55  (x en [47/103, 48/103], etc.) -> |f_q| < 1/NT;
  q = 7, 23, 31, 47: la celda central r = (q-1)/2 -> |f_q| < 1/NT;
  control: q = 103, celda r = 46 (no esta en R_q) -> |f_q| claramente mayor que 1/NT;
  control: q = 59, 83 (q = 3 mod 8), celda central -> no nula.
"""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
NT = 2_000_000


def leg(q):
    t = np.full(q, -1.0); t[0] = 0.0
    r = np.arange(1, q)
    t[(r * r) % q] = 1.0
    return t


n = np.arange(1, NT + 1, dtype=np.float64)


def fq(q, x):
    chi = leg(q)[np.arange(1, NT + 1) % q]
    return float(np.sum(chi * np.sin(2 * np.pi * n * x) / n ** 2))


def celda(q, r, etiqueta):
    vals = [fq(q, (r + t) / q) for t in (0.25, 0.5, 0.75)]
    print("q=%3d celda r=%3d %-28s f_q en puntos interiores: %s   (1/NT = %.1e)"
          % (q, r, etiqueta, ["%.2e" % v for v in vals], 1 / NT))


for r in (47, 51, 55):
    celda(103, r, "(en R_103)")
celda(103, 46, "(control: no en R_103)")
for q in (7, 23, 31, 47):
    celda(q, (q - 1) // 2, "(central, q = 7 mod 8)")
for q in (59, 83):
    celda(q, (q - 1) // 2, "(control: central, q = 3 mod 8)")
