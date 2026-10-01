# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: contraste de la formula cerrada para TODAS las filas y componentes (sin tope):
#        para 0 < J < p^a, con b = J mod (p-1) y u(J) = (J+1)/p^{v_p(J+1)}:   J ausente  <=>  u(J) < p * r_b,
#      r_b = ord_T f_b, f_b(T) = int_{Z_p^x} omega(u)^b (1+T)^{l(u)} dmu (rama omega^b; infinito si f_b = 0 mod P);
#      la fila 0 falta sii h0 = mu(Z_p) = 0 mod P.
#      Medido: pivotes de mahler_masas.matriz_B en todas las filas J < p^a. r_b: filas_negativas.serie_iwasawa.
#      Casos: (i) orbitas de la tabla de capas (n <= 50, p = 5, 7; pisos 2 y 3); (ii) cuadraticos con lambda de la rama
#      trivial >= p (tras el tope), p = 5, piso 4 (625 columnas); (iii) consecuencia k* = ceil(p r / (p-1)).
# FECHA: 2-oct-2026
# Uso (docker, /work, -e OPENBLAS_NUM_THREADS=4): sage -python formula_todas_filas.py
import math
import sys
import time

import numpy as np

sys.path.insert(0, '/work')
import filas_negativas as FN  # noqa: E402
import mahler_masas as MM  # noqa: E402
import mahler_padico as MP  # noqa: E402
from sage.all import ZZ, euler_phi, kronecker_character  # noqa: E402


def u(J, p):
    x = J + 1
    while x % p == 0:
        x //= p
    return x


def prediccion(p, a, r, h0nz):
    aus = set()
    if not h0nz:
        aus.add(0)
    for J in range(1, p ** a):
        rb = r[J % (p - 1)]
        if rb is None or u(J, p) < p * rb:
            aus.add(J)
    return aus


def ramas(mu, p, N, imax):
    return [FN.lam(FN.serie_iwasawa(mu, p, N, b, imax)) for b in range(p - 1)]


def h0_no_nulo(mu, p):
    return bool(((mu.sum(axis=0) % p) != 0).any())


tot = ok = 0
t0 = time.time()
for p, pisos in ((5, (2, 3)), (7, (2, 3))):
    for n in range(3, 51, 2):
        if n % p == 0 or euler_phi(n) % p == 0:
            continue
        obs, _ = MP.orbitas(n, p)
        for chi, orb, m, Km, chim in obs:
            ex = MP.exponentes(chim, Km, n)
            mu, g, f = FN.medida(ex, m, n, p, 4)
            r = ramas(mu, p, 4, 3 * p)
            h0 = h0_no_nulo(mu, p)
            for a in pisos:
                P = p ** a
                piv = set(MM.matriz_B(ex, ex[2], m, n, p, a, P - 1).pivot_rows())
                med = {J for J in range(P) if J not in piv}
                pre = prediccion(p, a, r, h0)
                tot += 1
                ok += (med == pre)
                if med != pre:
                    print("  DIFIERE n=%d p=%d a=%d r=%s: sobran %s faltan %s" % (
                        n, p, a, r, sorted(pre - med)[:8], sorted(med - pre)[:8]), flush=True)
    print("p=%d: %d de %d (orbita, piso) con TODAS las filas (%.0f s)" % (p, ok, tot, time.time() - t0), flush=True)
# (ii) tras el tope: cuadraticos p = 5 con r_0 >= 5, piso 4
p, a = 5, 4
tope = tope_ok = kst_ok = 0
for n in [75619] + list(range(3, 60000, 4)):
    if not ZZ(-n).is_fundamental_discriminant() or n % p == 0 or euler_phi(n) % p == 0:
        continue
    chi = kronecker_character(-n)
    ex = [-1 if v == 0 else (0 if v == 1 else 1) for v in chi.values()]
    mu, g, f = FN.medida(ex, 2, n, p, 4)
    r0 = FN.lam(FN.serie_iwasawa(mu, p, 4, 0, 3 * p))
    if r0 is None or r0 < p:
        continue
    r = ramas(mu, p, 4, 3 * p)
    h0 = h0_no_nulo(mu, p)
    P = p ** a
    piv = set(MM.cuadratico(n, p, a, P - 1)[1])
    med = {J for J in range(P) if J not in piv}
    pre = prediccion(p, a, r, h0)
    tope += 1
    tope_ok += (med == pre)
    # k*: primera vuelta presente del piso 3, leida en el piso 4 (J < 625)
    Jk = [p * p - 1 + (k - 1) * p * (p - 1) for k in range(1, 26) if p * p - 1 + (k - 1) * p * (p - 1) < P]
    kst = next((k for k, J in enumerate(Jk, 1) if J in piv), None)
    kst_ok += (kst == math.ceil(p * r0 / (p - 1)))
    if med != pre or kst != math.ceil(p * r0 / (p - 1)):
        print("  tras el tope n=%d r=%s: formula %s, k* medido %s predicho %d" % (
            n, r, 'IGUAL' if med == pre else 'DIFIERE', kst, math.ceil(p * r0 / (p - 1))), flush=True)
    if tope >= 40:
        break
assert tot > 0 and tope > 0
print("FORMULA u(J) < p r_b: tabla %d de %d; tras el tope (piso 4) %d de %d; k* = ceil(p r/(p-1)): %d de %d (%.0f s)" % (
    ok, tot, tope_ok, tope, kst_ok, tope, time.time() - t0), flush=True)
