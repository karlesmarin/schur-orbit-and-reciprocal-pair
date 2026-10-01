# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: barrido de caracteres CUADRATICOS impares chi = (-n/.), -n discriminante fundamental (n = 3 mod 4 libre de
#      cuadrados), p no divide n phi(n): vueltas del piso 2 (mahler_masas.cuadratico) frente a
#      min(lambda serie de Iwasawa, p-1) (filas_negativas.serie_iwasawa, rama <x>, N = 4: C(l(x), i) mod p con i < p^3).
#      Registra forzado = [chi(2) = 1] + [chi(p) = 1] para separar los casos no triviales. Rehace y GUARDA el barrido que el
#      auditor socratico hizo con un script borrado (2-oct-2026).
# FECHA: 2-oct-2026
# Uso (docker, /work): sage -python barrido_cuadratico_lambda.py p nmax
import sys
import time

sys.path.insert(0, '/work')
import filas_negativas as FN  # noqa: E402
import mahler_masas as MM  # noqa: E402
from sage.all import ZZ, euler_phi, kronecker_character  # noqa: E402

p, nmax = int(sys.argv[1]), int(sys.argv[2])
N = 4
Jk = [p - 1 + (k - 1) * (p - 1) for k in range(1, p)]
tot = ok = no_forz = no_forz_ok = 0
t0 = time.time()
for n in range(3, nmax + 1, 4):
    if not ZZ(-n).is_fundamental_discriminant() or n % p == 0 or euler_phi(n) % p == 0:
        continue
    chi = kronecker_character(-n)
    ex = [-1 if v == 0 else (0 if v == 1 else 1) for v in chi.values()]
    _, piv, c2, cp = MM.cuadratico(n, p, 2, Jk[-1])
    v = ''.join('P' if J in piv else 'A' for J in Jk)
    lt = v.find('P') if 'P' in v else p - 1
    mu, g, f = FN.medida(ex, 2, n, p, N)
    ls = FN.lam(FN.serie_iwasawa(mu, p, N, 0, p + 2))
    forz = int(c2 == 1) + int(cp == 1)
    bien = ls is not None and lt == min(ls, p - 1)
    tot += 1
    ok += bien
    if ls is not None and ls > forz:
        no_forz += 1
        no_forz_ok += bien
    if not bien or (ls or 0) > forz:
        print("n=%d: vueltas %s (lambda~ %d) | lambda serie %s | chi(2)=%s chi(p)=%s forzado %d%s" % (
            n, v, lt, ls, c2, cp, forz, '' if bien else '  <-- DIFIERE'), flush=True)
assert tot > 0
print("p=%d n<=%d: %d caracteres; vueltas == min(lambda, p-1): %d; por encima de lo forzado: %d (%d coinciden) (%.0f s)" % (
    p, nmax, tot, ok, no_forz, no_forz_ok, time.time() - t0), flush=True)
