# -*- coding: utf-8 -*-
r"""peso_conforme.py -- what the conformal weight of a column sees: the second moment, and no more.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

(a) M_2(T_mu) = |mu+rho|^2 mod p on every column of SU(5)_12, SU(6)_11, SU(4)_19; and
    M_2(T_mu) = 0  <=>  p does not divide the order of exp(2 pi i |mu+rho|^2 / 2p)
    = theta_mu * exp(pi i |rho|^2 / p), with theta_mu = exp(2 pi i Delta_mu).
(b) genus one: for h = 1, delta = 1 <=> M_2 M_3 != 0, on every orbit of centred spectra with
    p <= 31, 3 <= n <= min(9, p-2); and the orbits with h >= 2 and delta = 1, counted.
(c) SU(6)_11, mu = (0,0,0,0,1) and nu = (0,2,0,0,5): same weight mod 1, same charge, trivial
    stabiliser, same field degree, different semigroups; indices also by integer lattices.
(d) SU(5)_14, mu = (5,0,4,2): M_2 = 1, M_3 = 0, S_T = <2,5>, index 19^2, also by lattices.
Semigroups: Gamma_T (equal to S_T for n <= 12, theorem of moments).  Public domain (CC0).
"""
import sys
import time
from fractions import Fraction

from P3a import (centrados, columna, espectro, estabilizador, gamma_T, momento, orden_columna,
                 pesos, semigrupo, factor_str)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
t0 = time.time()
MAL = []


def check(ok, msg):
    print(("  ok     " if ok else "  FAILED ") + msg)
    if not ok:
        MAL.append(msg)
    sys.stdout.flush()


print("(a) the congruence")
for n, k in ((5, 12), (6, 11), (4, 19)):
    p = n + k
    malos = malos_orden = 0
    cols = pesos(n, k)
    for mu in cols:
        s = espectro(n, k, mu)
        M2 = momento(s["T"], 2, p)
        if M2 != s["Q"] * pow(n, -1, p) % p:
            malos += 1
        orden = (Fraction(s["Q"], 2 * n * p)).denominator      # order of exp(2 pi i |mu+rho|^2 / 2p)
        if (M2 == 0) != (orden % p != 0):
            malos_orden += 1
    check(malos == 0 and malos_orden == 0,
          "SU(%d)_%d, all %d columns: M_2 = |mu+rho|^2 mod p; M_2 = 0 <=> p does not divide the order"
          % (n, k, len(cols)))

print("(b) genus one")
fallos, hmayor, tot = 0, [], 0
for p in (7, 11, 13, 17, 19, 23, 29, 31):
    for n in range(3, min(9, p - 2) + 1):
        for T in centrados(p, n):
            tot += 1
            h = len(estabilizador(T, p))
            d = semigrupo(gamma_T(T, p))["genero"]
            if h == 1 and (d == 1) != (momento(T, 2, p) != 0 and momento(T, 3, p) != 0):
                fallos += 1
            if h > 1 and d == 1:
                hmayor.append((p, n, T, h))
    print("     p = %d done (%d orbits so far)" % (p, tot))
    sys.stdout.flush()
check(fallos == 0, "h = 1: delta = 1 <=> M_2 M_3 != 0, all %d orbits with p <= 31, n <= 9" % tot)
check(len(hmayor) == 18, "orbits with h >= 2 and delta = 1: %d; first %s" % (len(hmayor), hmayor[:1]))

print("(c) same weight, charge, stabiliser and degree; different semigroups")
a, b = columna(6, 11, (0, 0, 0, 0, 1)), columna(6, 11, (0, 2, 0, 0, 5))
ra, rb = orden_columna(6, 11, (0, 0, 0, 0, 1)), orden_columna(6, 11, (0, 2, 0, 0, 5))
check(a["peso"] == b["peso"] == Fraction(35, 204) and a["carga"] == b["carga"] == 5
      and a["h"] == b["h"] == 1 and a["grado"] == b["grado"] == ra["grado"] == rb["grado"] == 32,
      "weight %s / %s, charge %d / %d, h %d / %d, degree %d / %d"
      % (a["peso"], b["peso"], a["carga"], b["carga"], a["h"], b["h"], a["grado"], b["grado"]))
check(a["S"]["minimales"] == [2, 3] and b["S"]["minimales"] == [2, 5]
      and (a["indice"], b["indice"], ra["OA"], rb["OA"]) == (17 ** 2, 17 ** 4, 17 ** 2, 17 ** 4),
      "semigroups <2,3> / <2,5>, indices %s / %s (semigroup) and %s / %s (lattices)"
      % (factor_str(a["indice"]), factor_str(b["indice"]), factor_str(ra["OA"]), factor_str(rb["OA"])))

print("(d) the second moment is not enough")
c = columna(5, 14, (5, 0, 4, 2))
r = orden_columna(5, 14, (5, 0, 4, 2))
check(c["momentos"][1] == 1 and c["momentos"][2] == 0 and c["S"]["minimales"] == [2, 5]
      and c["indice"] == r["OA"] == 19 ** 2,
      "SU(5)_14 (5,0,4,2): M_2 = %d, M_3 = %d, S_T = <2,5>, index %s (semigroup), %s (lattices)"
      % (c["momentos"][1], c["momentos"][2], factor_str(c["indice"]), factor_str(r["OA"])))
print()
print("FAILED: %d %s ; time %.0f s" % (len(MAL), MAL, time.time() - t0))
