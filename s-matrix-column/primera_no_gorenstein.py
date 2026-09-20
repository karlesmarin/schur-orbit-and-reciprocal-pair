# -*- coding: utf-8 -*-
r"""primera_no_gorenstein.py -- where the Gorenstein property first fails, at prime height.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

(1) The semigroups S_T (exact formula) that occur for n = 3, 4 and every prime p <= 31: only N_0,
    <2,3>, <3,4>, all symmetric.
(2) The first non-symmetric S_T in the order of n and then of the prime height p (n <= 6, p <= 31).
(3) The column that realises it, SU(5)_12, mu = (0,0,2,1): spectrum, moments, the central order m,
    and, by integer lattices (P3a.orden_columna), [O:A], the conductor index [A:f] and [O:f]:
    Gorenstein <=> [A:f] = [O:A].  Control: mu = 0, maximal.
Public domain (CC0).
"""
import sys

from P3a import (centrados, columna, es_primo, orden_columna, s_exacto, semigrupo, factor_str,
                 estabilizador, momento)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
PR = [p for p in range(5, 32) if es_primo(p)]

print("(1) n = 3, 4, p <= 31: semigroups that occur (by their gaps)")
for n in (3, 4):
    vistos = set()
    for p in PR:
        if n <= p - 2:
            for T in centrados(p, n):
                S = semigrupo(s_exacto(T, p))
                vistos.add((tuple(S["huecos"]), tuple(S["minimales"]), S["simetrico"]))
    print("   n=%d:" % n, sorted(vistos))

print("(2) the first non-symmetric S_T, by n and then by p")
primero = None
for n in range(3, 7):
    for p in PR:
        if n > p - 2 or primero:
            continue
        for T in centrados(p, n):
            S = semigrupo(s_exacto(T, p))
            if not S["simetrico"]:
                primero = (n, p, T, S)
                break
n, p, T, S = primero
print("   n=%d, p=%d, orbit of T=%s: S_T = <%s>, gaps %s, h = %d, M_1..M_5 = %s"
      % (n, p, T, ",".join(map(str, S["minimales"])), S["huecos"], len(estabilizador(T, p)),
         [momento(T, r, p) for r in range(1, 6)]))

print("(3) the column SU(5)_12, mu = (0,0,2,1), and the control mu = 0")
for mu in ((0, 0, 2, 1), (0, 0, 0, 0)):
    c = columna(5, 12, mu)
    r = orden_columna(5, 12, mu)
    print("   mu=%s: T=%s, L=%d (L mod 5 = %d), m=%d, S_T gaps %s, conductor %d; semigroup index %s"
          % (mu, c["T"], c["L"], c["L"] % 5, c["m"], c["S"]["huecos"], c["S"]["conductor"],
             factor_str(c["indice"])))
    print("        lattices: degree %d, [O:A] = %s, [A:f] = %s, [O:f] = %s -> %s"
          % (r["grado"], factor_str(r["OA"]), factor_str(r["Af"]), factor_str(r["Of"]),
             "Gorenstein" if r["gorenstein"] else "NOT Gorenstein"))
