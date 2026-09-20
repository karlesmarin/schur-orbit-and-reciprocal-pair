# -*- coding: utf-8 -*-
r"""exponente.py -- the exponent of O_E/C_T: is p O_E contained in C_T?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Corollary 'the exponent': if c(S_T) <= e then p O_mu is contained in A_mu.  Two independent readings
per orbit of centred spectra:
  (i)  c(S_T) <= e, from the exact semigroup;
  (ii) p * (every basis vector of O_E) lies in the lattice C_T, solved on its Hermite basis
       (P3a.exponente_p), without the semigroup.
Control that must fail where it should: O_E itself lies in C_T exactly when delta = 0.
Reports the largest c/e.

  python exponente.py              the two ranges of the note: 3 <= n <= 6, p <= 31  and
                                   3 <= n <= 8, p <= 23
Public domain (CC0).
"""
import sys
from fractions import Fraction

from P3a import centrados, es_primo, exponente_p, delta_red, s_exacto, semigrupo, estabilizador

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def rango(nmax, pmax):
    tot = mal_i = mal_ii = ctrl = 0
    peor = (Fraction(0), None)
    for p in range(5, pmax + 1):
        if not es_primo(p):
            continue
        for n in range(3, min(nmax, p - 2) + 1):
            for T in centrados(p, n):
                tot += 1
                e = (p - 1) // len(estabilizador(T, p))
                S = semigrupo(s_exacto(T, p))
                c = S["conductor"]
                if c > e:
                    mal_i += 1
                    print("   c(S_T) > e:", p, n, T, c, e)
                pO, O = exponente_p(T, p)
                if not pO:
                    mal_ii += 1
                    print("   p O_E not in C_T:", p, n, T)
                if O != (delta_red(T, p) == 0):
                    ctrl += 1
                if Fraction(c, e) > peor[0]:
                    peor = (Fraction(c, e), (p, n, T, e, c))
    print("3 <= n <= %d, p <= %d: %d orbits | (i) c(S_T) > e: %d | (ii) p O_E not in C_T: %d | "
          "control (O_E in C_T <=> delta = 0) broken: %d" % (nmax, pmax, tot, mal_i, mal_ii, ctrl))
    print("   largest c/e = %s at (p, n, T, e, c) = %s" % (peor[0], peor[1]))
    sys.stdout.flush()
    return peor[0]


a = rango(6, 31)
b = rango(8, 23)
print("largest c/e over both ranges: %s  (<= 5/16: %s)" % (max(a, b), max(a, b) <= Fraction(5, 16)))
