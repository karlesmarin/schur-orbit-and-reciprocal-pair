# -*- coding: utf-8 -*-
r"""censo.py -- the census of centred spectra: delta(T) as a lattice index against g(Gamma_T), g(S_T).

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

For every prime PMIN <= p <= PMAX and NMIN <= n <= min(NMAX, p-2), one spectrum T per orbit of centred
n-subsets of F_p under F_p^x:
  delta = v_p[O_E : C_T], computed as the determinant of an integer lattice (P3a.delta_red,
          Hermite form over Z; it does not use the semigroup);
  g(Gamma_T), Gamma_T = < e, j : M_hj(T) != 0 >  (e always among the generators);
  g(S_T),     S_T = < e, j + e b_j >  with b_j from Teichmuller lifts (P3a.s_exacto).
Prints the number of orbits per (p, n), every disagreement, the largest delta per n, and whether
p <= (n-1)(n-2) (below the bound of the theorem of moments).

  python censo.py                  n = 3..6, 5 <= p <= 31            (note: 2414 orbits)
  python censo.py bajo             n = 7, 8, 9 with p <= 29, 29, 23  (note: 11 331 orbits,
                                   all below the bound of the theorem of moments)
  python censo.py NMIN NMAX PMIN PMAX
Public domain (CC0).
"""
import sys
import time

from P3a import centrados, delta_red, es_primo, gamma_T, s_exacto, semigrupo

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if sys.argv[1:] == ["bajo"]:
    RANGOS = [(7, 11, 29), (8, 11, 29), (9, 11, 23)]
    NMIN, NMAX, PMIN, PMAX = 7, 9, 11, 29
else:
    a = [int(x) for x in sys.argv[1:]]
    NMIN, NMAX, PMIN, PMAX = (a + [3, 6, 5, 31][len(a):])[:4]
    RANGOS = [(n, PMIN, PMAX) for n in range(NMIN, NMAX + 1)]

t0 = time.time()
total = mal = 0
maxd = {}
for p in range(PMIN, PMAX + 1):
    if not (p > 2 and es_primo(p)):
        continue
    for n in [n for n, lo, hi in RANGOS if lo <= p <= hi and n <= p - 2]:
        orbs = centrados(p, n)
        for T in orbs:
            d = delta_red(T, p)
            gG = semigrupo(gamma_T(T, p))["genero"]
            gS = semigrupo(s_exacto(T, p))["genero"]
            if not (d == gG == gS):
                mal += 1
                print("   DISAGREEMENT p=%d T=%s delta=%d g(Gamma)=%d g(S)=%d" % (p, T, d, gG, gS))
            maxd[n] = max(maxd.get(n, 0), d)
        total += len(orbs)
        print("p=%2d n=%d: %5d orbits%s" % (p, n, len(orbs), "   (below the bound (n-1)(n-2) = %d)"
                                            % ((n - 1) * (n - 2)) if p <= (n - 1) * (n - 2) else ""))
        sys.stdout.flush()
print("TOTAL orbits (%d <= n <= %d, %d <= p <= %d): %d ; disagreements delta / g(Gamma_T) / g(S_T): %d"
      % (NMIN, NMAX, PMIN, PMAX, total, mal))
print("largest delta per n:", dict(sorted(maxd.items())))
print("time %.0f s" % (time.time() - t0))
