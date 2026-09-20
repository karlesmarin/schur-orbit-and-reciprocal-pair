# -*- coding: utf-8 -*-
r"""nivel_rango.py -- level-rank SU(n)_k <-> SU(k)_n at prime height, column by column.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

For SU(3)_4 <-> SU(4)_3, SU(3)_8 <-> SU(8)_3 and SU(5)_6 <-> SU(6)_5 (662 columns in all):
  * the index [O_mu : A_mu] of EVERY column, computed as a determinant of integer lattices inside
    Z[zeta_(n(n+k))] (P3a.orden_columna: it does not use the spectrum or the semigroup), against
    p^(phi(m) delta) with delta the number of gaps of S_T (Corollary 'the index of a column');
  * the classes (h, gaps of S_T) on both sides: their counts must stand in the ratio n/k
    (Proposition 'level-rank').
Public domain (CC0).
"""
import sys
import time
from collections import Counter
from fractions import Fraction

from P3a import columna, orden_columna, pesos, phi

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
t0 = time.time()
total = mal = 0
clases = {}
for n, k in ((3, 4), (4, 3), (3, 8), (8, 3), (5, 6), (6, 5)):
    p = n + k
    cnt = Counter()
    malos = 0
    for mu in pesos(n, k):
        c = columna(n, k, mu)
        cnt[(c["h"], tuple(c["S"]["huecos"]))] += 1
        r = orden_columna(n, k, mu)
        total += 1
        if r["OA"] != c["indice"]:
            malos += 1
            print("   index by lattices != p^(phi(m) delta):", n, k, mu, r["OA"], c["indice"])
    mal += malos
    clases[(n, k)] = cnt
    print("SU(%d)_%d p=%d: %d columns; (h, gaps) -> columns: %s | index disagreements: %d"
          % (n, k, p, sum(cnt.values()), dict(sorted(cnt.items())), malos))
    sys.stdout.flush()
for n, k in ((3, 4), (3, 8), (5, 6)):
    a, b = clases[(n, k)], clases[(k, n)]
    ok = all(Fraction(a[x], b[x]) == Fraction(n, k) for x in set(a) | set(b) if b[x]) and set(a) == set(b)
    print("SU(%d)_%d vs SU(%d)_%d: every class in the ratio n/k = %d/%d: %s" % (n, k, k, n, n, k, ok))
print("columns: %d ; [O_mu:A_mu] = p^(phi(m) delta) failed in %d ; time %.0f s" % (total, mal, time.time() - t0))
