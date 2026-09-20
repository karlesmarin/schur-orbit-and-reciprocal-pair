# -*- coding: utf-8 -*-
r"""maximales.py -- the number of maximal columns: closed formula against enumeration.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Proposition 'the maximal columns':
   N_0(p, n) = n * sum_{h | p-1, h > 1, n mod h in {0,1}} [ C((p-1)/h, n//h) - Z_{p,(p-1)/h}(n//h) ],
Z_{p,e}(d) = number of d-subsets of the subgroup of order e of F_p^x with sum 0.  Against the
enumeration of every centred spectrum with S_T = N_0 (exact semigroup), each standing for n
columns, for p = 17 and 19 and every 2 <= n <= p-2.
Public domain (CC0).
"""
import sys

from P3a import maximales, maximales_enumeradas

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
bien = tot = 0
for p in (17, 19):
    for n in range(2, p - 1):
        f, en = maximales(p, n), maximales_enumeradas(p, n)
        tot += 1
        bien += f == en
        print("p=%d n=%2d: formula %5d  enumeration %5d  %s" % (p, n, f, en, "ok" if f == en else "DIFFERENT"))
        sys.stdout.flush()
print("formula = enumeration in %d of %d pairs (p, n)" % (bien, tot))
