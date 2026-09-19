# -*- coding: utf-8 -*-
r"""II_R_vacio.py -- R_q vacio para q = 3 mod 8?  (y tamanos de R_q en las otras clases), hasta QMAX.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

R_q = {0 <= r < q : T_q + q C(r) = 0, T(r) = 0},  C, T sumas parciales de (r/q) y r (r/q).
Para q = 3 mod 4: T_q = -q h, asi que R_q = {C(r) = h, T(r) = 0}.  Se registra tambien el 'casi':
min |T(r)| sobre los r con C(r) = h (cuanto falta para un elemento), para ver si el vacio es rigido.
Uso: python II_R_vacio.py [QMAX]
"""
import sys
import numpy as np
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
es = np.ones(QMAX + 1, dtype=bool); es[:2] = False
for p in range(2, int(QMAX ** .5) + 1):
    if es[p]:
        es[p * p::p] = False
tam = Counter(); no_vacio = []; casi = []
for q in np.nonzero(es)[0]:
    q = int(q)
    if q < 5:
        continue
    r = np.arange(q, dtype=np.int64)
    L = -np.ones(q, dtype=np.int64); L[0] = 0
    L[(r[1:] * r[1:]) % q] = 1
    C = np.cumsum(L); T = np.cumsum(r * L); Tq = int(T[-1])
    enR = (Tq + q * C == 0) & (T == 0)
    k = int(enR.sum())
    clase = "1mod4" if q % 4 == 1 else ("7mod8" if q % 8 == 7 else "3mod8")
    tam["%s |R|=%d" % (clase, k)] += 1
    if clase == "3mod8":
        if k:
            no_vacio.append((q, list(np.nonzero(enR)[0][:6])))
        sel = (Tq + q * C == 0)
        if sel.any():
            casi.append((q, int(np.abs(T[sel]).min())))
print("primos 5 <= q <= %d" % QMAX)
for c in sorted(tam):
    print("  %-18s %d" % (c, tam[c]))
print("q = 3 mod 8 con R_q NO vacio: %d  %s" % (len(no_vacio), no_vacio[:10]))
casi.sort(key=lambda x: x[1])
print("q = 3 mod 8: los r con C(r) = h existen en %d primos; menor |T(r)| alli: %s"
      % (len(casi), casi[:12]))
