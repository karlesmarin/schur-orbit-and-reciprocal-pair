# -*- coding: utf-8 -*-
r"""ley_2p_largo.py -- ¿emerge el -1/p^2, o se queda tapado por el error?

En ley_2p.py se comprobaron los cuatro pasos de su demostracion y salio que

    (F_n(p) - 2/p + 1/p^2) p^2  =  p^3 * E(n,p) / C(p,n),     E = e_{12} + e_{13} - e_{123},

con e_J el error signado de C_J respecto de C(p,n)/p^{|J|}.  Para n = 8 y p = 43 eso vale 0.887,
o sea DEL MISMO TAMANO que el propio termino -1 que se quiere ver.  No es que la formula falle:
es que en ese rango el -1 esta tapado.

La cuenta dice que deberia emerger: los errores se comportan como C p^{n/2} con C ~ 4e-4, luego
el termino de arriba es del orden de n! * C * p^{3-n/2}, que para n = 8 es ~16/p.  A p = 43 eso
predice ~0.37 y se midio 0.887; a p = 89 deberia bajar a ~0.18.  Aqui se mide de verdad en vez
de extrapolar, que es lo que distingue una comprobacion de una excusa.

El coste es el censo de tres momentos simultaneos, O(n p^4), asi que se corre solo para n = 8 y
en segundo plano.

    python ley_2p_largo.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from math import comb

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo                                              # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 89
N = 8


def cuenta(p, pots):
    dp = np.zeros((N + 1,) + (p,) * len(pots), dtype=np.int64)
    dp[(0,) + (0,) * len(pots)] = 1
    for t in range(p):
        des = tuple(pow(t, e, p) for e in pots)
        for k in range(N - 1, -1, -1):
            b = dp[k]
            for ax, d in enumerate(des):
                b = np.roll(b, d, axis=ax)
            dp[k + 1] += b
    return int(dp[(N,) + (0,) * len(pots)])


print("n = %d.  OJO AL NOMBRE: lo que sale de aqui es" % N)
print("   Ftilde = (C12 + C13 - C123)/(C(p,n)/p),  la proporcion en el LOCUS  M_2 M_3 = 0,")
print("NO es F_n(p), que es la proporcion con delta >= 2.  Las dos coinciden solo donde h = 1.")
print("La correccion exacta de los estabilizadores esta en ley_2p_estabilizadores.py, que ademas")
print("imprime la tabla buena.  La nota publico ESTA columna llamandola F_8, y era otra variable.")
print("%-5s %-10s %-10s %-9s %-13s %-13s %-12s"
      % ("p", "C12", "C13", "C123", "E firmado", "(F-2/p)p^2", "resto vs -1"))
for p in range(29, PMAX + 1):
    if not es_primo(p) or p <= N + 1:
        continue
    c12 = cuenta(p, (1, 2))
    c13 = cuenta(p, (1, 3))
    c123 = cuenta(p, (1, 2, 3))
    tot = comb(p, N)
    cent = tot // p
    E = (c12 - tot / p ** 2) + (c13 - tot / p ** 2) - (c123 - tot / p ** 3)
    F = (c12 + c13 - c123) / cent
    v = (F - 2 / p) * p * p
    print("%-5d %-10d %-10d %-9d %-13.1f %-13.4f %-12.4f"
          % (p, c12, c13, c123, E, v, v + 1))
    sys.stdout.flush()
print("")
print("Si la columna 'resto vs -1' tiende a 0, el -1/p^2 queda confirmado.")
print("Si se estanca lejos de 0, el refinamiento no se sostiene en este rango y hay que decirlo.")
