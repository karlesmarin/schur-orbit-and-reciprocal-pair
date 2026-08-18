# -*- coding: utf-8 -*-
# .DE DONDE SALEN LOS 313 PARES?   17 de agosto de 2026.
#
# El paper introduce "313 de 313" en rem:divided sin decir como se relaciona con las 333 fibras
# multiples (322 de dos terminos + 11 de cuatro).  Ningun numero enunciado da 313, y un lector que
# intente cuadrarlo no puede.  Aqui se DESGLOSA, por tamano de fibra, exactamente el mismo conteo
# que hace _probe_toggle.py: pares NO ordenados de terminos de una misma fibra que caen en puntos
# del soporte DISTINTOS (c1 != c2) y que se CANCELAN (s1.nu[c1] == -s2.nu[c2]).
#
# LO QUE SE MIDE
#   P1  el total, que tiene que reproducir el 313 de _probe_toggle.py   (FATAL)
#   P2  el desglose por tamano de fibra (2 terminos / 4 terminos)
#   P3  cuantas fibras de dos terminos caen en el MISMO punto canonico -- las 41 del texto
#   P4  y el reparto de pares por fibra de cuatro terminos
#
# CONTROL
#   C0  P1 es fatal: si no da 313, este guion no esta contando lo que cuenta _probe_toggle.py y el
#       desglose no vale para explicarlo.
#   C1  P2 tiene que sumar P1, y las cuentas de fibras tienen que reproducir 322, 11 y 333, que
#       vienen de divided_differences.py.  Si no, las dos poblaciones no son la misma.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_313.py
import itertools
import json
import sys
from collections import Counter

from divided_differences import CASOS
from _probe_toggle import objetos, progresiones

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

pares_por_tam = Counter()
fibras_por_tam = Counter()
mismo_punto = Counter()     # fibras cuyos DOS terminos caen en el mismo punto canonico
pares_por_fibra4 = Counter()

for (t, r, cota) in CASOS:
    Rp = (t - 1) // 2 + r
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu, obj, V, clases = objetos(list(Lam), t, r)
        if not nu:
            continue
        for X, lst in progresiones(nu, t, r).items():
            n = len(lst)
            if n < 2:
                continue
            fibras_por_tam[n] += 1
            if n == 2 and lst[0][2] == lst[1][2]:
                mismo_punto[n] += 1
            k = 0
            for p in range(n):
                for q in range(p + 1, n):
                    (k1, s1, c1) = lst[p]
                    (k2, s2, c2) = lst[q]
                    if c1 == c2:
                        continue
                    if s1 * nu[c1] != -s2 * nu[c2]:
                        continue
                    k += 1
            pares_por_tam[n] += k
            if n == 4:
                pares_por_fibra4[k] += 1

total = sum(pares_por_tam.values())
fib = sum(fibras_por_tam.values())

print("=" * 96)
print("DE DONDE SALEN LOS 313 PARES")
print("=" * 96)
print("")
print("  fibras con 2 o mas terminos           : %d   %s" % (fib, dict(sorted(fibras_por_tam.items()))))
print("  de las de 2, con los dos terminos en")
print("  el MISMO punto canonico (antisimetria): %d" % mismo_punto[2])
print("")
print("  PARES que cancelan en puntos DISTINTOS, por tamano de fibra:")
for n in sorted(pares_por_tam):
    print("     fibras de %d terminos : %4d pares" % (n, pares_por_tam[n]))
print("     TOTAL                : %4d" % total)
print("")
print("  reparto de pares por fibra de 4 terminos: %s" % dict(sorted(pares_por_fibra4.items())))
print("")
print("  P1  FATAL  el total reproduce el 313 de _probe_toggle.py : %s" % (total == 313))
print("  C1  las fibras reproducen 322 + 11 = 333                 : %s"
      % (fibras_por_tam.get(2) == 322 and fibras_por_tam.get(4) == 11 and fib == 333))
print("  C1  y las de dos en el mismo punto son las 41 del texto   : %s" % (mismo_punto[2] == 41))
assert total == 313, "*** no reproduce el 313: este conteo NO es el de _probe_toggle.py ***"
json.dump({"total": total, "pares_por_tam": dict(pares_por_tam),
           "fibras_por_tam": dict(fibras_por_tam), "mismo_punto_2": mismo_punto[2],
           "pares_por_fibra4": dict(pares_por_fibra4)},
          open("_probe_313_DUMP.json", "w"), indent=1)
print("")
print("=" * 96)
print("DONE")
