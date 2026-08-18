# -*- coding: utf-8 -*-
# AUDITORIA DE NUMEROS ANTES DE PUBLICAR.
#
# Dos peligros distintos y los dos reales hoy:
#   (1) un numero que se quedo obsoleto al ampliar una medida -- el barrido del determinante paso de
#       4 configuraciones a 12, y los senuelos cambiaron de denominador;
#   (2) DOS POBLACIONES MEZCLADAS: las estadisticas de cancelacion (1349 fibras, 333 multiples) son
#       de 4 configuraciones, y las del determinante son de 12.  Si el texto las pone juntas sin
#       decirlo, el lector suma peras y manzanas.
#
# Saca cada afirmacion "N de M" / "N/M" con su linea y su contexto, agrupadas por el par, para poder
# leer de un golpe si dos sitios que citan el mismo numero hablan de lo mismo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _numaudit.py

import io
import re
import sys
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
s = s.split(r"\begin{thebibliography}")[0]
lineas = s.split("\n")

PAR = re.compile(r"\$?(\d[\d\\,]*)\$?\s*(?:of|/)\s*\$?(\d[\d\\,]*)\$?")


def limpia(x):
    return x.replace("\\,", "").replace(",", "")


pares = defaultdict(list)
for i, L in enumerate(lineas, start=1):
    if L.strip().startswith("%"):
        continue
    for m in PAR.finditer(L):
        a, b = limpia(m.group(1)), limpia(m.group(2))
        if not (a.isdigit() and b.isdigit()):
            continue
        if int(b) == 0 or int(a) > int(b):
            continue
        ctx = re.sub(r"\s+", " ", L).strip()
        pares[(int(a), int(b))].append((i, ctx[:120]))

print("PARES N/M QUE APARECEN EN MAS DE UN SITIO  (candidatos a incongruencia)")
print("=" * 104)
for k in sorted(pares, key=lambda k: -len(pares[k])):
    if len(pares[k]) < 2:
        continue
    print("")
    print("  %d de %d   (%d apariciones)" % (k[0], k[1], len(pares[k])))
    for (i, c) in pares[k]:
        print("     L%-5d %s" % (i, c))

print("")
print("=" * 104)
print("DENOMINADORES QUE APARECEN CON VARIOS NUMERADORES (poblaciones que conviene comparar)")
print("=" * 104)
porden = defaultdict(set)
for (a, b) in pares:
    porden[b].add(a)
for b in sorted(porden, key=lambda b: -len(porden[b])):
    if len(porden[b]) < 3:
        continue
    print("  sobre %-8d : numeradores %s" % (b, sorted(porden[b])))
