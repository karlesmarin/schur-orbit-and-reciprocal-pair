# -*- coding: utf-8 -*-
# AUDITORIA DE AMBIGUEDAD.  .Usa el paper la misma palabra para dos cosas distintas?
#
# POR QUE.  Un paper de 62 paginas acumula polisemia sin que nadie lo note, y la polisemia mas
# peligrosa es la que uno mismo introduce tarde: hoy entraron "class" (de matrices salvo permutacion
# con signo) y "fibre", y el paper ya usaba "class" en dos sentidos.
#
# Lo que hace: para cada palabra de riesgo, saca el sustantivo que la precede o la sigue -- que es lo
# que fija el sentido -- y agrupa.  Si una palabra aparece con collocations muy distintas, es
# candidata a desambiguar.
#
# No decide nada por si mismo: imprime para que se lea.  Un control que decidiera "esto es ambiguo"
# seria peor que inutil.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _ambiguity.py

import io
import re
import sys
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RIESGO = ["class", "classes", "filter", "filters", "regular", "primitive", "unit", "units",
          "fibre", "fibres", "top", "weight", "sign", "order", "orbit", "point", "level"]

s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
s = s.split(r"\begin{thebibliography}")[0]
# fuera comentarios y matematicas, que no son prosa
s = re.sub(r"(?m)^%.*$", " ", s)
s = re.sub(r"\$[^$]*\$", " MATH ", s)
s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?", " ", s)
s = re.sub(r"\s+", " ", s)

pal = re.findall(r"[A-Za-z][A-Za-z-]*", s)
idx = defaultdict(Counter)
for i, w in enumerate(pal):
    lw = w.lower()
    if lw in RIESGO:
        antes = pal[i - 1].lower() if i else "^"
        desp = pal[i + 1].lower() if i + 1 < len(pal) else "$"
        idx[lw]["%s _ %s" % (antes, desp)] += 1

print("PALABRAS DE RIESGO Y SUS ACOMPANANTES  (sentido = quien va al lado)")
print("=" * 96)
for w in RIESGO:
    c = idx.get(w)
    if not c:
        continue
    tot = sum(c.values())
    print("")
    print("  %-10s %4d apariciones, %3d contextos distintos" % (w, tot, len(c)))
    for k, n in c.most_common(7):
        print("        %-40s %3d" % (k, n))

print("")
print("=" * 96)
print("FRASES QUE ARRANCAN CON UN PRONOMBRE SIN ANTECEDENTE VISIBLE")
print("=" * 96)
frases = re.split(r"(?<=[.!?]) +", s)
n = 0
for f in frases:
    m = re.match(r"^(It|This|That|These|Those|They)\s+(is|are|was|were|does|do|has|have)\b", f)
    if m:
        n += 1
        if n <= 12:
            print("   %s" % f[:110])
print("   total: %d" % n)
