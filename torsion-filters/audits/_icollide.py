# -*- coding: utf-8 -*-
# LA COLISION DE LA UNIDAD IMAGINARIA.
#
# El paper usa i como UNIDAD IMAGINARIA (zeta = e^{2 pi i / t}) y como INDICE DE FILA, y en dos
# formulas aparecian las dos cosas a la vez: "2i\sin(2\pi i a_j/t)" y "(2i)^n\det(\sin(2\pi ij/t))".
# El lector no puede saber cual es cual, y las dos formulas sostienen lemas.
#
# Este control busca formulas donde convivan un "i" suelto multiplicativo y un "i" de indice.
# Heuristica: una formula que contenga a la vez  \sin(...i...)  o  ^{i...}  o  _{i}  y un "i" pegado
# a un numero (2i, 4i) o precedido de ( .  No decide: senala para mirar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _icollide.py

import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
s = s.split(r"\begin{thebibliography}")[0]
s = re.sub(r"(?m)^%.*$", "", s)
texto = re.sub(r"\n(?!\s*\n)", " ", s)

IMAG = re.compile(r"(?<![A-Za-z\\])\d\s*i(?![A-Za-z])|\(\s*\d\s*i\s*\)")
INDICE = re.compile(r"[_^]\{?i\b|\bi\s*[a-z]?j\b|_\{i")

avisos = []
for m in re.finditer(r"\$([^$]+)\$", texto):
    f = m.group(1)
    if IMAG.search(f) and INDICE.search(f):
        avisos.append(f)

print("FORMULAS DONDE 'i' PODRIA SER A LA VEZ UNIDAD IMAGINARIA E INDICE")
print("=" * 96)
if not avisos:
    print("  ninguna")
for f in avisos:
    print("   %s" % f[:110])
print("")
print("  total: %d" % len(avisos))
print("")
print("  NOTA: es una heuristica y senala para MIRAR.  Las dos que encontro la lectura de fondo")
print("  --- '2i sin' en la prueba de lem:T y '(2i)^n' en la de cor:galoisquad --- ya estan")
print("  corregidas a \\sqrt{-1} con el indice nombrado.")
