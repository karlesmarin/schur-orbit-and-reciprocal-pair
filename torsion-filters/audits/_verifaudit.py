# -*- coding: utf-8 -*-
# Barre la tabla de verificacion fila a fila: busca filas duplicadas, filas sin estado, y filas cuyo
# enunciado nombra un resultado que el paper ya marca como probado (que es el desfase de hoy: una
# fila decia "verified" de algo que habia pasado a teorema por la manana).
import io
import re
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
i = s.find(r"\section{Verification}")
if i < 0:
    i = s.find(r"\label{sec:verif}")
j = s.find(r"\section{Attribution}")
tabla = s[i:j]

filas = []
for linea in tabla.split(r"\\"):
    if "&" not in linea:
        continue
    cols = [c.strip() for c in linea.split("&")]
    if len(cols) < 2:
        continue
    limpio = []
    for c in cols:
        c = re.sub(r"\\st(verif|proved|ext)\b", lambda m: "[" + m.group(1) + "]", c)
        c = re.sub(r"\\[a-zA-Z]+", " ", c)
        c = c.replace("{", "").replace("}", "").replace("$", "")
        limpio.append(" ".join(c.split()))
    filas.append(limpio)

print("filas con celdas: %d" % len(filas))
est = Counter()
sin_estado = []
for f in filas:
    m = [c for c in f if c.startswith("[")]
    if m:
        est[m[-1]] += 1
    else:
        sin_estado.append(f[0][:70])
print("estados: %s" % dict(est))
if sin_estado:
    print("filas SIN marcador de estado (%d):" % len(sin_estado))
    for f in sin_estado[:12]:
        print("   -", f)

print()
print("duplicados por enunciado:")
c = Counter(f[0][:60] for f in filas if f[0])
dup = [(k, v) for k, v in c.items() if v > 1 and len(k) > 12]
print("   %s" % (dup if dup else "ninguno"))

print()
print("filas, enunciado -> resultado:")
for n, f in enumerate(filas, 1):
    print("%3d  %-58s | %s" % (n, f[0][:58], " | ".join(f[1:])[:80]))
