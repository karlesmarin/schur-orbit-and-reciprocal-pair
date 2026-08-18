# -*- coding: utf-8 -*-
# Inventario de afirmaciones con marcador de estado, para el repaso de principio a fin.
# Saca cada \stverif / \stproved / \stext con su frase y con los numeros que arrastra, para poder
# cotejarlos contra los volcados de gates/ uno a uno.  Dos fallos de hoy fueron numericos y los dos
# los caza esto: un numero que no sale de su volcado, y un verbo de estado que ya no corresponde.
import io
import re
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
s = s.split(r"\begin{document}")[1]

# se trocea en "frases" por punto seguido, conservando el marcador
texto = " ".join(s.split())
frases = re.split(r"(?<=[.;])\s+", texto)

MARCA = re.compile(r"\\st(verif|proved|ext)\b")
NUM = re.compile(r"(?<![\w\\])(\d[\d\,\\]*\d|\d)(?![\w])")

cuenta = Counter()
filas = []
for f in frases:
    m = MARCA.search(f)
    if not m:
        continue
    cuenta[m.group(1)] += 1
    nums = [n.replace("\\,", "").replace(",", "") for n in NUM.findall(f)]
    limpia = re.sub(r"\\st(verif|proved|ext)\b", "", f)
    limpia = re.sub(r"\\[a-zA-Z]+", " ", limpia)
    limpia = limpia.replace("{", "").replace("}", "").replace("$", "")
    limpia = " ".join(limpia.split())
    filas.append((m.group(1), nums, limpia))

print("marcadores por tipo: %s   (total %d)" % (dict(cuenta), sum(cuenta.values())))
print("=" * 110)
for i, (tipo, nums, f) in enumerate(filas, 1):
    print("%3d [%-7s] %s" % (i, tipo, f[:250]))
    if nums:
        print("        numeros: %s" % ", ".join(nums[:14]))
