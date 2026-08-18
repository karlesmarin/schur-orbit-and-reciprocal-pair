# -*- coding: utf-8 -*-
# Inventario de lo que queda: toda fila de Attribution que NO diga "proved", y los problemas abiertos.
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
k = s.find(r"\subsection*{What is ours}")
m = s.find(r"\section{Open problems}")


def limpia(x):
    x = re.sub(r"\\[a-zA-Z]+", " ", x)
    return " ".join(x.replace("{", "").replace("}", "").replace("$", "").split())


print("=== FILAS DE 'What is ours' QUE NO SON 'proved' ===")
for linea in s[k:m].split("\\\\"):
    if "&" not in linea:
        continue
    txt = " ".join(linea.split())
    if re.search(r"\bproved\b", txt):
        continue
    t2 = limpia(txt)
    if len(t2) > 15:
        print("  -", t2[:160])

print()
print("=== PROBLEMAS ABIERTOS ===")
for mm in re.finditer(r"\\begin\{problem\}\[([^\]]*)\]\\label\{([^}]*)\}", s):
    print("  - %-22s %s" % (mm.group(2), mm.group(1)))
print()
print("=== CONJETURAS ===")
for mm in re.finditer(r"\\begin\{conjecture\}\[([^\]]*)\]\\label\{([^}]*)\}", s):
    print("  - %-22s %s" % (mm.group(2), mm.group(1)))
