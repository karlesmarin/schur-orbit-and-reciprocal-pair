# -*- coding: utf-8 -*-
# CRUCE DE LAS DOS TABLAS.  La tabla de verificacion dice el estado con un marcador de color; la de
# Attribution lo dice con una palabra.  Nada garantizaba que coincidieran, y hoy no coincidian en
# tres filas: cor:oddsign y rem:gkrs figuraban como "verified" en una y "proved" en la otra.
# Este control empareja las dos por la etiqueta \ref y compara.
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()

i = s.find(r"\label{sec:verif}")
j = s.find(r"\section{Attribution}")
k = s.find(r"\subsection*{What is ours}")
m = s.find(r"\section{Open problems}")
verif, ours = s[i:j], s[k:m]


def filas(bloque):
    out = {}
    for linea in bloque.split(r"\\"):
        if "&" not in linea:
            continue
        etiquetas = re.findall(r"\\(?:eq)?ref\{([a-z]+:[A-Za-z0-9]+)\}", linea)
        if not etiquetas:
            continue
        txt = " ".join(linea.split())
        if re.search(r"\b(open|conditional)\b", txt):
            est = "open"          # una conjetura y su evidencia: la tabla mide, Attribution la deja abierta
        # "not proved" contiene "proved": contarlo como probado convierte una fila HONESTA en una
        # discrepancia, y tres filas de conjeturas salieron marcadas asi.  Un control que cria lobos
        # esconde el fallo real, que es justo lo que este cruce existe para encontrar.
        elif r"\stproved" in txt or (re.search(r"\bproved\b", txt)
                                     and not re.search(r"\bnot\s+proved\b", txt)):
            est = "proved"
        elif r"\stext" in txt:
            est = "external"
        elif r"\stverif" in txt or re.search(r"\b(verified|measured|computed)\b", txt):
            est = "verified"
        else:
            est = "?"
        if est == "open":         # el estado "open" convive con una fila de evidencia medida
            out.setdefault("__open__", set()).add(est)
            for e in etiquetas:
                out.setdefault(e, set()).update({"open", "verified"})
            continue
        for e in etiquetas:
            out.setdefault(e, set()).add(est)
    return out


V, O = filas(verif), filas(ours)
comunes = sorted(set(V) & set(O))
print("etiquetas en la tabla de verificacion : %d" % len(V))
print("etiquetas en What is ours             : %d" % len(O))
print("comunes a las dos                     : %d" % len(comunes))
print()
malas = []
for e in comunes:
    a, b = V[e], O[e]
    if not (a & b):
        malas.append((e, sorted(a), sorted(b)))
if malas:
    print("DISCREPANCIAS entre las dos tablas:")
    for e, a, b in malas:
        print("   %-24s verificacion=%s   attribution=%s" % (e, a, b))
else:
    print("sin discrepancias: toda etiqueta comun tiene el mismo estado en las dos tablas")

print()
solo_o = sorted(set(O) - set(V))
print("en Attribution y NO en la tabla de verificacion (%d):" % len(solo_o))
for e in solo_o:
    print("   %-24s %s" % (e, sorted(O[e])))
