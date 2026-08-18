# -*- coding: utf-8 -*-
# .CITAN LAS DOS EDICIONES LO MISMO, EN LOS MISMOS SITIOS?   18 de agosto de 2026.
#
# POR QUE.  _transcheck.py compara CONJUNTOS de claves citadas y de \bibitem, y con eso da 0
# discrepancias aunque una edicion cite una fuente tres veces y la otra dos.  La pregunta es si las
# referencias son las mismas entre las dos ediciones; esto lo mide de verdad:
#
#   B1  la LISTA de \bibitem --- etiqueta y clave, EN ORDEN --- tiene que ser identica.
#   B2  el MULTICONJUNTO de claves citadas en el cuerpo tiene que ser identico: no basta el conjunto,
#       porque una cita de menos es una atribucion que falta en una de las dos ediciones.
#   B3  y donde difiera, se nombra la clave y se imprime el contexto de cada aparicion, que es lo
#       unico que permite decidir cual de las dos esta bien.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _bibparity.py
import collections
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CITA = re.compile(r"\\cite\{([^}]*)\}")
ITEM = re.compile(r"\\bibitem\[([^\]]*)\]\{([^}]*)\}")


def lee(f):
    s = io.open(f, encoding="utf-8").read()
    cuerpo = s.split(r"\begin{thebibliography}")[0]
    claves = []
    for m in CITA.finditer(cuerpo):
        for k in m.group(1).split(","):
            claves.append((k.strip(), cuerpo.count("\n", 0, m.start()) + 1))
    return collections.Counter(k for k, _ in claves), claves, ITEM.findall(s), cuerpo


cEN, lEN, bEN, txEN = lee("orbit_pair_ii.tex")
cES, lES, bES, txES = lee("orbit_pair_ii_es.tex")

print("=" * 96)
print("PARIDAD DE REFERENCIAS ENTRE LAS DOS EDICIONES")
print("=" * 96)
print("")
print("  B1  \\bibitem: EN %d, ES %d, lista identica y en el mismo orden : %s"
      % (len(bEN), len(bES), bEN == bES))
if bEN != bES:
    for i, (a, b) in enumerate(zip(bEN, bES)):
        if a != b:
            print("        primera diferencia en la posicion %d: EN %s / ES %s" % (i, a, b))
            break

print("  B2  citas en el cuerpo: EN %d, ES %d" % (sum(cEN.values()), sum(cES.values())))
dif = sorted(k for k in set(cEN) | set(cES) if cEN[k] != cES[k])
print("      claves con distinto numero de citas: %d" % len(dif))
for k in dif:
    print("")
    print("      %-12s  EN %d   ES %d" % (k, cEN[k], cES[k]))
    for etq, txt, lst in (("EN", txEN, lEN), ("ES", txES, lES)):
        for kk, ln in lst:
            if kk == k:
                i = txt.find("\n".join(txt.split("\n")[ln - 1:ln]))
                ctx = " ".join(txt.split("\n")[ln - 1].split())
                print("         %s l.%-5d %s" % (etq, ln, ctx[:110]))

print("")
print("  VEREDICTO: las dos ediciones citan lo mismo : %s" % (not dif and bEN == bES))
print("")
print("=" * 96)
print("DONE")
