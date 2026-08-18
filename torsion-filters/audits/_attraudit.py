# -*- coding: utf-8 -*-
# AUDITORIA DE ATRIBUCION.  .Nos estamos colgando medallas de otros?
#
# POR QUE.  No atribuirnos formulas que ya existen y son de otros.  El paper tiene dos tablas -- lo que NO es nuestro y lo que SI --
# y nada garantizaba que cada fuente citada apareciera en la primera.  Esto lo comprueba.
#
# CONTROLES
#   A1  toda clave con \bibitem se cita en el cuerpo (lo mira _citeaudit.py, aqui solo se recuerda).
#   A2  toda clave citada aparece en la tabla de lo que NO es nuestro, o bien en una fila de "lo que
#       si" que declara explicitamente que parte se toma prestada.  Una clave que solo vive en el
#       cuerpo es una deuda sin declarar.
#   A3  ninguna fila de "lo que si es nuestro" dice "proved" a secas si su enunciado nombra una
#       fuente: si hay cita en la fila, tiene que haber una clausula de prestamo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _attraudit.py

import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()

CITE = re.compile(r"\\cite\{([^}]+)\}")
BIB = re.compile(r"\\bibitem\[[^\]]*\]\{([^}]+)\}")

i_ours = s.index(r"\subsection*{What is ours}")
i_attr = s.index(r"\section{Attribution}")
prestado = s[i_attr:i_ours]          # todo lo anterior a "What is ours" dentro de Attribution
nuestro = s[i_ours:]

claves_bib = sorted(set(BIB.findall(s)))
en_prestado = set()
for m in CITE.finditer(prestado):
    en_prestado.update(k.strip() for k in m.group(1).split(","))
en_nuestro = set()
for m in CITE.finditer(nuestro):
    en_nuestro.update(k.strip() for k in m.group(1).split(","))

print("fuentes con bibitem                       : %d" % len(claves_bib))
print("citadas en la parte de PRESTADO           : %d" % len(en_prestado))
print("citadas en filas de 'What is ours'        : %d" % len(en_nuestro))

huerfanas = [k for k in claves_bib if k not in en_prestado and k not in en_nuestro]
print("")
print("A2  fuentes que NO aparecen en NINGUNA de las dos tablas de Attribution (%d):" % len(huerfanas))
for k in huerfanas:
    n = len(re.findall(r"\\cite\{[^}]*\b" + re.escape(k) + r"\b[^}]*\}", s))
    print("   %-18s  citada %d vez/veces en el cuerpo" % (k, n))
if not huerfanas:
    print("   ninguna: toda fuente citada esta declarada en Attribution")

print("")
print("A3  filas de 'What is ours' que citan una fuente, con su clausula:")
filas = [f for f in nuestro.split("\\\\") if "\\cite{" in f]
for f in filas:
    txt = " ".join(f.split())
    amp = txt.rfind("&")
    izq, der = (txt[:amp], txt[amp + 1:]) if amp > 0 else (txt, "")
    ok = any(w in der.lower() for w in ["is ", "are ", "after", "ours", "not ours", "classical",
                                        "standard", "elementary", "correspondent"])
    print("   [%s] %s" % ("ok " if ok else "!! ", txt[:150]))
