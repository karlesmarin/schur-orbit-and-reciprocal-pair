# -*- coding: utf-8 -*-
# AUDITORIA TIPOGRAFICA DE LAS FORMULAS, ANTES DE PUBLICAR.
#
# Busca los fallos que compilan sin protestar y por eso no salen en el log.  La PRIMERA version de
# este control dio 38 avisos y 37 eran falsos positivos mios:
#   - \rightsquigarrow y \leftrightarrow contienen "left"/"right" como subcadena;
#   - una formula partida en dos lineas del fuente tiene $ impares en cada una, y es legitima;
#   - 2^rr! se compone BIEN (2^r por r!), que es justo lo que se quiere;
#   - \frac12 y \frac kt son TeX valido.
# Un control que cria lobos esconde el fallo real, asi que se afinaron los cuatro.
#
# Lo que se comprueba ahora:
#   F1  \left sin \right, contando SOLO los que van seguidos de delimitador.
#   F2  llaves descompensadas en una formula.
#   F3  $ impares en el PARRAFO (no en la linea).
#   F4  indice o exponente de dos o mas DIGITOS sin llaves: x^12 se compone x^1 2.
#   F6  operadores pegados.
#   F7  \bigl/\bigr descompensados.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _formaudit.py

import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
s = s.split(r"\begin{thebibliography}")[0]
s = re.sub(r"(?m)^%.*$", "", s)
lineas = s.split("\n")

avisos = []


def avisa(i, que, ctx):
    avisos.append((i, que, re.sub(r"\s+", " ", ctx).strip()[:118]))


# F3 por PARRAFO: una formula puede ocupar varias lineas del fuente.
inicio, dolares, buf = 1, 0, []
for i, L in enumerate(lineas, start=1):
    if not L.strip():
        if dolares % 2:
            avisa(inicio, "F3 dolares impares en el parrafo", " ".join(buf))
        inicio, dolares, buf = i + 1, 0, []
        continue
    dolares += L.replace(r"\$", "").count("$")
    buf.append(L)
if dolares % 2:
    avisa(inicio, "F3 dolares impares en el parrafo", " ".join(buf))

# formulas en linea, ya unidas por parrafo
texto = re.sub(r"\n(?!\s*\n)", " ", s)
for m in re.finditer(r"\$([^$]+)\$", texto):
    f = m.group(1)
    i = s[:s.find(f)].count("\n") + 1 if f in s else 0
    izq = len(re.findall(r"\\left[({\[|.]", f))
    der = len(re.findall(r"\\right[)}\]|.]", f))
    if izq != der:
        avisa(i, "F1 left/right descompensados", f)
    if f.count(r"\bigl") != f.count(r"\bigr"):
        avisa(i, "F7 bigl/bigr descompensados", f)
    if f.count("{") != f.count("}"):
        avisa(i, "F2 llaves descompensadas", f)
    for mm in re.finditer(r"[_^]\d{2,}", f):
        avisa(i, "F4 indice de varios digitos sin llaves: %s" % mm.group(0), f)
    if re.search(r"(?<![+\-])\+\+|(?<!-)--(?!-)", f):
        avisa(i, "F6 operadores pegados", f)

print("AVISOS TIPOGRAFICOS EN FORMULAS")
print("=" * 100)
if not avisos:
    print("  ninguno")
for (i, que, ctx) in avisos:
    print("  L%-5d %-40s %s" % (i, que[:40], ctx))
print("")
print("  total: %d" % len(avisos))
