# -*- coding: utf-8 -*-
# Longitud del abstract, contada COMO LA CUENTA arXiv: caracteres LITERALES del texto que se pega
# en el campo, no del texto renderizado.  El tope es 1920 y ya nos mordio una vez -- ver la nota de
# trampas de arXiv.
#
# CALIBRADO contra la v2 del Paper I, que entro: el abstract pegado de PARA_SUBIR_v2.txt mide 1908
# caracteres contando `\mu_t`, `\lambda`, `\ge` y demas TAL CUAL.  arXiv NO expande nada: cuenta la
# cadena.  Una version anterior de este guion borraba todo `\xxx` antes de medir y por eso daba
# ~160 caracteres de menos -- es decir, decia que cabia cuando no cabia.
#
# Las dos unicas transformaciones legitimas son las que hay que hacer ANTES de pegar:
#   * `\emph{X}` -> `X`.  arXiv no interpreta macros fuera de $...$ y lo imprimiria literal.
#   * las macros NUESTRAS (`\zt`, `\Sp`, `\spc`, ...) -> su definicion.  arXiv no las conoce.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _abslen.py
import io
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# macros propias del preambulo que hay que deletrear antes de pegar
NUESTRAS = {r"\zt": r"\mu_t", r"\Sp": r"\operatorname{Sp}", r"\spc": r"\mathrm{sp}",
            r"\ZZ": r"\mathbb{Z}", r"\CC": r"\mathbb{C}", r"\Newt": r"\operatorname{Newt}",
            r"\sgn": r"\operatorname{sgn}", r"\core": r"\operatorname{core}",
            r"\quo": r"\operatorname{quot}", r"\rk": r"\operatorname{rank}"}

TOPE = 1920

s = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
a = s.split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]

crudo = " ".join(a.split())
print("en el .tex, tal cual        : %d" % len(crudo))

t = re.sub(r"\\emph\{([^}]*)\}", r"\1", crudo)
for k in sorted(NUESTRAS, key=len, reverse=True):
    t = re.sub(re.escape(k) + r"(?![a-zA-Z])", NUESTRAS[k].replace("\\", "\\\\"), t)
t = " ".join(t.split())

margen = TOPE - len(t)
print("LITERAL, como cuenta arXiv  : %d   (tope %d, margen %d)%s"
      % (len(t), TOPE, margen, "" if margen >= 0 else "   !! SE PASA"))

restantes = sorted(set(re.findall(r"\\[a-zA-Z]+", t)) & set(NUESTRAS))
if restantes:
    print("  macros nuestras sin expandir: %s" % ", ".join(restantes))
if r"\emph" in t:
    print("  QUEDA un \\emph: arXiv lo imprimira literal")

print("")
print("--- el texto exacto que hay que pegar en el campo abstract ---")
print(t)
