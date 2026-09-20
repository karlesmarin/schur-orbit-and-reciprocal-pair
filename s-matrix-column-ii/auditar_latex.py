# -*- coding: utf-8 -*-
r"""auditar_latex.py -- comandos de LaTeX a los que se les ha comido la contrabarra.

En la pagina 17 llego a imprimirse, literalmente, "extsfA262568".  La causa: en algun parche
\textsf perdio la barra y el \t se convirtio en un TABULADOR, de modo que el fuente decia
"<TAB>extsf{A262568}".  LaTeX compila eso sin una sola queja --- imprime el tabulador como
espacio, la palabra "extsf" como texto y el argumento como texto --- y sale con codigo 0.

Es el fallo mas barato de cometer y el mas caro de ver: no hay aviso, no hay error, y la
paridad EN/ES no lo caza porque compara formulas, no palabras.  Aqui se caza por tres vias:

  (1) un TABULADOR en el fuente: no hay ninguno legitimo en esta nota;
  (2) un nombre de comando conocido pegado a una llave SIN su contrabarra delante;
  (3) las secuencias que delatan la barra comida por el escape de Python o del shell:
      un tabulador, un salto de linea o un retorno seguidos de un nombre de comando.

    python auditar_latex.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))

# los que empiezan por una letra que es tambien un escape de Python: \t \n \r \f \v \b \a
PELIGROSOS = ["textsf", "textbf", "textit", "texttt", "text", "times", "to", "tfrac", "top",
              "newline", "noindent", "nonumber", "ref", "rangle", "right",
              "frac", "footnote", "bigl", "binom", "alpha", "approx", "varepsilon", "vspace"]
CONOCIDOS = PELIGROSOS + ["emph", "cite", "eqref", "label", "begin", "end", "section",
                          "subsection", "medido", "marco", "includegraphics", "caption"]

fallos = []
hechas = []


def afirma(ok, que):
    print("   %-4s %s" % ("ok" if ok else "MAL", que))
    hechas.append(que)
    if not ok:
        fallos.append(que)


for fich in ("P3b.tex", "P3b_es.tex"):
    s = io.open(os.path.join(AQUI, fich), encoding="utf-8").read()
    print("")
    print("%s" % fich)

    # (1) tabuladores: en esta nota no hay ninguno legitimo
    tabs = [i + 1 for i, l in enumerate(s.splitlines()) if "\t" in l]
    afirma(not tabs, "sin tabuladores en el fuente%s"
           % ("" if not tabs else " --- los hay en las lineas %s" % tabs))

    # (2) un comando conocido pegado a { sin contrabarra delante
    sueltos = []
    for m in re.finditer(r"(.)(%s)\{" % "|".join(CONOCIDOS), s):
        if m.group(1) not in "\\" and (m.group(1).isalpha() is False):
            sueltos.append((s[:m.start()].count("\n") + 1, m.group(2), repr(m.group(1))))
    afirma(not sueltos, "ningun comando conocido sin su contrabarra%s"
           % ("" if not sueltos else " --- %s" % sueltos[:6]))

    # (3) la firma exacta de la barra comida por un escape
    escapes = []
    for c, nombre in (("\t", "t"), ("\n", "n"), ("\r", "r"), ("\f", "f"), ("\v", "v")):
        for cmd in CONOCIDOS:
            if not cmd.startswith(nombre):
                continue
            aguja = c + cmd[len(nombre):] + "{"
            for m in re.finditer(re.escape(aguja), s):
                escapes.append((s[:m.start()].count("\n") + 1, "\\" + cmd))
    afirma(not escapes, "ninguna contrabarra comida por un escape%s"
           % ("" if not escapes else " --- %s" % escapes))

    # (4) el control: el criterio caza un caso plantado
    sucio = s[:200] + "\textsf{prueba}" + s[200:]
    caza = "\t" in sucio and re.search(r"\textsf\{", sucio)
    afirma(bool(caza), "control: el criterio caza un \\textsf plantado con la barra comida")

print("")
print("=" * 82)
print("comprobaciones: %d ; fallos: %d" % (len(hechas), len(fallos)))
if not hechas:
    sys.exit("FATAL: cero comprobaciones. No poder medir no es aprobar")
print("VEREDICTO:", "ningun comando roto" if not fallos
      else "COMANDOS ROTOS:\n   - " + "\n   - ".join(fallos))
sys.exit(1 if fallos else 0)
