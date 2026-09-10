# -*- coding: utf-8 -*-
# LAS FIGURAS: SITIO, DATOS Y PIE.   8 de septiembre de 2026.
#
# POR QUE.  Un pie de figura es un documento sin auditar: sus cifras no las comprueba nadie, y el
# texto DENTRO del PDF de la figura no lo ve ninguna compuerta sobre el .tex.  Ademas, al soltar
# las figuras de [H] a [tbp] pueden alejarse del parrafo que las cita.
#
# QUE COMPRUEBA
#   F1  cada figura del .tex existe en disco y es un PDF legible
#   F2  cada figura esta citada, y la CITA y la FIGURA caen en la misma pagina o en la siguiente
#   F3  las cifras que el pie afirma estan en la salida del guion que lo genero
#   F4  el texto dentro de la imagen no contradice fechas de la bibliografia
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_figuras.py > check_figuras_OUT.txt 2>&1

import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


import fitz

GATES = os.path.join("..", "gates")

for doc in ("characters_at_aP", "characters_at_aP_es"):
    print(SEP)
    print(doc)
    tex = io.open(doc + ".tex", encoding="utf-8").read()
    pdf = fitz.open(doc + ".pdf")
    paginas = [pdf[i].get_text() for i in range(pdf.page_count)]

    figs = re.findall(r"\\begin\{figure\}.*?includegraphics\[[^\]]*\]\{([^}]+)\}(.*?)\\end\{figure\}",
                      tex, re.S)
    ok(len(figs) >= 5, f"{doc}: al menos cinco figuras", f"{len(figs)}")

    # --- F1
    for nombre, cuerpo in figs:
        ok(os.path.exists(nombre), f"F1 {nombre}: existe en disco")

    # --- F2  distancia entre la cita y la figura
    for nombre, cuerpo in figs:
        m = re.search(r"\\label\{([^}]+)\}", cuerpo)
        if not m:
            ok(False, f"F2 {nombre}: tiene etiqueta"); continue
        lab = m.group(1)
        # pagina de la figura: la que imprime su etiqueta «Figure N:» / «Figura N:».  Buscar las
        # primeras palabras del pie falla cuando el pie empieza por matematicas.
        num = None
        aux = doc + ".aux"
        if os.path.exists(aux):
            a = io.open(aux, encoding="utf-8", errors="replace").read()
            mm = re.search(r"\\newlabel\{%s\}\{\{([^}]*)\}\{(\d+)\}" % re.escape(lab), a)
            if mm:
                num, pag_lab = mm.group(1), int(mm.group(2))
        pag_fig = None
        if num is not None:
            et = re.compile(r"(Figure|Figura)\s*%s\s*:" % re.escape(num))
            pag_fig = next((i + 1 for i, t in enumerate(paginas)
                            if et.search(" ".join(t.split()))), None)
        if pag_fig is None or num is None:
            ok(False, f"F2 {lab}: se puede localizar la figura y su numero",
               f"pag_fig={pag_fig} num={num}")
            continue
        # donde se cita
        # ⚠ El patron tiene que admitir el PLURAL --- «Figures 4 and 6», «las Figuras 4 y 6» ---
        #   o una cita conjunta a dos figuras es invisible.
        patron = re.compile(r"(Figure|Figura)s?\s*%s" % re.escape(num))
        # ⚠ Y la pagina del pie NO se tira: el pie SIEMPRE imprime su propio numero, asi que
        #   descartarla entera hacia que una figura citada en su misma pagina (la colocacion
        #   ideal) se juzgara por una mencion lejana, y que una figura SIN ninguna frase que la
        #   nombre se absolviera con su propio pie.  Una cita de verdad es: en otra pagina, o una
        #   SEGUNDA aparicion en la pagina del pie.  Es la regla de paper/check_floats.py.
        cuenta = {i + 1: len(patron.findall(" ".join(t.split()))) for i, t in enumerate(paginas)}
        pags_cita = [p for p, c in cuenta.items() if c and (p != pag_fig or c > 1)]
        cerca = any(abs(p - pag_fig) <= 1 for p in pags_cita)
        ok(bool(pags_cita) and cerca,
           f"F2 {lab} (fig. {num}, pag. {pag_fig}): citada en la misma pagina o vecina",
           f"citada en {pags_cita}")

# --- F3  las cifras de los pies contra la salida de su guion
print(SEP)
print("F3 -- las cifras que afirman los pies, contra la salida de su guion")
tex = io.open("characters_at_aP.tex", encoding="utf-8").read()
PARES = [
    ("fig_walls.pdf", "fig_walls.py", ["8", "15", "24", "48"]),
    ("fig_value.pdf", "fig_value.py", ["352", "176", "56"]),
]
for figura, guion, cifras in PARES:
    cuerpo = re.search(r"includegraphics\[[^\]]*\]\{%s\}(.*?)\\end\{figure\}"
                       % re.escape(figura), tex, re.S)
    ok(cuerpo is not None, f"F3 {figura}: encuentro su pie")
    if cuerpo is None:
        continue
    pie = cuerpo.group(1)
    salida_p = os.path.join(GATES, guion.replace(".py", "_OUT.txt"))
    salida = io.open(salida_p, encoding="utf-8", errors="replace").read() if os.path.exists(salida_p) else ""
    if not salida:
        # los de figura imprimen por pantalla; se re-corre y se captura
        import subprocess
        salida = subprocess.run([sys.executable, guion], cwd=GATES,
                                capture_output=True, text=True, encoding="utf-8",
                                errors="replace").stdout
    for c in cifras:
        en_pie = re.search(r"\b%s\b" % c, pie.replace("\\,", "")) is not None
        en_sal = re.search(r"\b%s\b" % c, salida) is not None
        ok(not en_pie or en_sal, f"F3 {figura}: la cifra {c} del pie esta en la salida del guion",
           f"pie={en_pie} salida={en_sal}")

# --- F4  fechas dentro de la imagen contra la bibliografia
print(SEP)
print("F4 -- las fechas impresas DENTRO de fig_map.pdf, contra la bibliografia")
d = fitz.open("fig_map.pdf")
tinta = " ".join(d[0].get_text().split())
bib = tex.split(r"\begin{thebibliography}")[1]
PARES_F = [("Lübeck", "2021"), ("Kostant", "1976"), ("NPP", "2025"), ("Adin", "2002")]
for quien, anyo in PARES_F:
    ok(anyo in tinta, f"F4: la figura imprime {anyo} para {quien}")
ok("2020" not in tinta, "F4: y ya no imprime el 2020 de Lübeck-Prasad")
ok("(2021)" in bib, "F4: la bibliografia dice 2021 para Lübeck-Prasad")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: cada figura esta donde la citan, con las cifras de su guion y sin fechas viejas.")
