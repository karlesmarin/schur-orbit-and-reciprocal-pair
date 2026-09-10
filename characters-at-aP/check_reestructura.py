# -*- coding: utf-8 -*-
# AUDITORIA DE LA REESTRUCTURACION.   8 de septiembre de 2026.
#
# POR QUE.  El articulo se ha reordenado entero, dos veces y en dos lenguas, moviendo bloques
# verbatim.  Un movimiento de bloques no rompe la compilacion pero SI deja: frases huerfanas que
# remiten a un sitio que ya no esta antes, afirmaciones viejas que el resultado nuevo desmiente,
# parrafos duplicados, figuras cuyo pie ya no describe lo que hay, y referencias colgando.
# Nada de eso lo ve pdflatex.
#
# QUE COMPRUEBA
#   R1  CITAS: toda \cite tiene \bibitem y todo \bibitem se cita  (EN y ES)
#   R2  FIGURAS: cada \includegraphics tiene \caption y \label, y cada \label se referencia
#   R3  DUPLICADOS: parrafos largos repetidos --- la firma de un bloque pegado dos veces
#   R4  HUERFANAS: remisiones hacia atras a secciones que ahora van despues
#   R5  OBSOLETAS: frases que declaran abierto lo que el articulo ya cierra
#   R6  PREGUNTAS: ninguna de las abiertas puede estar contestada dentro del propio articulo
#   R7  HUECOS: paginas del PDF con poca tinta (floats mal colocados)
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_reestructura.py > check_reestructura_OUT.txt 2>&1

import io
import re
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 96
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


TEX = {n: io.open("%s.tex" % n, encoding="utf-8").read() for n in ("characters_at_aP", "characters_at_aP_es")}

# ---------------------------------------------------------------- R1
print(SEP)
print("R1 -- CITAS: toda \\cite tiene entrada, y toda entrada se cita")
for n, t in TEX.items():
    cuerpo = t.split(r"\begin{thebibliography}")[0]
    citadas = {k.strip() for g in re.findall(r"\\cite\{([^}]+)\}", cuerpo) for k in g.split(",")}
    entradas = set(re.findall(r"\\bibitem\{([^}]+)\}", t))
    ok(not (citadas - entradas), f"{n}: ninguna cita sin entrada", str(sorted(citadas - entradas)))
    ok(not (entradas - citadas), f"{n}: ninguna entrada sin citar", str(sorted(entradas - citadas)))

# ---------------------------------------------------------------- R2
print(SEP)
print("R2 -- FIGURAS: pie, etiqueta y al menos una referencia en el texto")
for n, t in TEX.items():
    figs = re.findall(r"\\begin\{figure\}(.*?)\\end\{figure\}", t, re.S)
    print(f"   {n}: {len(figs)} figuras")
    for f in figs:
        img = re.search(r"includegraphics\[[^\]]*\]\{([^}]+)\}", f)
        lab = re.search(r"\\label\{([^}]+)\}", f)
        cap = "\\caption" in f
        nombre = img.group(1) if img else "??"
        if not (img and lab and cap):
            ok(False, f"{n}: figura {nombre} sin pie o sin etiqueta")
            continue
        refs = len(re.findall(r"\\ref\{%s\}" % re.escape(lab.group(1)), t))
        ok(refs >= 1, f"{n}: {nombre} ({lab.group(1)}) referenciada en el texto",
           f"{refs} referencias")

# ---------------------------------------------------------------- R3
print(SEP)
print("R3 -- DUPLICADOS: parrafos largos repetidos (firma de un bloque pegado dos veces)")
for n, t in TEX.items():
    parr = [" ".join(p.split()) for p in t.split("\n\n")]
    parr = [p for p in parr if len(p) > 220 and not p.startswith("%")]
    c = Counter(parr)
    rep = [p for p, k in c.items() if k > 1]
    ok(not rep, f"{n}: ningun parrafo largo duplicado", f"{len(rep)} duplicados")
    for p in rep[:3]:
        print("      ...", p[:110])

# ---------------------------------------------------------------- R4
print(SEP)
print("R4 -- HUERFANAS: remisiones que apuntan a una seccion posterior con 'above/arriba'")
ORDEN = ["sec:map", "sec:norm", "sec:crit", "sec:alph", "sec:cap", "sec:valor", "sec:dims",
         "sec:cuenta", "sec:g2", "sec:types", "sec:open", "sec:cierre", "sec:preguntas",
         "sec:controles"]
for n, t in TEX.items():
    pos = {}
    for s in ORDEN:
        m = re.search(r"\\label\{%s\}" % re.escape(s), t)
        pos[s] = m.start() if m else None
    malas = []
    for m in re.finditer(r"(above|arriba|earlier|antes)[^.]{0,60}?\\S?\\ref\{(sec:[a-z0-9]+)\}", t):
        s = m.group(2)
        if pos.get(s) and m.start() < pos[s]:
            malas.append((m.group(0)[:60], s))
    ok(not malas, f"{n}: ninguna remision 'arriba' a una seccion posterior", str(malas[:3]))

# ---------------------------------------------------------------- R5
print(SEP)
print("R5 -- OBSOLETAS: el articulo no puede seguir declarando abierto lo que ya cierra")
PROHIBIDAS_EN = [
    ("the general value is open", "el valor general ya no esta abierto"),
    ("even the absolute value is open", "el valor absoluto en d>=3 ya no esta abierto"),
    ("its sign, open", "el signo ya no esta abierto: Teorema del valor con signo"),
    ("what its sign is, I do not know", "idem"),
    ("the sign I have not determined", "idem"),
]
PROHIBIDAS_ES = [
    ("el valor general está abierto", "el valor general ya no esta abierto"),
    ("está abierto hasta el valor absoluto", "idem"),
    ("el signo, abierto", "el signo ya no esta abierto"),
    ("cuál es su signo, no lo sé", "idem"),
    ("su signo no lo he determinado", "idem"),
]
# 8-sep: R5 tenia «the sign I have not determined» en su lista y NO lo vio, porque en el fichero la
#   frase esta partida por un salto de linea y aqui se buscaba con un espacio simple.  Una compuerta
#   sobre LaTeX que no normaliza los blancos es ciega a cualquier frase de dos lineas, o sea a casi
#   todas.  Se normaliza antes de buscar.
for n, lista in (("characters_at_aP", PROHIBIDAS_EN), ("characters_at_aP_es", PROHIBIDAS_ES)):
    t = re.sub(r"\s+", " ", TEX[n]).lower()
    for frase, motivo in lista:
        ok(frase.lower() not in t, f"{n}: no dice «{frase[:42]}»", motivo)

# ---------------------------------------------------------------- R6
print(SEP)
print("R6 -- PREGUNTAS: las abiertas no pueden estar contestadas en el propio articulo")
for n, t in TEX.items():
    m = re.search(r"\\label\{sec:preguntas\}(.*?)\\appendix", t, re.S)
    if not m:
        ok(False, f"{n}: encuentro la seccion de preguntas"); continue
    bloque = m.group(1)
    items = re.findall(r"\\item \\textbf\{([^}]+)\}", bloque)
    # ⚠ 9-sep-2026: eran cinco.  La clasificacion del Teorema thm:clasif abrio una sexta ---
    #   cuanto vale en las tres familias que NO son las potencias --- y el titulo de la seccion
    #   dejo de decir un numero, para que no vuelva a quedarse rancio al anyadir o quitar una.
    #   Lo que la compuerta cuida no es el numero sino que las DOS ediciones tengan el mismo.
    ok(len(items) == 6, f"{n}: seis preguntas", f"{len(items)}")
    for it in items:
        print("      -", " ".join(it.split())[:88])
    # ninguna puede citar un teorema del propio articulo como respuesta ya dada
    ok("Theorem~\\ref{thm:signo} settles" not in bloque
       and "queda resuelto por el Teorema" not in bloque,
       f"{n}: ninguna pregunta se contesta a si misma")

# ---------------------------------------------------------------- R7
print(SEP)
print("R7 -- HUECOS: paginas con poca tinta")
try:
    import fitz
    for n in ("characters_at_aP", "characters_at_aP_es"):
        d = fitz.open("%s.pdf" % n)
        flojas = []
        for i, p in enumerate(d):
            txt = len(p.get_text().strip())
            imgs = len(p.get_images())
            if txt < 900 and imgs == 0 and i < d.page_count - 1:
                flojas.append((i + 1, txt))
        print(f"   {n}: {d.page_count} paginas; con menos de 900 caracteres y sin figura: {flojas}")
except Exception as e:
    print("   (sin PyMuPDF: %s)" % e)

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: la reestructuracion no ha dejado huerfanas, duplicados ni obsoletas.")
