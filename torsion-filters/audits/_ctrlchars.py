# -*- coding: utf-8 -*-
# NI UN CARACTER DE CONTROL EN EL FUENTE, Y NADA QUE LA EXTRACCION SE COMA EN UN ENUNCIADO.
# 18 de agosto de 2026.
#
# POR QUE EXISTE ESTO.  Dos fallos reales del mismo dia, ninguno de los cuales veia ninguna de las
# veinticinco auditorias que ya habia:
#
#   (1) Un heredoc de bash interpreto `\b` como RETROCESO y escribio 0x08 dentro de `$M(\beta)$` en
#       la edicion castellana.  LaTeX compilo sin una queja --- el byte es invisible --- y el PDF
#       salio con la macro comida.  Un fuente con caracteres de control es un fuente roto que no
#       avisa.
#   (2) Una lectura externa vio conj:H terminando en "=1" cuando pone "|...|=1": el extractor de
#       texto convierte el glifo de la barra grande en 0x0c y se lo come.  El PDF esta bien; lo que falla es leerlo por extraccion.  137 delimitadores en
#       26 de las 69 paginas.
#
# LO QUE COMPRUEBA
#   C1  FATAL  cero caracteres de control en los dos .tex.  Ni uno.
#   C2         los delimitadores grandes que caen DENTRO de un enunciado numerado --- theorem,
#              proposition, lemma, corollary, conjecture --- y que por tanto un lector-por-extraccion
#              puede perder.  No es fatal: es una lista para decidir si ese enunciado deberia
#              escribirse de forma que no dependa del glifo, como se hizo con conj:H.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _ctrlchars.py
import io
import re
import sys
import unicodedata

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FUENTES = ["orbit_pair_ii.tex", "orbit_pair_ii_es.tex"]
ENUNCIADOS = ("theorem", "proposition", "lemma", "corollary", "conjecture")
GRANDES = re.compile(r"\\(?:B|b)ig[lr]?\||\\left\||\\right\||\\(?:B|b)ig[lr]\\?[{(]|\\left[{(]")

print("=" * 96)
print("CARACTERES DE CONTROL EN EL FUENTE, Y DELIMITADORES DE RIESGO EN LOS ENUNCIADOS")
print("=" * 96)

fatal = 0
for f in FUENTES:
    s = io.open(f, encoding="utf-8").read()
    malos = []
    for i, ch in enumerate(s):
        if ch in "\n\t":
            continue
        if unicodedata.category(ch) in ("Cc", "Cf"):
            malos.append((s.count("\n", 0, i) + 1, hex(ord(ch)),
                          " ".join(s[max(0, i - 45):i + 45].split())))
    print("")
    print("  C1  %-22s caracteres de control: %d" % (f, len(malos)))
    for ln, c, ctx in malos:
        fatal += 1
        print("        !! linea %d  %s  ...%s..." % (ln, c, ctx))

print("")
print("  C2  delimitadores grandes DENTRO de enunciados numerados (riesgo de extraccion):")
for f in FUENTES:
    s = io.open(f, encoding="utf-8").read()
    n = 0
    for env in ENUNCIADOS:
        for m in re.finditer(r"\\begin\{%s\}(.*?)\\end\{%s\}" % (env, env), s, re.S):
            cuerpo = m.group(1)
            hallados = GRANDES.findall(cuerpo)
            if not hallados:
                continue
            etq = re.search(r"\\label\{([^}]*)\}", cuerpo)
            n += len(hallados)
            print("        %-22s %-18s %d  (%s)"
                  % (f, etq.group(1) if etq else env, len(hallados),
                     ", ".join(sorted(set(hallados))[:4])))
    print("        %-22s total en enunciados: %d" % (f, n))

print("")
print("  VEREDICTO  C1 (fatal): %s" % ("PASA" if fatal == 0 else "FALLA, %d caracteres" % fatal))
print("             C2 es una lista para mirar, no un fallo.")
if fatal:
    raise SystemExit("*** hay caracteres de control en el fuente ***")
print("")
print("=" * 96)
print("DONE")
