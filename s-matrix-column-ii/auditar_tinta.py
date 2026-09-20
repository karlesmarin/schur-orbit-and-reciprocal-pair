# -*- coding: utf-8 -*-
r"""auditar_tinta.py -- lo que esta ESCRITO DENTRO de una figura tambien es una afirmacion.

auditar_figuras.py compara cada figura con su PIE y con la salida de su compuerta.  No mira el
texto que va dibujado dentro del PDF de la figura, y eso es un agujero: la figura 2 llego a
rotular la ley como la seccion 5, a Gauss como la 4 y a Bring como la 3, que era el ORDEN
ANTERIOR.  El pie ya usa el orden nuevo.  Corregir la leyenda no arregla nada:
hay que regenerar la imagen.

Un numero de seccion dentro de una imagen es tinta: no se renumera solo, ningun \ref lo alcanza,
y ninguna de las otras auditorias lo lee.  Aqui se lee.

    python auditar_tinta.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import io
import os
import re
import sys

import fitz

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))

fallos = []
hechas = []


def afirma(ok, que):
    print("   %-4s %s" % ("ok" if ok else "MAL", que))
    hechas.append(que)
    if not ok:
        fallos.append(que)


# el numero que LaTeX imprime para cada seccion, desde el .aux, que es la unica fuente de verdad
numero = {}
for m in re.finditer(r"\\newlabel\{(sec:\w+)\}\{\{(\d+)\}",
                     io.open(os.path.join(AQUI, "P3b.aux")).read()):
    numero[m.group(1)] = m.group(2)

print("secciones segun el .aux: %s" % numero)

# Prohibir el numero seria peor que el fallo: la figura 2 lo necesita, porque reparte el eje del
# defecto entre las secciones que lo explican.  Lo que hay que exigir es que COINCIDA.  Cada
# numero dibujado se identifica por la sena que lleva al lado, declarada aqui a proposito: la
# maquina no puede adivinar a que seccion se refiere un "5" suelto.
ESPERADO = {
    "fig_mapa": [("2/p", "sec:ley"), ("Gau", "sec:centro"), ("Bring", "sec:bring")],
    "fig_mapa_es": [("2/p", "sec:ley"), ("semisistemas", "sec:centro"), ("Bring", "sec:bring")],
}
patron = re.compile(r"(?:\u00a7|\bS\b|\bSec(?:tion|ci\u00f3n)?\.?)\s*(\d)[^\w]{0,3}([^\u00a7]{0,26})")

print("")
print("CADA NUMERO DE SECCION DIBUJADO COINCIDE CON EL QUE LA NOTA IMPRIME")
for nombre in sorted(f for f in os.listdir(AQUI) if f.startswith("fig_") and f.endswith(".pdf")):
    base = nombre[:-4]
    texto = " ".join(" ".join(p.get_text().split())
                     for p in fitz.open(os.path.join(AQUI, nombre)))
    marcas = patron.findall(texto)
    if base not in ESPERADO:
        afirma(not marcas, "%-22s no dibuja ningun numero de seccion%s"
               % (nombre, "" if not marcas else " --- pero lleva %s" % marcas))
        continue
    afirma(len(marcas) == len(ESPERADO[base]),
           "%-22s dibuja %d numeros, como se espera" % (nombre, len(marcas)))
    for (dibujado, cola), (sena, etiqueta) in zip(marcas, ESPERADO[base]):
        afirma(sena in cola and dibujado == numero[etiqueta],
               "%-22s %-13s -> dibuja %s, la nota imprime %s"
               % (nombre, sena, dibujado, numero[etiqueta]))

print("")
print("CONTROL: el criterio reconoce un numero de seccion plantado")
for prueba in ("la ley de S 5x", "ver Section 3x", "Secci\u00f3n 4x"):
    afirma(bool(patron.search(prueba)), "reconoce %r" % prueba)

print("")
print("=" * 82)
print("comprobaciones: %d ; fallos: %d" % (len(hechas), len(fallos)))
if not hechas:
    sys.exit("FATAL: cero comprobaciones. No poder medir no es aprobar")
print("VEREDICTO:", "ninguna figura lleva tinta que el texto pueda desmentir" if not fallos
      else "TINTA CADUCADA:\n   - " + "\n   - ".join(fallos))
sys.exit(1 if fallos else 0)
