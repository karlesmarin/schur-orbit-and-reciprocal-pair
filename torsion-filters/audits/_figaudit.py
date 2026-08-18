# -*- coding: utf-8 -*-
# Auditoria de figuras: cuales tienen \label y cuales se CITAN en el texto, y el hueco de pagina.
# Un pie no es una cita: una figura que nadie referencia no la mira nadie.
# Authors: Carles Marin, Claude (AI assistant).
import io
import re
import subprocess
import sys

s = io.open('orbit_pair_ii.tex', encoding='utf-8').read()
cuerpo = s.split(r'\begin{thebibliography}')[0]
labels = re.findall(r'\\label\{(fig:[^}]*)\}', s)
refs = set(re.findall(r'\\ref\{(fig:[^}]*)\}', cuerpo))
incl = re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]*)\}', s)
print("figuras incluidas (%d): %s" % (len(incl), incl))
print("con \\label   (%d): %s" % (len(labels), labels))
print("CITADAS      (%d): %s" % (len(refs), sorted(refs)))
sin = [l for l in labels if l not in refs]
print("SIN CITAR    (%d): %s" % (len(sin), sin))
if sin:
    sys.exit(1)
print("todas las figuras se citan en el texto.")
