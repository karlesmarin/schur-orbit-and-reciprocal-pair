# -*- coding: utf-8 -*-
# Auditoria de citas: entradas de la bibliografia que nunca se citan, y citas sin entrada.
# Authors: Carles Marin, Claude (AI assistant).
import io
import re

s = io.open('orbit_pair_ii.tex', encoding='utf-8').read()
cuerpo = s.split(r'\begin{thebibliography}')[0]
listadas = set(re.findall(r'\\bibitem\[[^\]]*\]\{([^}]*)\}', s))
citadas = set()
for m in re.findall(r'\\cite\{([^}]*)\}', cuerpo):
    for k in m.split(','):
        citadas.add(k.strip())
print("entradas en la bibliografia :", len(listadas))
print("claves citadas en el cuerpo :", len(citadas))
print("NUNCA citadas               :", sorted(listadas - citadas))
print("citadas y SIN entrada       :", sorted(citadas - listadas))
