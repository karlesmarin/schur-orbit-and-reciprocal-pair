# -*- coding: utf-8 -*-
# Auditoria de citas: entradas de la bibliografia que nunca se citan, y citas sin entrada.
# Authors: Carles Marin, Claude (AI assistant).
import io
import re

s = io.open('orbit_pair_ii.tex', encoding='utf-8').read()
cuerpo = s.split(r'\begin{thebibliography}')[0]
listadas = set(re.findall(r'\\bibitem\[[^\]]*\]\{([^}]*)\}', s))
citadas = set()
# El argumento opcional de \cite[...] {...} lo tenia ciego: el 20 de agosto declaro Hall48 «nunca
# citada» estando citada dos veces como \cite[Theorem 2]{Hall48}.  Un control que no ve una forma
# valida del comando no es estricto, es incompleto --- y su falso positivo cuesta lo mismo que un
# falso negativo, porque invita a «arreglar» lo que no esta roto.
for m in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]*)\}', cuerpo):
    for k in m.split(','):
        citadas.add(k.strip())
print("entradas en la bibliografia :", len(listadas))
print("claves citadas en el cuerpo :", len(citadas))
print("NUNCA citadas               :", sorted(listadas - citadas))
print("citadas y SIN entrada       :", sorted(citadas - listadas))
