# -*- coding: utf-8 -*-
# Mueve la subseccion del Lema de regularidad de la seccion del caso impar a la del filtro,
# que es donde pertenece: es el enunciado invariante del filtro y vale para las dos paridades.
# La Observacion B/C se queda en el impar, porque esa SI es sobre la paridad.
#
# Verifica que el bloque se mueve entero y que el fichero no pierde ni gana nada mas.
# Authors: Carles Marin, Claude (AI assistant).
import io

p = 'orbit_pair_ii.tex'
s = io.open(p, encoding='utf-8').read()
n0 = len(s)

ini = s.index(r'\subsection{The filter is regularity in the group itself}')
fin = s.index(r'\begin{remark}[the parity dichotomy, second switch]')
bloque = s[ini:fin]
assert r'\label{lem:regular}' in bloque, "el bloque no lleva el lema"
resto = s[:ini] + s[fin:]

# destino: justo antes de la Verification de la seccion del filtro (la primera que aparece)
ancla = r'\subsection{Verification}' + '\n' + \
        r'Lemma~\ref{lem:T} was checked against two independent computations'
assert ancla in resto, "no encuentro el ancla de la Verification del filtro"
resto = resto.replace(ancla, bloque + ancla, 1)

assert len(resto) == n0, "el fichero cambio de tamaño: %d -> %d" % (n0, len(resto))
io.open(p, 'w', encoding='utf-8').write(resto)
print("movido: %d caracteres, el fichero conserva el tamaño" % len(bloque))
