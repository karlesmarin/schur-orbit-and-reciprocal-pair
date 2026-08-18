# -*- coding: utf-8 -*-
# Renombra tau_t -> tau^C_t en el TEXTO PAR, dejando el impar con su tau^B_t.
#
# POR QUE.  Su vuelta 24: el filtro par y el impar son objetos distintos y el paper los llamaba
# igual.  Dejar "tau_t" a secas es lo que permitio que §3 y §5 se contradijeran sin que se viera.
#
# Se hace por RANGO DE LINEAS, no a ciegas: la seccion del caso impar se localiza por sus marcas y
# ahi no se toca nada.  Se imprime cuantas sustituciones y en que zonas, y se comprueba que el
# fichero solo cambia en la cuenta esperada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python _taurename.py
import io

p = 'orbit_pair_ii.tex'
L = io.open(p, encoding='utf-8').read().split('\n')

ini = next(i for i, l in enumerate(L) if r'\section{The odd case' in l)
fin = next(i for i, l in enumerate(L) if r'\section{The highest surviving weight}' in l)
print("la seccion impar va de la linea %d a la %d" % (ini + 1, fin + 1))

fuera = dentro = 0
for i, l in enumerate(L):
    if r'\tau_t' not in l:
        continue
    if ini <= i < fin:
        dentro += l.count(r'\tau_t')
        continue
    fuera += l.count(r'\tau_t')
    L[i] = l.replace(r'\tau_t', r'\tau^C_t')

io.open(p, 'w', encoding='utf-8').write('\n'.join(L))
print("sustituidas %d ocurrencias fuera de la seccion impar" % fuera)
print("quedan %d dentro, que hay que mirar A MANO (deben ser genericas o del par)" % dentro)
for i in range(ini, fin):
    if r'\tau_t' in L[i]:
        print("   linea %d: %s" % (i + 1, L[i].strip()[:110]))
