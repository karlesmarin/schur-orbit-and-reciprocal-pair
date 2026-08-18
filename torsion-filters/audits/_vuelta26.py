# -*- coding: utf-8 -*-
# Auditoria de los puntos de su vuelta 26 contra el .tex, uno por uno.
# OJO: buscar la palabra 'odd' a secas da falso positivo -- el lema dice "odd in a_j",
# que es funcion impar.  Se buscan las frases, no la palabra.
# Cada linea es (que pidio, patron que TIENE que estar, patron que NO puede estar).
# Authors: Carles Marin, Claude (AI assistant).
import io
import re

s = io.open('orbit_pair_ii.tex', encoding='utf-8').read()
lem = s[s.index(r'\begin{lemma}[the filter, with its sign]'):]
lem = lem[:lem.index(r'\end{proof}')]

CHECKS = [
    ("1 · Lema 3.1 sin restos de impar",
     lambda: (r'Let $t=2m+2$ be even' in lem) and ('odd $t$' not in lem) and ('For odd' not in lem) and ("2m'+1" not in lem)),
    ("1b · §3.2 cita la Prop 5.3 y no 'el caso impar'",
     lambda: r'Lemma \ref{lem:T} lists three species of wall, and Proposition \ref{prop:oddfilter}' in s),
    ("2 · pie de la Figura 2 distingue las dos reglas",
     lambda: 'Proposition \\ref{prop:oddfilter}, odd order of the same rank' in s
             and 'The two filters drawn rather than described' in s),
    ("3 · fuera la contradiccion affine/finite de la Obs 3.3",
     lambda: 'affine wall of the level rather than a root' not in s),
    ("4 · NPP no presentado como cota universal",
     lambda: 'which is more general than a bound' in s),
    ("4b · la constante en tipo B, atribuida bien",
     lambda: 'it is already $1$ in\n\\cite{NPP}' in s or 'already $1$ in' in s),
    ("5 · el Corolario 3.5 cubre las dos paridades",
     lambda: 'when $t$ is odd the same holds in the coordinates' in s),
    ("6 · la Prop 11.1 es teorema con 'Let t be even'",
     lambda: r'Let $t$ be even. Then \eqref{eq:square} holds for every $\mu$.' in s),
    ("7 · PROBAR la formula local (13)",
     lambda: r'\begin{lemma}[the parity is local]' in s
             and r'Lemma \ref{lem:local}: the local form of the $r_i$ & box $\le8$, $R=3$ &' in s),
    ("8 · la tabla lleva el Teorema 6.1",
     lambda: r'Theorem~\ref{thm:fusion}: both filters are the minimal fusion projections' in s),
    ("8b · el desplazamiento ya no es 'verified'",
     lambda: r'Proposition~\ref{prop:newtden}: $\mathrm{top}\,\Newt(N_\delta)$, hence the shift & proved' in s),
    ("8c · la etiqueta de la Figura 5",
     lambda: 'the sum of the positive roots of the free factor' not in
             io.open('fig_law3d.py', encoding='utf-8').read()),
    ("9 · el Problema 17.1 reconoce lo que §6 ya hace",
     lambda: 'Two of the three things that would have to hold are now available' in s),
    ("10 · 'by triangularity' explicitado",
     lambda: 'by\nthe usual triangularity of a product of characters' in s
             or 'usual triangularity of a product of characters' in s),
    ("11 · la fila impar de verificacion cita la Prop 5.3",
     lambda: r'Proposition \ref{prop:oddfilter}, odd $t$, both directions' in s),
    ("12 · la conclusion habla de las TRES manifestaciones",
     lambda: r'\section{Closing: one parity, seen three times}' in s
             and 'third \\emph{manifestation} rather than a third independent mechanism' in s),
    ("13 · el esquema branching -> fusion como enunciado estructural",
     lambda: 'minimal fusion projection on the frozen factor' in s),
]

print("=" * 84)
print("SUS PETICIONES DE LA VUELTA 26, CONTRA EL FICHERO")
print("=" * 84)
hechos = 0
for nombre, f in CHECKS:
    try:
        ok = bool(f())
    except Exception as e:
        ok = False
    hechos += ok
    print("  [%s]  %s" % ("HECHO    " if ok else "PENDIENTE", nombre))
print("=" * 84)
print("  %d de %d" % (hechos, len(CHECKS)))
