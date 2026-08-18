# -*- coding: utf-8 -*-
# Auditoria de los cinco puntos de su vuelta 27, contra el .tex.
#
# Sus cinco, tal como los enuncio:
#   1  Remark 5.6 citaba "Lemma 3.1" para el menor congelado, que ahora es el resultado PAR.
#   2  las referencias a la localizacion decian 5.4; con el nuevo 5.5 hay que apuntar a sec:oddlocal.
#   3  lenguaje residual de 'unit': el titulo de la seccion 8, la Observacion 5.7 y el Problema 17.5.
#   4  la frase de las cintas en el cierre era MATEMATICAMENTE FALSA: el signo de Littlewood es el
#      del teselado (alturas), no la cantidad de cintas; y el checkerboard no usa el signo, usa
#      |lambda|-|nu| = tk y |nu|-|mu| = 2j.
#   5  "las dos paredes del Lema 3.1 se funden en el lado impar": el Lema 3.1 ya no tiene lado impar.
#
# Y ademas, de su propuesta nueva: que GKRS quede citado en la zona de la localizacion.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python _vuelta27.py
import io

s = io.open('orbit_pair_ii.tex', encoding='utf-8').read()

CHECKS = [
    ("1 · el menor congelado ya no invoca el Lema 3.1",
     lambda: 'by the type-$B$ Weyl-numerator criterion underlying Proposition' in s
             and r'the frozen minor is a $B_{m'"'"'}$ Weyl numerator at the torsion' in s),
    ("1b · y dice CUAL es el criterio, no solo su nombre",
     lambda: r'two columns coincide as soon as $A_i\equiv\pm A_j$' in s),
    ("2 · la localizacion se cita como sec:oddlocal, no sec:oddverif",
     lambda: s.count(r'\ref{sec:oddlocal}') >= 3),
    ("2b · y §5.1 manda a las dos: medida y localizada",
     lambda: r'\ref{sec:oddverif}, and localised in \S\ref{sec:oddlocal}' in s),
    ("3 · la Observacion 5.7 habla de la MITAD IMPAR de la conjetura",
     lambda: 'for the odd half of Conjecture' in s
             and 'for the odd analogue rather than for the conjecture itself' not in s),
    ("3b · la seccion 8 ya no se llama 'The unit conjecture'",
     lambda: r'\section{The extremal primitivity conjecture}' in s
             and r'\section{The unit conjecture}' not in s),
    ("3c · y el Problema 17.5 tampoco dice 'the unit is the part'",
     lambda: 'the unit is the part of the structure' not in s
             and 'extremal primitivity is the part of the structure' in s),
    ("4 · el cierre ya NO dice que el signo sea el numero de cintas",
     lambda: 'the parity of the number of $t$-ribbons removed' not in s),
    ("4b · y dice lo correcto: tamaños, kt celdas y un par",
     lambda: r'differs from $\lambda$ by $kt$ cells' in s
             and 'changes the size by an even number' in s),
    ("4c · y avisa explicitamente de que no interviene ningun signo",
     lambda: 'no \\emph{sign} enters here' in s),
    ("5 · el cierre no dice que el Lema 3.1 tenga lado impar",
     lambda: 'The two walls of\nLemma \\ref{lem:T} collapse to one on the odd side' not in s
             and r'$a_i\equiv0$ and' in s and r'$2a_i\equiv0$ are different at an even modulus' in s),
    ("+ · GKRS citado en la zona de la localizacion",
     lambda: r'\cite{GKRS}' in s and r'\bibitem[GKRS98]{GKRS}' in s),
    ("+ · y Landweber, para el multiplete de Dirac",
     lambda: r'\cite{Landweber}' in s and r'\bibitem[Lan01]{Landweber}' in s),
    ("+ · el factor cruzado enunciado y probado como denominador relativo",
     lambda: r'\begin{proposition}[the cross factor is the relative denominator of an equal-rank pair]' in s),
    ("+ · Springer citado para la conjugacion de los elementos regulares",
     lambda: r'\cite{Springer74}' in s and r'\bibitem[Spr74]{Springer74}' in s),
    ("+ · y la sugerencia de GKRS atribuida a un correspondiente",
     lambda: 'is due to a\ncorrespondent, and we are grateful for it' in s),
]

print("=" * 84)
print("SUS PETICIONES DE LA VUELTA 27, CONTRA EL FICHERO")
print("=" * 84)
hechos = 0
for nombre, f in CHECKS:
    try:
        ok = bool(f())
    except Exception:
        ok = False
    hechos += ok
    print("  [%s]  %s" % ("HECHO    " if ok else "PENDIENTE", nombre))
print("=" * 84)
print("  %d de %d" % (hechos, len(CHECKS)))
