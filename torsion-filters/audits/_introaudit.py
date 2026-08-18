# -*- coding: utf-8 -*-
# Auditoria de la INTRODUCCION contra el cuerpo del paper.
#
# POR QUE.  Hemos ido escribiendo encima de lo que habia, y la introduccion y el abstract son
# justo lo que no se relee al insertar una seccion.  Esto busca frases que el propio cuerpo ya
# contradice, y comprueba que la introduccion mencione cada seccion que existe.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python _introaudit.py
import io
import re

s = io.open('orbit_pair_ii.tex', encoding='utf-8').read()
ini = s.index(r'\section{Introduction}')
fin = s.index(r'\section{\texorpdfstring{The reduction')
intro = s[ini:fin]
# La rebanada empieza DESPUES del marcador: cogiendolo dentro, los 16 caracteres de
# `\begin{abstract}` entraban en la cuenta y el auditor anunciaba 1921 --- "se pasa del tope de
# arXiv" --- cuando el texto que se pega mide 1908.  Un falso positivo sobre el numero que decide
# si el envio se rechaza es peor que no medirlo: la cuenta buena la da `_abslen.py`, que ademas
# expande `\emph` y las macros nuestras, y esta calibrado contra la v2 del Paper I que entro.
ab = s[s.index(r'\begin{abstract}') + len(r'\begin{abstract}'):s.index(r'\end{abstract}')]


def solo_prosa(txt):
    """quita el contenido de \\ref, \\label y \\cite: la palabra 'unit' dentro de \\ref{sec:unit}
       no es prosa desfasada, y una primera version de este auditor la conto como tal."""
    # y se colapsan los saltos de linea: si no, una frase buscada a mano no casa cuando el .tex la
    # parte en dos lineas, y el auditor da un falso positivo que cuesta media hora localizar.
    return " ".join(re.sub(r"\\(ref|eqref|label|cite)\{[^}]*\}", " ", txt).split())


intro_p, ab_p = solo_prosa(intro), solo_prosa(ab)

SOSPECHOSAS = [
    ("afirma branching simplectico para todo t", "symplectic branching"),
    ("dice 'for every t'", r"for every \$t\$"),
    # 'unit' solo es error cuando nombra NUESTRA clase extremal; la constante de NPP si es una
    # unidad de verdad ("collapses to a unit only when that group is a torus") y esa frase es buena.
    ("dice 'unit' de nuestra clase (deberia ser primitive)",
     r"\bunit\b(?!\s+only when that group is a torus)"),
    ("usa eta+rho sin doblar en tipo B", r"\(\\eta\+\\rho\)\(\\xi\)"),
    ("atribuye el shift al rho del factor libre", "positive roots of the free"),
    ("atribuye el 0,pm1 a Kostant", r"\{0,\\pm1\}.{0,80}Kostant"),
]
print("=" * 78)
print("AUDITORIA DE ABSTRACT + INTRODUCCION")
print("=" * 78)
mal = 0
for zona, txt in (("abstract", ab_p), ("introduccion", intro_p)):
    print("\n  %s (%d caracteres)" % (zona, len(txt)))
    for nombre, pat in SOSPECHOSAS:
        n = len(re.findall(pat, txt))
        if n:
            mal += 1
            print("     !! %-42s : %d" % (nombre, n))
    if not any(re.findall(p, txt) for _, p in SOSPECHOSAS):
        print("     sin frases desfasadas")
    # La longitud que decide el envio NO se mide aqui: `\zt` cuenta como cuatro caracteres y
    # `\emph{...}` hay que quitarlo antes de pegar, y esas dos reglas viven en `_abslen.py`.  Este
    # numero es solo el tamano de la prosa con las referencias fuera, para ver si crece.
    if zona == "abstract":
        n = len(" ".join(txt.split()))
        print("     %-42s : %d   (el tope de arXiv lo mide _abslen.py)"
              % ("prosa del abstract, sin \\ref/\\cite", n))

# toda seccion del cuerpo tiene que aparecer en el recorrido de la introduccion
secs = re.findall(r"\\section\{[^}]*\}\\label\{(sec:[a-z0-9]+)\}", s)
secs += re.findall(r"\\section\{\\texorpdfstring\{[^}]*\}\{[^}]*\}\}\\label\{(sec:[a-z0-9]+)\}", s)
citadas = set(re.findall(r"\\ref\{(sec:[a-z0-9]+)\}", intro))
faltan = [x for x in secs if x not in citadas]
print("\n  secciones del cuerpo: %d | citadas en la introduccion: %d" % (len(secs), len(citadas)))
print("  NO mencionadas en el recorrido: %s" % (faltan if faltan else "ninguna"))
print("=" * 78)
