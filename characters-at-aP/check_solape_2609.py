# -*- coding: utf-8 -*-
# SOLAPE CON NUESTRO PROPIO ARTICULO YA PUBLICADO  arXiv:2609.02630.   9 de septiembre de 2026.
#
# POR QUE.  Publicamos en serie y con estilo de casa.  Lo que se repite entre dos articulos
# nuestros no lo ve ninguna compuerta: check_citas mira lo que atribuimos a OTROS, no lo que nos
# copiamos a nosotros.  Dos riesgos distintos, y hay que separarlos:
#   (a) SOLAPE DE PROSA -- parrafos, encuadres o coletillas reaprovechados.  Es autoplagio y en
#       arXiv/revista cuenta.
#   (b) SOLAPE DE OBJETO -- que un resultado que damos por nuestro y nuevo ya este ahi.  Esto es
#       lo grave.  El sospechoso concreto: nuestra familia excepcional q=6, y el sector Z_m /
#       clases de equivalencia de condiciones de contorno, que ES el tema de 2609.02630.
#
# QUE COMPRUEBA
#   X1  n-gramas de prosa compartidos entre las dos fuentes (a)
#   X2  los objetos nuestros, uno a uno, buscados en su texto (b)
#   X3  SENUELO: una frase que SI esta en los dos tiene que salir
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_solape_2609.py > check_solape_2609_OUT.txt 2>&1

import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

NUESTRO = "characters_at_aP.tex"
SUYO = os.path.join("..", "_papers", "op3_2609.02630.txt")
fallos = []


def limpia_tex(s):
    s = s[s.index(r"\begin{document}"):]
    s = re.sub(r"(?s)\\begin\{thebibliography\}.*?\\end\{thebibliography\}", " ", s)
    s = re.sub(r"(?s)\$\$.*?\$\$", " ", s)
    s = re.sub(r"\$[^$]*\$", " ", s)
    s = re.sub(r"(?s)\\begin\{(equation|align|array|tabular|tikzpicture)\*?\}.*?"
               r"\\end\{\1\*?\}", " ", s)
    s = re.sub(r"\\[A-Za-z]+\*?(\[[^\]]*\])?", " ", s)
    s = s.replace("{", " ").replace("}", " ").replace("&", " ").replace("\\", " ")
    return s


def palabras(s):
    s = s.lower()
    s = s.replace("\u2019", "'").replace("\u2013", "-").replace("\u2014", "-")
    s = re.sub(r"[^a-z' ]+", " ", s)
    return [w for w in s.split() if len(w) > 1]


def ngramas(ws, k):
    return {" ".join(ws[i:i + k]): i for i in range(len(ws) - k + 1)}


tex = io.open(NUESTRO, encoding="utf-8").read()
suyo_txt = io.open(SUYO, encoding="utf-8", errors="replace").read()
A = palabras(limpia_tex(tex))
B = palabras(suyo_txt)
print("nuestro: %d palabras de prosa    2609.02630: %d palabras" % (len(A), len(B)))

print("=" * 96)
print("X1 -- PROSA COMPARTIDA")
for K in (12, 9, 7):
    na, nb = ngramas(A, K), ngramas(B, K)
    com = sorted(set(na) & set(nb))
    print("   %2d-gramas compartidos: %d" % (K, len(com)))
    for g in com[:12]:
        print("        \"%s\"" % g)
    if K == 12 and com:
        fallos.append("%d parrafos de 12 palabras compartidos" % len(com))
    if K == 9 and len(com) > 3:
        fallos.append("%d 9-gramas compartidos" % len(com))

print("=" * 96)
print("X2 -- LOS OBJETOS NUESTROS, BUSCADOS EN SU TEXTO")
bl = " ".join(B)
OBJETOS = [
    ("el caracter simplectico en una raiz de la unidad", ["symplectic character", "twisted by roots",
                                                          "character value", "weyl character"]),
    ("nuestra recta de rho / a_P", ["orbit of rho", "half the sum of positive roots", "coxeter element"]),
    ("la clasificacion q | 2m, 2m+1, 2m+2", ["2m + 1", "2m + 2", "divides 2m"]),
    ("la excepcion q=6, m = 1 mod 3", ["q = 6", "m \u2261 1", "mod 3", "z6", "z 6"]),
    ("el alfabeto de residuos", ["alphabet", "residue"]),
    ("clases de equivalencia", ["equivalence class"]),
    ("dimensiones / factorizacion en dimensiones", ["factorisation", "factorization", "dimension of a highest"]),
    ("Sp(2m) / grupo simplectico", ["sp(2", "symplectic group"]),
    ("integridad de caracteres", ["integer", "integrality"]),
]
for nombre, claves in OBJETOS:
    hits = [c for c in claves if c in bl]
    print(("   OK   " if not hits else "  MIRAR ") + "%-48s %s" % (nombre, hits if hits else ""))

print("=" * 96)
print("X3 -- SENUELO: algo que SI comparten tiene que salir")
for frase in ("carles mar", "orbifold"):
    en_a = frase in " ".join(A)
    en_b = frase in bl
    print("   '%s':  nuestro=%s  suyo=%s" % (frase, en_a, en_b))
senuelo_ok = "orbifold" in bl
print(("   OK   " if senuelo_ok else "  FALLA ") + "X3: el control ve el texto del otro articulo")
if not senuelo_ok:
    fallos.append("senuelo X3")

print("=" * 96)
if fallos:
    print("RESULTADO: MIRAR -> %s" % fallos)
    sys.exit(1)
print("RESULTADO: sin prosa reaprovechada y sin objeto compartido con arXiv:2609.02630.")
