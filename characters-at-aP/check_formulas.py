# -*- coding: utf-8 -*-
# LAS FORMULAS QUE ATRIBUIMOS, TOKEN A TOKEN CONTRA LA FORMULA DE SU FUENTE.
# 9 de septiembre de 2026.
#
# POR QUE.  check_citas compara FRASES y check_numeracion comprueba que el NUMERO del resultado
# exista y tenga cerca la frase que le atribuimos.  Ninguna de las dos mira las FORMULAS.  Y el
# articulo reproduce formulas ajenas: el Delta de Lecouvey, el factor ortogonal de Ayyer-Kumari,
# el criterio y el signo de Cellini-Moseneder-Papi.  Ahi un solo token cambia la matematica.
# [[audit-the-formula-not-the-words]] --- Carles, 9-sep: "revisar las referencias, tambien por
# formulas".
#
# EL CASO QUE LO ILUSTRA, Y NO ES HIPOTETICO.  El keybox de la §7 dice que la clase 0 de Lecouvey
# lleva  prod_{r<=s} (1 - x_r x_s)  y que POR ESO es de tipo C.  Si fuera  r < s  seria tipo D:
# el termino r = s es  (1 - x_r^2),  que es la raiz LARGA  2e_r,  y es justo lo que separa C de D.
# Un "<" por un "<=" convertiria una atribucion correcta en una falsa sin que compile distinto,
# sin que ninguna otra compuerta lo viera, y sin que la frase cambiara.
#
# COMO SE COMPARA.  No se normaliza LaTeX contra la capa de texto de un PDF -- eso no es fiable.
# Se descompone cada formula en sus TOKENS ESENCIALES: los operandos, los exponentes y la
# CONDICION DE INDICE.  Cada token tiene que estar en la ventana de la fuente; y se declaran
# tambien los tokens PROHIBIDOS, los que indicarian la formula equivocada.
#
# QUE COMPRUEBA
#   F1  cada formula atribuida: sus tokens esenciales, en la ventana de su fuente
#   F2  ...y ningun token prohibido en esa misma ventana
#   F3  que la formula que el ARTICULO imprime lleva esos mismos tokens
#   F4  SENUELOS: cambiar un token TIENE que romper la comprobacion
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_formulas.py > check_formulas_OUT.txt 2>&1

import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
PAP = os.path.join("..", "_papers")
fallos = []


def norm(s):
    s = s.replace("’", "'").replace("ﬁ", "fi").replace("ﬂ", "fl").replace("ﬃ", "ffi")
    s = s.replace("−", "-").replace("¨", "").replace("´", "")
    return re.sub(r"\s+", " ", s).strip().lower()


FUENTES = {
    "Lec07a": "lecouvey_math_0607038.txt",
    "AK22": "arxiv_2109.11310.txt",
    "CMP06": "cellini_moseneder_papi_math_0507610.txt",
    "Kar24": "karmakar_2412.17324.txt",
}
TXT = {}
for k, v in FUENTES.items():
    p = os.path.join(PAP, v)
    if not os.path.exists(p):
        print(f"  FALLA  falta {v}")
        fallos.append(v)
        continue
    TXT[k] = norm(io.open(p, encoding="utf-8", errors="replace").read())
    print(f"   fuente {k:<8} {v}  ({len(TXT[k])} caracteres)")

TEX = norm(io.open("characters_at_aP.tex", encoding="utf-8").read())

# (clave, ancla en la fuente, ventana, tokens que TIENEN que estar, tokens PROHIBIDOS, desc)
FORMULAS = [
    ("Lec07a", "∆i(0) =", 400,
     ["r≤s", "(1 -xrxs)", "(1 -xj xi )"],
     ["r<s r,s∈i(0)"],
     "Lecouvey, Delta de la clase 0: lleva r<=s, o sea la raiz LARGA --- tipo C, no D"),
    ("Lec07a", "∆x(k) =", 300,
     ["i<j", "(1 -xj xi )"],
     ["(1 -xrxs)"],
     "Lecouvey, Delta de los pares plegados: solo i<j --- tipo A"),
    ("Lec07a", "the half sum of positive roots is equal to", 260,
     ["(1 2, ..., rp -1 2)", "type brp"],
     [],
     "Lecouvey, la semisuma del bloque de tipo B: (1/2, ..., r_p-1/2)"),
    ("AK22", "if coret(λ) is a symplectic t-core with rank r", 700,
     ["so λ( t 2 -1)(xt)", "t even", "1 t odd"],
     [],
     "AK22 Thm 2.11(2): el factor ortogonal so_{lambda(t/2-1)}(X^t), en la rama t PAR"),
    ("CMP06", "corollary 3.2", 500,
     ["/∈1 2z", "⌊2(λ+ρ,α)⌋"],
     [],
     "CMP Cor. 3.2: el criterio con 1/2 Z, y el signo con el SUELO DE 2(lambda+rho,alpha)"),
    ("Kar24", "λ0 + ρn-k = η0(λ)/2", 200,
     ["η1(λ)/2"],
     [],
     "Karmakar Thm 5.1: las dos recetas van con /2"),
]

print("=" * 100)
print("F1/F2 -- CADA FORMULA ATRIBUIDA, TOKEN A TOKEN, EN LA VENTANA DE SU FUENTE")
for clave, ancla, vent, deben, prohibidos, desc in FORMULAS:
    if clave not in TXT:
        continue
    t = TXT[clave]
    pos = [m.start() for m in re.finditer(re.escape(norm(ancla)), t)]
    if not pos:
        print("  FALLA " + desc + f"   <-- el ancla '{ancla}' no aparece")
        fallos.append(desc)
        continue
    faltan, sobran = [], []
    for tok in deben:
        if not any(norm(tok) in t[i:i + vent] for i in pos):
            faltan.append(tok)
    for tok in prohibidos:
        if any(norm(tok) in t[i:i + vent] for i in pos):
            sobran.append(tok)
    bien = not faltan and not sobran
    print(("   OK   " if bien else "  FALLA ") + desc)
    if faltan:
        print(f"           faltan: {faltan}")
        fallos.append(desc + " (faltan)")
    if sobran:
        print(f"           PROHIBIDOS presentes: {sobran}")
        fallos.append(desc + " (prohibidos)")

print("=" * 100)
print("F3 -- Y LA FORMULA QUE EL ARTICULO IMPRIME LLEVA ESOS MISMOS TOKENS")
NUESTRAS = [
    (r"\prod_{r\le s}(1-x_rx_s)", "r\\le s",
     "el articulo escribe el Delta de la clase 0 con r\\le s, no con r<s"),
    (r"(\tfrac12,\dots,r_p-\tfrac12)", "r_p-\\tfrac12",
     "y la semisuma del bloque B con r_p - 1/2"),
]
for frag, tok, desc in NUESTRAS:
    hay = norm(frag) in TEX
    print(("   OK   " if hay else "  FALLA ") + desc)
    if not hay:
        fallos.append(desc)

print("=" * 100)
print("F4 -- SENUELOS: cambiar un token TIENE que romperlo")
_t = TXT.get("Lec07a", "")
_p = [m.start() for m in re.finditer(re.escape(norm("∆i(0) =")), _t)]
_falso = any(norm("r<s r,s∈i(0)") in _t[i:i + 400] for i in _p)
print(("   OK   " if not _falso else "  FALLA ")
      + "F4: 'r<s' NO esta en la ventana de la clase 0 --- si estuviera, seria tipo D")
if _falso:
    fallos.append("senuelo F4a")
_falso2 = norm(r"\prod_{r<s}(1-x_rx_s)") in TEX
print(("   OK   " if not _falso2 else "  FALLA ")
      + "F4: y el articulo NO imprime la version con r<s")
if _falso2:
    fallos.append("senuelo F4b")
# y que la comprobacion de F3 sabe fallar
_falso3 = norm(r"\prod_{r\le s}(1-x_rx_t)") in TEX
print(("   OK   " if not _falso3 else "  FALLA ")
      + "F4: un operando cambiado (x_t por x_s) tampoco esta")
if _falso3:
    fallos.append("senuelo F4c")
# ⚠ los tres de arriba solo prueban AUSENCIA.  Falta probar que el mecanismo POSITIVO de F1 sabe
#   fallar: si la busqueda en ventana encontrara cualquier cosa, F1 daria verde siempre.
_inventado = any(norm("(1 -xrxt)") in _t[i:i + 400] for i in _p)
print(("   OK   " if not _inventado else "  FALLA ")
      + "F4: y un token INVENTADO no se encuentra en esa ventana --- F1 sabe fallar")
if _inventado:
    fallos.append("senuelo F4d")

print("=" * 100)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: las formulas que atribuimos coinciden, token a token, con las de su fuente.")
