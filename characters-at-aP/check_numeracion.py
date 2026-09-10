# -*- coding: utf-8 -*-
# CADA RESULTADO NUMERADO QUE ATRIBUIMOS: QUE EXISTA, Y QUE DIGA LO QUE DECIMOS.
# 9 de septiembre de 2026.
#
# POR QUE.  check_citas comprueba FRASES literales; no comprueba NUMEROS. El articulo atribuye
# unos treinta resultados numerados a diez articulos ajenos --- "Corollary 3.2 of CMP", "Theorem
# 2.11 of AK22", "Prop. 3.2.5 of Lec07a" --- y ninguna compuerta miraba que ese numero exista en
# esa fuente ni que sea el que dice lo que le atribuimos.  Un numero equivocado manda al lector a
# otro teorema y parece un error nuestro de matematica, no de tecleo.
# [[audit-the-formula-not-the-words]] --- Carles, 9-sep: "hemos cambiado tanto que hay que revisar
# las referencias, tambien por formulas".
#
# ⚠ LO PRIMERO QUE ENCONTRO.  Kum24 (arXiv:2212.12477) NO ESTABA en _papers/, y el articulo le
#   atribuye DOS teoremas numerados, 2.2 y 2.13.  Se citaba de memoria una fuente que no estaba
#   en la biblioteca.  Bajado el 9-sep.
#
# QUE COMPRUEBA
#   M1  cada (fuente, resultado numerado) que citamos EXISTE literalmente en esa fuente
#   M2  ...y cerca de el aparece la frase que caracteriza lo que le atribuimos
#   M3  SENUELOS: un numero que no existe, y un numero que existe pero dice otra cosa
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_numeracion.py > check_numeracion_OUT.txt 2>&1

import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
PAP = os.path.join("..", "_papers")
fallos = []


def norm(s):
    s = s.replace("’", "'").replace("“", '"').replace("”", '"')
    s = s.replace("ﬁ", "fi").replace("ﬂ", "fl").replace("ﬃ", "ffi")
    s = s.replace("−", "-")
    return re.sub(r"\s+", " ", s).strip().lower()


FUENTES = {
    "CMP06": "cellini_moseneder_papi_math_0507610.txt",
    "NPP25": "NPP25_2504.14684.txt",
    "Kar24": "karmakar_2412.17324.txt",
    "Prasad16b": "_prasad_israe.txt",
    "AK22": "arxiv_2109.11310.txt",
    "Kum24": "kumari_2212.12477.txt",
    "Lec07a": "lecouvey_math_0607038.txt",
    "Lec07b": "lecouvey_math_0703514.txt",
    "AK25": "arxiv_2501.00275.txt",
}
TXT = {}
for k, v in FUENTES.items():
    p = os.path.join(PAP, v)
    if not os.path.exists(p):
        print(f"  FALLA  falta la fuente {v} para {k}")
        fallos.append(v)
        continue
    TXT[k] = norm(io.open(p, encoding="utf-8", errors="replace").read())
    print(f"   fuente {k:<10} {v}  ({len(TXT[k])} caracteres)")

# (clave, etiqueta literal del resultado, frase que debe estar CERCA, que le atribuimos)
CITAS = [
    # ⚠ la sonda 'eigenvalue' era MIA y estaba mal: su Cor. 3.2 no la usa.  Dice, literal,
    #   "A weight lambda belongs to Palc if and only if it is dominant and (lambda+rho,alpha) not
    #   in 1/2 Z for any alpha", y luego el signo con el suelo de 2(lambda+rho,alpha).
    ("CMP06", "corollary 3.2", "if and only if it is dominant",
     "CMP Cor. 3.2 --- el criterio, literal"),
    ("CMP06", "corollary 3.2", "belongs to the root lattice",
     "CMP Cor. 3.2 --- ...y que ahi lambda cae en el reticulo de raices"),
    ("CMP06", "proposition 3.5", None, "CMP Prop. 3.5 --- nuestro conj:cap con d=1"),
    ("NPP25", "corollary 3.7", None, "NPP Cor. 3.7 --- el recuento en el estadistico de CORRAICES"),
    ("NPP25", "proposition 3.8", None, "NPP Prop. 3.8 --- los dos miembros como centralizadores"),
    ("Kar24", "theorem 5.1", "sp(2n, c)", "Karmakar Thm 5.1 --- es de Sp(2n,C)"),
    ("AK22", "theorem 2.11", "sp2tn", "AK22 Thm 2.11 --- es de Sp_{2tn}"),
    ("AK22", "theorem 2.11", "symplectic t-core", "AK22 Thm 2.11 --- su rama que factoriza"),
    ("Kum24", "theorem 2.2", None, "Kumari Thm 2.2 --- el alfabeto con una tirada parcial"),
    # ⚠ 9-sep-2026.  Aqui decia "theorem 2.13", y EL ARTICULO CITABA ESE NUMERO.  No existe: los
    #   teoremas de su §2 son 2.2, 2.5, 2.8, 2.11, 2.14, 2.15 y 2.17.  Lo que hay con ese numero
    #   es una OBSERVACION 2.13, que dice otra cosa --- que NO se recuperan los resultados de
    #   [AK22] poniendo y=0.  El que le atribuiamos es el 2.8, que es el simplectico con una
    #   letra de mas.  Este bloque nacio de encontrar eso.
    ("Kum24", "theorem 2.8", "symplectic character",
     "Kumari Thm 2.8 --- el alfabeto con una letra de mas, en el grupo SIMPLECTICO"),
    ("Kum24", "theorem 2.8", "ωt´1x, yq",
     "Kumari Thm 2.8 --- y su alfabeto lleva la letra suelta y"),
    ("Lec07a", "theorem 3.2.3", "2p -1", "Lec07a Thm 3.2.3 --- es el caso ell IMPAR"),
    ("Lec07a", "proposition 3.2.5", "sp2n", "Lec07a Prop 3.2.5 --- es de Sp_2n"),
    ("Lec07b", "theorem 4.5.1", "symplectic or orthogonal", "Lec07b Thm 4.5.1"),
]

print("=" * 96)
print("M1/M2 -- CADA RESULTADO NUMERADO, EN SU FUENTE, Y CON SU CONTENIDO CERCA")
VENTANA = 1400
for clave, etiqueta, cerca, desc in CITAS:
    if clave not in TXT:
        continue
    t = TXT[clave]
    # ⚠ La primera version usaba t.find(), o sea la PRIMERA aparicion --- que casi siempre es la
    #   mencion en la introduccion ("probamos el Teorema 2.11..."), no el enunciado.  Con eso
    #   fallaban AK22 2.11 y Lec07b 4.5.1, que son correctos.  Hay que mirar TODAS las apariciones.
    posiciones = [m.start() for m in re.finditer(re.escape(norm(etiqueta)), t)]
    if not posiciones:
        print("  FALLA " + desc + f"   <-- '{etiqueta}' NO APARECE en {FUENTES[clave]}")
        fallos.append(desc)
        continue
    if cerca is None:
        print("   OK   " + desc + f"   ('{etiqueta}' aparece {len(posiciones)} vez/veces)")
        continue
    hay = any(norm(cerca) in t[i:i + VENTANA] for i in posiciones)
    print(("   OK   " if hay else "  FALLA ") + desc
          + ("" if hay else f"   <-- '{cerca}' no esta en las {VENTANA} letras siguientes"))
    if not hay:
        fallos.append(desc)

print("=" * 96)
print("M3 -- SENUELOS")
_a = "theorem 9.9" in TXT.get("Lec07a", "x")
print(("   OK   " if not _a else "  FALLA ") + "M3: un numero inventado NO se encuentra")
if _a:
    fallos.append("senuelo M3a")
_i = TXT.get("Kar24", "").find("theorem 5.1")
_b = _i >= 0 and "so2n+1" not in TXT["Kar24"][_i:_i + 300]
print(("   OK   " if _b else "  FALLA ")
      + "M3: y la ventana discrimina --- el Thm 5.1 de Karmakar no abre hablando de SO(2n+1)")
if not _b:
    fallos.append("senuelo M3b")

print("=" * 96)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: los resultados numerados que atribuimos existen y dicen lo que les atribuimos.")
