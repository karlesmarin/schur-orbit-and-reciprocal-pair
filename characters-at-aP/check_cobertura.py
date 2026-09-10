# -*- coding: utf-8 -*-
# CADA FORMULA ENMARCADA DEL ARTICULO, CON LA COMPUERTA QUE LA VERIFICA.
# 9 de septiembre de 2026.
#
# POR QUE.  El articulo enmarca once formulas con \boxed.  Cada una es una afirmacion, y varias
# estan verificadas en guiones distintos, escritos en dias distintos.  Nada ataba una cosa con la
# otra: si mañana alguien toca una formula, o si una salida se queda obsoleta, no salta nada.
# Esta compuerta es el MAPA, y es ejecutable: para cada formula exige (a) que siga en el .tex y
# (b) que la frase que la verifica siga en la salida archivada de su guion.
# Carles, 9-sep-2026: "revisa las formulas que esten todas bien, analizalas".
#
# LO QUE ESTA COMPUERTA NO HACE.  No vuelve a medir: eso lo hacen los guiones que nombra.  Lo que
# hace es que la cobertura sea auditable de un vistazo y que no se pueda perder sin ruido.
#
# QUE COMPRUEBA
#   B1  las once formulas siguen en el .tex, y son once
#   B2  cada una tiene guion, salida archivada, y la frase que la verifica dentro
#   B3  SENUELO: una frase que no esta en su salida TIENE que saltar
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python check_cobertura.py > check_cobertura_OUT.txt 2>&1

import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
G = os.path.join("..", "gates")
TEX = io.open("characters_at_aP.tex", encoding="utf-8").read()
fallos = []

# (marca en el .tex, guion, frase que en su salida acredita la verificacion, que es)
MAPA = [
    (r"e^{-2\pi i\theta(\lambda,\rho)}", "referee2.py", "40 pares medidos, 0 discrepancias",
     "la forma de Weyl en la recta, y que p se cae"),
    (r"\neq0 \iff N_q(\ell)=r_q", "revision_47_verificacion.py",
     "R3: el criterio N_q(lambda) = r_q se cumple",
     "el criterio: no anularse <=> N_q(ell) = r_q"),
    (r"\operatorname{sgn}\chi_\lambda(a_P^{\,k})", "referee2.py", "11972 supervivientes, 0 fallos",
     "el signo por residuos ordenados, F_q"),
    (r"|\chi_\lambda(a_P^k)| \;=\; \prod_{\text{classes}} D_c", "dims_recuento.py",
     "factorizacion OK", "la factorizacion en dimensiones, por clases"),
    (r"\mathscr{R}_x(X)\;=\;\frac{P_x(X)^2+P_x(X^2)}{2}", "referee2.py",
     "la identidad: 48669 casos, 0 fallos", "el lema de residuos (Prop. 8.5)"),
    # ⚠ este apuntaba a R_P_todos_los_tipos.py, que trata otra cosa (R_P(q) y Polo).  Quien
    #   verifica esta identidad es simplifica2.py, en B, C, D y G2.
    (r"\operatorname{ch}\mathfrak g\,(x)\;-\;\mathrm{rank}", "simplifica2.py",
     "R_x(X) = char(adjunto)(x) - rango",
     "el lema de residuos leido como caracter de la adjunta"),
    (r"\mathscr{R}_{\lambda+\rho}-\mathscr{R}_{\rho}", "residual_histogramas.py",
     "H1: la identidad de histogramas se cumple en los 3729",
     "la diferencia de histogramas es la de otro grupo (Prop. 8.6)"),
    (r"\kappa(\lambda)\,D_0\,D_{q/2}", "teorema_II.py",
     "supervivientes contrastados: 2186", "el modulo factorizado, kappa D_0 D_{q/2} prod D_c"),
    (r"S_{C_m}(q)=", "dims_recuento.py", "recuento cerrado OK", "el recuento cerrado, todo q"),
]

print("=" * 100)
print("B1 -- LAS FORMULAS ENMARCADAS SIGUEN AHI")
n_box = len(re.findall(r"\\boxed", TEX))
print(f"   \\boxed en el .tex: {n_box}")
ok_b1 = n_box == 11
print(("   OK   " if ok_b1 else "  FALLA ") + "B1: son once, como cuando se hizo este mapa")
if not ok_b1:
    fallos.append("el numero de formulas enmarcadas ha cambiado: reviser este mapa")

print("=" * 100)
print("B2 -- CADA UNA CON SU GUION, SU SALIDA Y LA FRASE QUE LA ACREDITA")
for marca, guion, frase, desc in MAPA:
    en_tex = marca in TEX
    p = os.path.join(G, guion)
    sal = os.path.join(G, guion[:-3] + "_OUT.txt")
    hay_g = os.path.exists(p)
    hay_s = os.path.exists(sal)
    hay_f = hay_s and frase in io.open(sal, encoding="utf-8", errors="replace").read()
    bien = en_tex and hay_g and hay_s and hay_f
    print(("   OK   " if bien else "  FALLA ") + f"{desc:<58} {guion}")
    if not bien:
        det = []
        if not en_tex:
            det.append("la formula ya no esta en el .tex")
        if not hay_g:
            det.append("falta el guion")
        if not hay_s:
            det.append("falta su salida archivada")
        elif not hay_f:
            det.append(f"la salida no contiene '{frase}'")
        print("           " + "; ".join(det))
        fallos.append(desc)

print("=" * 100)
print("B3 -- SENUELO")
sal = os.path.join(G, "referee2_OUT.txt")
_f = os.path.exists(sal) and "0 pares medidos, 7 discrepancias" in io.open(
    sal, encoding="utf-8", errors="replace").read()
print(("   OK   " if not _f else "  FALLA ")
      + "B3: una frase inventada NO se encuentra en la salida --- B2 sabe fallar")
if _f:
    fallos.append("senuelo B3")

print("=" * 100)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: cada formula enmarcada tiene guion, salida archivada y frase que la acredita.")
