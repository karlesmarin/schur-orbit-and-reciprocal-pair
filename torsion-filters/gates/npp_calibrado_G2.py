# -*- coding: utf-8 -*-
# CALIBRAR EL INSTRUMENTO CONTRA EL EJEMPLO DE G2 DE NPP.   20 de agosto de 2026.
#
# POR QUE, Y POR QUE VA ANTES.  Lo primero es la traduccion a los terminos de NPP, calibrada
# contra un ejemplo DE ELLOS --- el G2 de la §7 de arXiv:2504.14684 --- y ANTES de mirar nuestra
# familia.  Sin esto, los 72/72 y los 22 casos de igualdad con caracter cero no son una medida:
# son una traduccion sin comprobar.
#
# LA §7 DE NPP, LEIDA LITERAL DEL PDF (arXiv:2504.14684, pp. 19-20).  Raices positivas de G2 en la base
# {alpha1, alpha2} de Bourbaki, con sus CORRAICES, y los emparejamientos con lambda+rho para
# lambda = k.omega1 + l.omega2:
#     <l+r, a1^> = k+1        <l+r, a2^> = l+1        <l+r, a3^> = k+3l+4
#     <l+r, a4^> = 2k+3l+5    <l+r, a5^> = k+l+2      <l+r, a6^> = k+2l+3
# Y su Proposicion 7.1, con Theta_{k,l}(C2) en la unica clase de orden 2:
#     0                          si k, l impares
#     (k+l+2)(3l+k+4)/8          si k, l pares
#     -(k+1)(k+2l+3)/8           si k impar, l par
#     -(l+1)(3l+2k+5)/8          si k par, l impar
#
# QUE SE CALIBRA
#   G1  los seis emparejamientos, y el reparto de PARIDAD por casos, contra su analisis literal:
#       los seis pares si k,l impares; i=3,5 si ambos pares; i=2,4 si k par l impar; i=1,6 si k
#       impar l par.  Si esto falla, he leido mal su tabla y todo lo demas sobra.
#   G2  dim Z_Ghat((lambda.rho)(-1)) = rango + #{raices con <l+r,alpha^> par}, contando las raices
#       en los dos signos.  Y dim Z_G(C2) = 6, que su propio rho reproduce: 2 + 2*2 = 6.
#   G3  FATAL: su dicotomia contra su propia Proposicion 7.1 ---
#           Theta = 0   <=>   dim Z_Ghat > dim Z_G = 6.
#       Es el enunciado de la Question 8.1 medido dentro de un caso que ellos ya resolvieron.
#   G4  y la mitad de IGUALDAD: cuando las dimensiones coinciden, .es el valor no nulo, y es el
#       producto de los emparejamientos pares dividido por el de rho (que vale 8)?
#   D1  SENUELO: emparejar con la RAIZ en vez de con la CORRAIZ.  En G2 las dos difieren por el
#       factor 3 en las raices largas, asi que tiene que romper --- y si no rompe, mi lectura no
#       distingue raiz de corraiz, que es justo el eje del enunciado.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python npp_calibrado_G2.py > npp_calibrado_G2_OUT.txt 2>&1

import json
import sys
from collections import Counter
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

KMAX = 12
RANGO_G2 = 2
DIM_G2 = 14


def pares_corraiz(k, l):
    """<lambda+rho, alpha_i^> para las seis raices positivas, en el orden de su §7."""
    return [k + 1, l + 1, k + 3 * l + 4, 2 * k + 3 * l + 5, k + l + 2, k + 2 * l + 3]


def pares_raiz(k, l):
    """SENUELO: emparejar con la RAIZ.  alpha2, alpha5, alpha6 son las que llevan el 3."""
    c = pares_corraiz(k, l)
    return [c[0], 3 * c[1], c[2], c[3], 3 * c[4], 3 * c[5]]


def theta_suya(k, l):
    """Su Proposicion 7.1, transcrita."""
    if k % 2 == 1 and l % 2 == 1:
        return 0
    if k % 2 == 0 and l % 2 == 0:
        return Fraction((k + l + 2) * (3 * l + k + 4), 8)
    if k % 2 == 1 and l % 2 == 0:
        return Fraction(-(k + 1) * (k + 2 * l + 3), 8)
    return Fraction(-(l + 1) * (3 * l + 2 * k + 5), 8)


print("=" * 100)
print("CALIBRADO CONTRA EL G2 DE SU §7")
print("=" * 100)

# ---------------------------------------------------------------- G1
print("")
print("G1  el reparto de paridad por casos, contra su analisis literal")
ESPERADO = {(1, 1): {0, 1, 2, 3, 4, 5},   # los seis
            (0, 0): {2, 4},               # i = 3, 5  (1-indexado)
            (0, 1): {1, 3},               # i = 2, 4
            (1, 0): {0, 5}}               # i = 1, 6
g1_ok, g1_mal = 0, []
tot = 0
for k in range(0, KMAX + 1):
    for l in range(0, KMAX + 1):
        tot += 1
        v = pares_corraiz(k, l)
        pares = {i for i, x in enumerate(v) if x % 2 == 0}
        esp = ESPERADO[(k % 2, l % 2)]
        if pares == esp:
            g1_ok += 1
        elif len(g1_mal) < 4:
            g1_mal.append((k, l, sorted(pares), sorted(esp)))
print("     %d de %d   %s" % (g1_ok, tot, "COINCIDE con su tabla" if g1_ok == tot
                              else "*** FALLA: %s ***" % (g1_mal,)))

# ---------------------------------------------------------------- G2
v_rho = pares_corraiz(0, 0)
h_rho = sum(1 for x in v_rho if x % 2 == 0)
dim_ZG = RANGO_G2 + 2 * h_rho
print("")
print("G2  dim Z_G(C2) desde rho : rango 2 + 2*%d = %d   %s"
      % (h_rho, dim_ZG, "= 6, el centralizador SL2xSL2/mu2  OK" if dim_ZG == 6 else "*** != 6 ***"))
print("     y el producto de los emparejamientos pares de rho : %d   %s"
      % (__import__("math").prod([x for x in v_rho if x % 2 == 0]),
         "= 8, el denominador de su Prop 7.1" if
         __import__("math").prod([x for x in v_rho if x % 2 == 0]) == 8 else "*** != 8 ***"))

# ---------------------------------------------------------------- G3 y G4
print("")
print("G3  FATAL  la dicotomia de NPP contra su propia Proposicion 7.1")
c3 = Counter()
c4_ok, c4_tot, c4_mal = 0, 0, []
d1 = Counter()
mal3 = []
for k in range(0, KMAX + 1):
    for l in range(0, KMAX + 1):
        v = pares_corraiz(k, l)
        h = sum(1 for x in v if x % 2 == 0)
        dimZGhat = RANGO_G2 + 2 * h
        pred_cero = dimZGhat > dim_ZG
        th = theta_suya(k, l)
        c3[(pred_cero, th == 0)] += 1
        if pred_cero != (th == 0) and len(mal3) < 5:
            mal3.append((k, l, dimZGhat, th))
        if dimZGhat == dim_ZG:
            c4_tot += 1
            prod = __import__("math").prod([x for x in v if x % 2 == 0])
            if abs(th) == Fraction(prod, 8) and th != 0:
                c4_ok += 1
            elif len(c4_mal) < 4:
                c4_mal.append((k, l, str(th), prod))
        vs = pares_raiz(k, l)
        hs = sum(1 for x in vs if x % 2 == 0)
        d1[(RANGO_G2 + 2 * hs > dim_ZG, th == 0)] += 1

n = sum(c3.values())
ac3 = c3[(True, True)] + c3[(False, False)]
ad1 = d1[(True, True)] + d1[(False, False)]
print("     predice cero y ES cero : %4d      predice cero y NO lo es : %4d"
      % (c3[(True, True)], c3[(True, False)]))
print("     predice no-cero y es 0 : %4d      predice no-cero y no lo es: %4d"
      % (c3[(False, True)], c3[(False, False)]))
print("     acierta %d de %d   %s" % (ac3, n, "EQUIVALENCIA -- el instrumento calibra"
                                      if ac3 == n else "*** FALLA en %d: %s ***" % (n - ac3, mal3)))
print("")
print("G4  la mitad de IGUALDAD: |Theta| == (producto de los pares)/8")
print("     %d de %d   %s" % (c4_ok, c4_tot, "COINCIDE" if c4_ok == c4_tot
                              else "*** falla en %d: %s ***" % (c4_tot - c4_ok, c4_mal)))
print("")
print("D1  SENUELO emparejar con la RAIZ en vez de la CORRAIZ : %d de %d   %s"
      % (ad1, n, "*** NO DISCRIMINA: mi lectura no distingue raiz de corraiz ***" if ad1 == n
         else "falla en %d --- DISCRIMINA, la corraiz es la buena" % (n - ad1)))

json.dump(dict(G1=[g1_ok, tot], dim_ZG=dim_ZG, G3={str(a): int(b) for a, b in c3.items()},
               G3_ac=int(ac3), G4=[int(c4_ok), int(c4_tot)], D1=int(ad1), n=int(n)),
          open("npp_calibrado_G2_DUMP.json", "w"), indent=1)
print("")
print("=" * 100)
print("DONE")
