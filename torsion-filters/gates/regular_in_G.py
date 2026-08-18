# -*- coding: utf-8 -*-
# (R_C): EL FILTRO ES REGULARIDAD EN EL GRUPO ORIGINAL, NO EN EL DUAL.   16 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 22 de la consulta.  Nosotros habiamos leido la pared que sobra en el caso
# par como "la pared AFIN del nivel".  El corrige, y mejor:
#
#   * el elemento que hay que mirar es  x_{eta,t} = (xi^{a_1},...,xi^{a_m}) en el grupo ORIGINAL,
#     con  a_i = eta_i + rho_i;
#   * las raices de C_m son  e_i - e_j,  e_i + e_j,  2 e_i,  y evaluan a  xi^{a_i-a_j},
#     xi^{a_i+a_j},  xi^{2 a_i};
#   * luego x no es regular  <=>  a_i = +-a_j  o  2 a_i = 0  (mod t)  --  que son EXACTAMENTE las
#     tres paredes del Lema 3.1.
#
#       (R)   tau_t(eta) != 0   <=>   (eta + rho)(xi)  es regular semisimple en G.
#
#   * NPP mide la regularidad en el DUAL.  Al pasar C_m -> B_m la raiz LARGA 2e_i se vuelve la
#     CORTA e_i, y la condicion 2a_i = 0 se degrada a a_i = 0.  Con t impar da igual porque 2 es
#     invertible; con t par aparece la solucion extra a_i = t/2.  Esa es la discrepancia entera.
#
# POR QUE ESTE GUION NO USA SAGE.  Porque puede, y porque asi es una verificacion INDEPENDIENTE de
# los numeros que dio Sage: aritmetica exacta en Z[x]/(Phi_t(x)) con sympy, cero librerias de
# teoria de representaciones.  Si los dos instrumentos coinciden, el dato es del objeto.
#
# LO QUE SE MIDE
#   R1  tau != 0  (por el determinante de Weyl, numerador exacto)  contra  "(eta+rho)(xi) regular
#       en G", raiz por raiz del sistema ORIGINAL.  Las dos direcciones por separado.
#   R2  lo mismo contra "regular en el DUAL" (la lectura de NPP).  Tiene que coincidir en t impar
#       y fallar en t par, y los fallos tienen que ser exactamente los a_i = t/2.
#   R3  el numero de fallos, y el testigo minimo de cada t par.
#
# CONTROLES
#   C0  el denominador de Weyl (eta = 0) tiene que ser NO NULO en todos los t: si fuera cero el
#       cociente no existiria y todo lo demas seria ruido.
#   C1  SEÑUELO: la misma equivalencia con las raices de A_{m-1} (solo e_i - e_j).  Tiene que
#       fallar, porque se deja dos de las tres familias.
#   C2  n impreso siempre; y para t par se cuenta cuantos eta caen SOLO por la pared 2a=0.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python regular_in_G.py

import json
import sys
from itertools import combinations

from sympy import Poly, cyclotomic_poly, symbols, ZZ

x = symbols('x')


def make_ring(t):
    """aritmetica exacta en Z[x]/(Phi_t(x)):  xi^k -> x^k reducido."""
    phi = Poly(cyclotomic_poly(t, x), x, domain=ZZ)
    pow_cache = {}

    def xp(k):
        k %= t
        if k not in pow_cache:
            pow_cache[k] = Poly(x ** k, x, domain=ZZ).rem(phi)
        return pow_cache[k]
    return phi, xp


def det_poly(M, phi):
    """determinante por eliminacion de Laplace, exacto, matrices pequeñas."""
    n = len(M)
    if n == 1:
        return M[0][0]
    tot = Poly(0, x, domain=ZZ)
    for j in range(n):
        if M[0][j].is_zero:
            continue
        sub = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        term = (M[0][j] * det_poly(sub, phi)).rem(phi)
        tot = tot + term if j % 2 == 0 else tot - term
    return tot.rem(phi)


def numerador(b, t, phi, xp):
    """det( zeta^{i b_j} - zeta^{-i b_j} ),  i,j = 1..m."""
    m = len(b)
    M = [[(xp((i + 1) * b[j]) - xp(-(i + 1) * b[j])).rem(phi) for j in range(m)]
         for i in range(m)]
    return det_poly(M, phi)


def particiones(k, maxlen):
    if maxlen == 0:
        if k == 0:
            yield ()
        return
    if k == 0:
        yield (0,) * maxlen
        return
    for first in range(k, 0, -1):
        for resto in particiones(k - first, maxlen - 1):
            if not resto or first >= resto[0]:
                yield (first,) + tuple(resto)


RES = []
print("=" * 122)
print("(R): EL FILTRO ES REGULARIDAD EN EL GRUPO ORIGINAL.  Verificacion independiente, sin Sage.")
print("=" * 122)
print("  t | G     | rk | eta probados | tau!=0 | (R) en G ORIGINAL | (R) en el DUAL | señuelo A_{m-1}")
print("  " + "-" * 116)

for t in range(3, 13):
    if t % 2 == 0:
        tipo, rk = "C", (t - 2) // 2
        # a_j = eta_j + rho_j,  rho_C = (m, m-1, ..., 1);  zeta = xi;  b = a
        doble = 1
        f_orig, f_dual = 2, 1
    else:
        tipo, rk = "B", (t - 1) // 2
        # A_j = 2 a_j = 2 eta_j + 2(m'-j) + 1;  zeta = xi^{(t+1)/2}, primitiva;  b = 2a
        doble = 2
        f_orig, f_dual = 1, 1
    if rk < 1:
        continue
    phi, xp = make_ring(t)

    # C0: el denominador
    if doble == 1:
        b0 = [rk - j for j in range(rk)]
    else:
        b0 = [2 * (rk - j - 1) + 1 for j in range(rk)]
    den = numerador(b0, t, phi, xp)
    if den.is_zero:
        print("  %2d | %s%-4d| C0 FALLA: el denominador de Weyl es CERO" % (t, tipo, rk))
        continue

    n = ok_G = ok_dual = ok_sen = 0
    vivos = 0
    fallos_dual = []
    for k in range(0, 3 * t + 1):
        for e in particiones(k, rk):
            eta = tuple(list(e) + [0] * (rk - len(e)))[:rk]
            if len(eta) != rk:
                continue
            n += 1
            if doble == 1:
                b = [eta[j] + rk - j for j in range(rk)]
            else:
                b = [2 * eta[j] + 2 * (rk - j - 1) + 1 for j in range(rk)]
            vivo = not numerador(b, t, phi, xp).is_zero
            if vivo:
                vivos += 1
            pares = all((b[i] - b[j]) % t != 0 and (b[i] + b[j]) % t != 0
                        for i, j in combinations(range(rk), 2))
            # (R) en el grupo ORIGINAL.  C_m: la raiz LARGA es 2e_i  ->  2 a_i = 0, y b = a.
            #                            B_m': la raiz corta es e_i ->   a_i = 0, y b = 2a,
            #                            que con 2 invertible mod t impar es b_i = 0.
            reg_G = pares and all((f_orig * b[i]) % t != 0 for i in range(rk))
            # (R) en el DUAL.  C_m^v = B_m: la larga 2e_i se vuelve la corta e_i -> a_i = 0.
            reg_D = pares and all((f_dual * b[i]) % t != 0 for i in range(rk))
            # SEÑUELO: solo las raices de A_{rk-1}
            reg_A = all((b[i] - b[j]) % t != 0 for i, j in combinations(range(rk), 2))
            if reg_G == vivo:
                ok_G += 1
            if reg_D == vivo:
                ok_dual += 1
            else:
                fallos_dual.append((eta, [v % t for v in b]))
            if reg_A == vivo:
                ok_sen += 1
    print("  %2d | %s_%-3d | %2d | %12d | %6d | %17s | %14s | %s"
          % (t, tipo, rk, rk, n, vivos,
             "%d / %d" % (ok_G, n), "%d / %d" % (ok_dual, n), "%d / %d" % (ok_sen, n)))
    if fallos_dual:
        eta, bm = fallos_dual[0]
        print("       el dual falla en %d; testigo minimo eta=%s con b mod %d = %s"
              % (len(fallos_dual), str(eta), t, str(bm)))
    sys.stdout.flush()
    RES.append({"t": int(t), "tipo": tipo, "rango": int(rk), "n_eta": int(n), "vivos": int(vivos),
                "R_original": int(ok_G), "R_dual": int(ok_dual), "senuelo_A": int(ok_sen),
                "n_fallos_dual": int(len(fallos_dual)),
                "testigo_dual": {"eta": [int(v) for v in fallos_dual[0][0]],
                                 "b_mod_t": [int(v) for v in fallos_dual[0][1]]} if fallos_dual else None})

print("")
print("=" * 122)
print("  LECTURA, escrita ANTES de correr:")
print("   * si (R) en el ORIGINAL acierta el 100 % en las dos paridades, el Lema 3.1 se reenuncia")
print("     como 'regularidad de (eta+rho)(xi) en G' y deja de necesitar tres casos.")
print("   * si el DUAL acierta en t impar y falla en t par, la discrepancia B/C esta medida.")
print("   * si el señuelo A_{m-1} acierta parecido, las otras dos familias de raices no decidian.")
json.dump(RES, open("regular_in_G_RESULT.json", "w"), indent=1)
print("=" * 122)
print("DONE")
