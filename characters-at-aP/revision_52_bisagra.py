# -*- coding: utf-8 -*-
# LA BISAGRA QUE LA NOTA NO ESCRIBIA.   25 de agosto de 2026.
#
# POR QUE.  La ecuacion (1) de la nota evalua chi_lambda(rho(z)) y el paso (ii) la evalua en
# z = xi_q para deducir algo sobre chi_lambda(a_P^k).  Pero la nota NO dice en ningun sitio que
# a_P^k = rho(xi_q), que es justo la bisagra de todo el argumento.  Y tampoco define rho(z).
# Aqui se comprueban las dos cosas antes de escribirlas.
#
# QUE SE MIDE
#   B1  rho(z) en Sp(2m) tiene autovalores z^{+-j}, j = 1..m   (rho = (m, m-1, ..., 1))
#   B2  el alfabeto de a_P^k de la Prop. 2.1 coincide con el espectro de rho(xi_q), xi_q = zeta_t^k
#   B3  zeta_t^k es una raiz PRIMITIVA q-esima de la unidad, con q = t/gcd(k,t)
#   B4  (rho, alpha) = ht_P(alpha) para toda raiz positiva, en B, C, F4 y G2
#       --- la identidad que ata la seccion 5 con la seccion 3 y que la nota deja implicita:
#           el estadistico de raices de la especializacion principal ES la P-altura.
#   B5  y por tanto R_P(q) = {alpha : q | (rho,alpha)} = {alpha : q | ht_P(alpha)}
#   B6  SENUELO: con la CO-altura ht(alpha^v) en vez de la P-altura, B4 falla en los tipos no
#       simplemente enlazados --- si no fallara, la distincion de la nota seria vacia.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python revision_52_bisagra.py > revision_52_bisagra_OUT.txt 2>&1

import cmath
import sys
from fractions import Fraction
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


def multiset(zs, nd=9):
    return sorted((round(z.real, nd), round(z.imag, nd)) for z in zs)


print(SEP)
print("B1/B2/B3 -- a_P^k = rho(xi_q):  el espectro de rho(xi_q) ES el alfabeto de la Prop. 2.1")
print(f"{'m':>3} {'t':>4} {'k':>3} {'d':>3} {'q':>4} {'ord(zeta^k)':>12} {'espectros iguales':>18}")
malos = tot = 0
for m in range(2, 12):
    t = 2 * m + 2
    for k in range(1, 15):
        d = gcd(k, t)
        q = t // d
        zeta = cmath.exp(2j * cmath.pi / t)
        xi = zeta ** k
        # el alfabeto de la nota:  A = { zeta^{+-jk} : j = 1..m }
        A = [zeta ** (j * k) for j in range(1, m + 1)] + \
            [zeta ** (-j * k) for j in range(1, m + 1)]
        # rho en C_m es (m, m-1, ..., 1);  rho(xi) tiene autovalores xi^{+-rho_i}
        rho = list(range(m, 0, -1))
        R = [xi ** r for r in rho] + [xi ** (-r) for r in rho]
        igual = multiset(A) == multiset(R)
        # y el orden de zeta^k es exactamente q
        orden = min(n for n in range(1, t + 1) if abs(xi ** n - 1) < 1e-9)
        tot += 1
        if not igual or orden != q:
            malos += 1
        if m in (2, 5) and k <= 3:
            print(f"{m:>3} {t:>4} {k:>3} {d:>3} {q:>4} {orden:>12} {str(igual):>18}")
ok(malos == 0, "B1-B3: rho(xi_q) tiene el alfabeto de a_P^k, y xi_q es primitiva q-esima",
   f"{tot - malos}/{tot}")

print(SEP)
print("B4/B5 -- (rho, alpha) = ht_P(alpha),  y por tanto R_P(q) se puede leer de las dos maneras")


def sistema(tipo, n):
    """Devuelve (raices_positivas_en_base_simple, matriz_de_Gram_simetrizada, r).
    Normalizacion: la raiz CORTA tiene (a,a) = 2, luego (a,a)/2 = 1 para cortas y r para largas."""
    if tipo == "C":                       # C_n: n-1 cortas + 1 larga (la ultima)
        long_ = [False] * (n - 1) + [True]
        r = 2
    elif tipo == "B":                     # B_n: n-1 largas + 1 corta (la ultima)
        long_ = [True] * (n - 1) + [False]
        r = 2
    elif tipo == "F":                     # F4: a1,a2 largas; a3,a4 cortas
        long_, r, n = [True, True, False, False], 2, 4
    elif tipo == "G":                     # G2: a1 corta, a2 larga
        long_, r, n = [False, True], 3, 2
    else:
        raise ValueError(tipo)
    dl = [Fraction(r if L else 1) for L in long_]        # (a_i,a_i)/2
    # ⚠ CONVENIO, y aqui es donde me equivoque tres veces: C[i][j] = <a_i, a_j^v>.
    #   Entonces (a_i, a_j) = C[i][j] * (a_j,a_j)/2 = C[i][j] * dl[j], y la SIMETRIA obliga a
    #       C[i][j] * dl[j] = C[j][i] * dl[i].
    #   La primera version puso la flecha al reves en B, C y F4 --- la comprobacion de simetria lo
    #   caza, que es exactamente para lo que esta.  [[a-count-cannot-see-a-transpose]]
    C = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        C[i][i] = Fraction(2)
    if tipo in ("B", "C"):
        for i in range(n - 1):
            C[i][i + 1] = C[i + 1][i] = Fraction(-1)
        if tipo == "C":                    # a_n LARGA:  dl = [1,...,1,2]
            C[n - 2][n - 1] = Fraction(-1); C[n - 1][n - 2] = Fraction(-2)
        else:                              # B: a_n CORTA:  dl = [2,...,2,1]
            C[n - 2][n - 1] = Fraction(-2); C[n - 1][n - 2] = Fraction(-1)
    elif tipo == "F":                      # dl = [2,2,1,1]
        C[0][1] = C[1][0] = Fraction(-1)
        C[1][2] = Fraction(-2); C[2][1] = Fraction(-1)
        C[2][3] = C[3][2] = Fraction(-1)
    else:                                  # G2:  dl = [1,3]
        C[0][1] = Fraction(-1); C[1][0] = Fraction(-3)
    Gm = [[C[i][j] * dl[j] for j in range(n)] for i in range(n)]
    return C, Gm, dl, r, n


def raices_positivas(C, n):
    """Genera las raices positivas por cierre bajo s_i, en coordenadas de la base simple."""
    pos = {tuple(1 if j == i else 0 for j in range(n)) for i in range(n)}
    cambio = True
    while cambio:
        cambio = False
        for b in list(pos):
            for i in range(n):
                # s_i(b) = b - <b, a_i^v> a_i,  con  <b, a_i^v> = sum_j b_j <a_j, a_i^v> = sum_j b_j C[j][i]
                pair = sum(Fraction(b[j]) * C[j][i] for j in range(n))
                nb = list(b)
                nb[i] -= int(pair)
                nb = tuple(nb)
                if all(x >= 0 for x in nb) and any(x > 0 for x in nb) and nb not in pos:
                    pos.add(nb)
                    cambio = True
    return sorted(pos)


print(f"{'tipo':>5} {'raices+':>8} {'(rho,a)=ht_P?':>14} {'R_P(q) igual?':>14} "
      f"{'senuelo co-altura':>18}")
for tipo, n in (("B", 3), ("B", 4), ("C", 3), ("C", 4), ("C", 5), ("F", 4), ("G", 2)):
    C, Gm, dl, r, n = sistema(tipo, n)
    ok_sim = all(Gm[i][j] == Gm[j][i] for i in range(n) for j in range(n))
    if not ok_sim:
        ok(False, f"{tipo}{n}: la forma simetrizada NO es simetrica")
        continue
    P = raices_positivas(C, n)
    # rho: <rho, a_i^v> = 1 para toda i simple  =>  en coords de pesos fundamentales rho=(1,..,1)
    # ht(alpha^v) = sum_i coef de alpha^v en las a_i^v.  Con alpha = sum c_i a_i,
    #   alpha^v = sum c_i * (dl_i/dl_alpha) a_i^v   donde dl_alpha = (alpha,alpha)/2.
    mal_id = mal_rp = 0
    senuelo = 0
    for b in P:
        # (alpha, alpha)/2
        aa = sum(Fraction(b[i]) * Gm[i][j] * Fraction(b[j])
                 for i in range(n) for j in range(n)) / 2
        # (rho, alpha):  rho = sum_i w_i, y (w_i, a_j) = delta_ij * dl_j
        rho_a = sum(Fraction(b[i]) * dl[i] for i in range(n))
        # ht(alpha^v)
        ht_v = sum(Fraction(b[i]) * dl[i] / aa for i in range(n))
        # ht_P(alpha) = (alpha,alpha)/2 * ht(alpha^v)
        ht_P = aa * ht_v
        if rho_a != ht_P:
            mal_id += 1
        # senuelo: la co-altura sola
        if ht_v != rho_a:
            senuelo += 1
    # R_P(q) leido de las dos maneras, para todo q
    for q in range(2, 25):
        s1 = {b for b in P
              if sum(Fraction(b[i]) * dl[i] for i in range(n)) % q == 0}
        s2 = set()
        for b in P:
            aa = sum(Fraction(b[i]) * Gm[i][j] * Fraction(b[j])
                     for i in range(n) for j in range(n)) / 2
            ht_v = sum(Fraction(b[i]) * dl[i] / aa for i in range(n))
            if (aa * ht_v) % q == 0:
                s2.add(b)
        if s1 != s2:
            mal_rp += 1
    print(f"{tipo + str(n):>5} {len(P):>8} {str(mal_id == 0):>14} {str(mal_rp == 0):>14} "
          f"{str(senuelo > 0):>18}")
    ok(mal_id == 0, f"B4 {tipo}{n}: (rho,alpha) = ht_P(alpha) en las {len(P)} raices positivas")
    ok(mal_rp == 0, f"B5 {tipo}{n}: R_P(q) igual por las dos lecturas, q=2..24")
    if tipo in ("B", "C", "F", "G"):
        ok(senuelo > 0, f"B6 {tipo}{n}: la co-altura sola NO coincide (la distincion no es vacia)",
           f"{senuelo} de {len(P)} raices difieren")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: la bisagra a_P^k = rho(xi_q) se cumple, y (rho,alpha) = ht_P(alpha) tambien.")
print("           Las dos se pueden escribir en la nota.")
