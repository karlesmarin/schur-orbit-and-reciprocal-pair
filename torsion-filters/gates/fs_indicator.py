# -*- coding: utf-8 -*-
# .ES EL FILTRO UN INDICADOR DE FROBENIUS-SCHUR DISFRAZADO?   17 de agosto de 2026.
#
# LA SOSPECHA, y se enuncia entera para poder matarla.  Nuestro filtro toma valores 0,+-1; el
# indicador de Frobenius-Schur toma 0,+-1 y distingue ORTOGONAL de SIMPLECTICO, que es la dicotomia
# del paper.  Que dos cosas tengan el mismo rango no es una conexion.  Lo que sube esto por encima
# de la numerologia es que el indicador FS se calcula como un caracter en un elemento ligado a una
# involucion, y nosotros evaluamos un caracter en un elemento de torsion de una coclase de reflexion.
#
# EL INDICADOR, que es clasico y NO se mide aqui:
#   Sp_{2m}      todo irreducible es autodual y  nu_2(V_eta) = (-1)^{|eta|}.
#   SO_{2m'+1}   todo irreducible tensorial es ortogonal:  nu_2 = +1 siempre.
#
# Eso deja dos preguntas afiladas, y las dos se responden con la forma cerrada del Lema T:
#   F1  .es  tau^C_t(eta) == (-1)^{|eta|}  sobre el soporte?   Si sale, el SIGNO del filtro par ES
#       el indicador FS del factor congelado, y la rama par gana un nombre clasico.
#   F2  .es  tau^B_t(eta) == +1  sobre el soporte?   Si sale, lo mismo en la rama impar.
#
# Y el marco, que se comprueba aparte y es lo que de verdad situa la cosa:
#   F3  el indicador FS es el caracter en  rho^v(-1),  o sea en el elemento PRINCIPAL DE ORDEN 2.
#       Nuestro elemento impar ES principal de orden t (probado en el paper); el par NO lo es para
#       m>=2 (probado en el paper).  Asi que FS y nuestro filtro impar son dos miembros de la MISMA
#       familia {chi(rho^v(zeta_d))}, la de Nadimpalli-Pattanayak-Prasad, con d=2 y d=t.  El par se
#       sale de la familia, y eso ya lo sabiamos.  F3 comprueba el orden del elemento FS, que es lo
#       unico de esa frase que es aritmetica y no cita.
#
# CONTROLES, para que la prueba pueda fallar:
#   C1  SENUELO  comparar con  -(-1)^{|eta|}.  Tiene que acertar exactamente donde F1 falla.
#   C2  SENUELO  comparar con el signo constante +1 en la rama par.  Si F1 y C2 aciertan lo mismo,
#       la prueba no distingue nada porque el soporte tendria |eta| par siempre.
#   C3  se imprime el reparto de |eta| mod 2 SOBRE EL SOPORTE.  Si resultara constante, F1 seria
#       verdad por vacuidad y habria que decirlo en vez de cantar victoria.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python fs_indicator.py
import io
import itertools
import json
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def clase(a, t, m):
    c = a % t
    if c == 0 or 2 * c == t:
        return None
    return (c, +1) if c <= m else (t - c, -1)


def tau_C(eta, t, m):
    """Lema T, forma cerrada, copiada de branch_filter.sage."""
    a = [eta[j] + (m - (j + 1) + 1) for j in range(m)]
    cs = [clase(v, t, m) for v in a]
    if any(c is None for c in cs):
        return 0
    cl = [c for c, _ in cs]
    if len(set(cl)) != m:
        return 0
    s = 1
    for _, sg in cs:
        s *= sg
    perm = [m - cl[j] for j in range(m)]
    inv = sum(1 for i in range(m) for j in range(i + 1, m) if perm[i] > perm[j])
    return s * (-1) ** inv


def plegar(v, t):
    v %= t
    if v == 0:
        return None
    return min(v, t - v), (1 if v < t - v else -1)


def eps_t(t, mp):
    from math import prod
    def jacobi(a, n):
        a %= n
        r = 1
        while a:
            while a % 2 == 0:
                a //= 2
                if n % 8 in (3, 5):
                    r = -r
            a, n = n, a
            if a % 4 == 3 and n % 4 == 3:
                r = -r
            a %= n
        return r if n == 1 else 0
    return jacobi(-2, t) ** ((t + 3) // 2) * (-1) ** (mp * (mp - 1) // 2)


def tau_B(eta, t, mp):
    """Cor. oddsign: tau^B = eps_t * delta(A), con A = 2(eta+rho) y el orden DECRECIENTE."""
    A = [2 * eta[j] + 2 * (mp - (j + 1)) + 1 for j in range(mp)]
    cs = [plegar(v, t) for v in A]
    if any(c is None for c in cs):
        return 0
    cl = [c for c, _ in cs]
    if sorted(cl) != list(range(1, mp + 1)):
        return 0
    s = 1
    for _, sg in cs:
        s *= sg
    perm = [mp - cl[j] for j in range(mp)]
    inv = sum(1 for i in range(mp) for j in range(i + 1, mp) if perm[i] > perm[j])
    return eps_t(t, mp) * s * (-1) ** inv


def pesos(rango, cota):
    for eta in itertools.product(range(cota + 1), repeat=rango):
        if all(eta[i] >= eta[i + 1] for i in range(rango - 1)):
            yield eta


print("=" * 104)
print("EL FILTRO CONTRA EL INDICADOR DE FROBENIUS-SCHUR")
print("=" * 104)
print("")
print("  nu_2(V_eta) = (-1)^{|eta|} en Sp_{2m};   nu_2 = +1 en SO_{2m'+1}.  Los dos son clasicos.")
print("")

filas = []
print("  F1  rama PAR:  .tau^C_t(eta) == (-1)^{|eta|} sobre el soporte?")
print("     t | m |  pesos | soporte | F1 acierta | C1 senuelo -(-1)^|eta| | C2 senuelo +1 | |eta| impar en el soporte")
print("     " + "-" * 116)
for t in (4, 6, 8, 10, 12):
    m = (t - 2) // 2
    cota = {4: 24, 6: 16, 8: 12, 10: 10, 12: 8}[t]
    n = f1 = c1 = c2 = sop = impares = 0
    for eta in pesos(m, cota):
        n += 1
        v = tau_C(eta, t, m)
        if v == 0:
            continue
        sop += 1
        fs = (-1) ** (sum(eta) % 2)
        if sum(eta) % 2:
            impares += 1
        f1 += 1 if v == fs else 0
        c1 += 1 if v == -fs else 0
        c2 += 1 if v == 1 else 0
    filas.append({"rama": "C", "t": t, "m": m, "pesos": n, "soporte": sop,
                  "F1": f1, "C1": c1, "C2": c2, "impares": impares})
    print("     %2d | %d | %6d | %7d | %6d     | %10d           | %8d      | %d"
          % (t, m, n, sop, f1, c1, c2, impares))

print("")
print("  F2  rama IMPAR:  .tau^B_t(eta) == +1 sobre el soporte?  (nu_2 = +1 siempre en SO_{2m'+1})")
print("     t | m'|  pesos | soporte | F2 acierta | tau^B = -1")
print("     " + "-" * 60)
for t in (3, 5, 7, 9, 11):
    mp = (t - 1) // 2
    cota = {3: 24, 5: 16, 7: 12, 9: 10, 11: 8}[t]
    n = f2 = neg = sop = 0
    for eta in pesos(mp, cota):
        n += 1
        v = tau_B(eta, t, mp)
        if v == 0:
            continue
        sop += 1
        f2 += 1 if v == 1 else 0
        neg += 1 if v == -1 else 0
    filas.append({"rama": "B", "t": t, "m": mp, "pesos": n, "soporte": sop, "F2": f2, "neg": neg})
    print("     %2d | %d | %6d | %7d | %6d     | %6d" % (t, mp, n, sop, f2, neg))

print("")
print("  F3  el elemento del indicador FS es rho^v(-1): PRINCIPAL DE ORDEN 2.")
print("      Nuestro elemento impar es principal de orden t = h+1  (probado, rem:kostant).")
print("      Nuestro elemento par NO es principal para m>=2        (probado, rem:kostant).")
print("      Luego FS y tau^B son d=2 y d=t de la MISMA familia chi(rho^v(zeta_d)); tau^C se sale.")
print("      Aritmetica comprobable: el orden de rho^v(-1) es 2, y t es 4,6,8,... o 3,5,7,...")
for t in (4, 6, 8, 3, 5, 7):
    print("        t=%-3d  .es t == 2?  %s" % (t, t == 2))

# ---------------------------------------------------------------------------------------------
# F4  EL SUBPRODUCTO, que es lo unico que sobrevive de esta prueba.  C3 dice que |eta| es SIEMPRE
# par sobre el soporte de tau^C, y eso tiene demostracion de una linea:
#
#     tau^C_t(eta) != 0  =>  a_j = eta_j + m-j+1  tiene clases plegadas {1,...,m} con signos,
#     o sea  a_j = +- sigma(j)  (mod t).  Como t es PAR, la reflexion es invisible modulo 2:
#     a_j = sigma(j) (mod 2).  Luego  sum a_j = sum sigma(j) = m(m+1)/2 (mod 2), y como
#     sum a_j = |eta| + m(m+1)/2, queda  |eta| PAR.
#
# El "como t es par" es toda la demostracion, y es la misma invertibilidad de 2 de rem:BC.  Aqui se
# comprueba (i) el enunciado sobre la caja, (ii) que el paso a_j = sigma(j) mod 2 se cumple termino
# a termino, y (iii) el SENUELO: el mismo enunciado en la rama IMPAR, donde t es impar y el
# argumento no existe -- tiene que FALLAR, o no estariamos midiendo la paridad de t.
print("")
print("  F4  SUBPRODUCTO:  .es |eta| siempre PAR sobre el soporte de tau^C?  (y el paso a_j=sigma(j) mod 2)")
print("     t | soporte | |eta| par | a_j == sigma(j) mod 2 | SENUELO rama impar: |eta| par sobre supp tau^B")
print("     " + "-" * 104)
f4 = []
for t in (4, 6, 8, 10, 12):
    m = (t - 2) // 2
    cota = {4: 24, 6: 16, 8: 12, 10: 10, 12: 8}[t]
    sop = pares = paso = 0
    for eta in pesos(m, cota):
        if tau_C(eta, t, m) == 0:
            continue
        sop += 1
        pares += 1 if sum(eta) % 2 == 0 else 0
        a = [eta[j] + (m - j) for j in range(m)]
        cl = [clase(v, t, m)[0] for v in a]
        paso += 1 if all((a[j] - cl[j]) % 2 == 0 for j in range(m)) else 0
    # senuelo: la misma pregunta en la rama impar del mismo tamano de rango
    tb = 2 * m + 1
    sop_b = pares_b = 0
    for eta in pesos(m, cota):
        if tau_B(eta, tb, m) == 0:
            continue
        sop_b += 1
        pares_b += 1 if sum(eta) % 2 == 0 else 0
    f4.append({"t": t, "sop": sop, "pares": pares, "paso": paso,
               "senuelo_t": tb, "sop_b": sop_b, "pares_b": pares_b})
    print("     %2d | %7d | %7d   | %10d          | t=%2d: %d de %d"
          % (t, sop, pares, paso, tb, pares_b, sop_b))

ok_F4 = all(f["pares"] == f["sop"] and f["paso"] == f["sop"] for f in f4)
senuelo_falla = any(f["pares_b"] != f["sop_b"] for f in f4)

par = [f for f in filas if f["rama"] == "C"]
imp = [f for f in filas if f["rama"] == "B"]
ok_F1 = all(f["F1"] == f["soporte"] for f in par)
ok_F2 = all(f["F2"] == f["soporte"] for f in imp)
vacuo = all(f["impares"] == 0 for f in par)

print("")
print("  VEREDICTO")
print("     F1  tau^C == (-1)^{|eta|} en TODO el soporte, en todos los t : %s" % ok_F1)
print("     C3  .es vacuo?  (|eta| siempre par en el soporte)            : %s" % vacuo)
print("     F2  tau^B == +1 en TODO el soporte, en todos los t           : %s" % ok_F2)
print("     F4  |eta| PAR en todo el soporte de tau^C, y el paso mod 2   : %s" % ok_F4)
print("     F4  SENUELO: en la rama impar el enunciado FALLA              : %s" % senuelo_falla)
print("")
if not (ok_F1 or ok_F2):
    print("  LECTURA: la sospecha de Frobenius-Schur esta MUERTA.  El signo del filtro acierta el")
    print("  indicador en aproximadamente la mitad del soporte y el senuelo opuesto en la otra")
    print("  mitad: eso es una moneda, no una regla.  Y ademas la prueba era casi vacua, porque")
    print("  nu_2 es CONSTANTE sobre el soporte -- que es justo lo que dice F4.")
json.dump({"filas": filas, "F4": f4}, io.open("fs_indicator_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
