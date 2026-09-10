# -*- coding: utf-8 -*-
# EL SISTEMA R_P(q) EN TODOS LOS TIPOS, Y SU RELACION EXACTA CON EL R(m) DE NPP-POLO.
# 25 de agosto de 2026.
#
# POR QUE.  Polo (arXiv:2504.09204) clasifica, a peticion de Prasad, el sistema
#     R(m) = { beta : m | ht(beta) },
# que es el del centralizador de las potencias del elemento PRINCIPAL a_Q.  El segundo elemento de
# Kostant, a_P, da OTRO sistema, y en tipo C lo hemos medido:  R_P(q) = { alpha : q | <alpha,rho> }.
# Esta gate hace las dos cosas que faltan:
#   (1) definir R_P(q) en TODOS los tipos, y
#   (2) comprobar la relacion exacta con R(q) en el grupo DUAL.
#
# LA DEFINICION UNIFORME, deducida de la normalizacion de Kostant.  a_P = exp(x_P) con
# <alpha, x_P> proporcional a (alpha, alpha).  Luego, sobre las simples, el exponente vale 1 en las
# CORTAS y r en las LARGAS (r = razon de longitudes).  Para alpha = sum c_i alpha_i:
#
#     ht_P(alpha) := sum_i c_i * w_i ,     w_i = 1 si alpha_i corta,  r si larga
#     R_P(q)      := { alpha : q | ht_P(alpha) }
#
# y a_Q da la altura ordinaria ht(alpha) = sum_i c_i, o sea R(q) de NPP-Polo.
#
# LA PREDICCION, que sale de una linea de algebra y NO se postula:
#     alpha = sum c_i alpha_i  =>  ht_P(alpha) = ((alpha,alpha)/2) * ht(alpha^v)
# porque alpha^v = sum_i c_i ((alpha_i,alpha_i)/(alpha,alpha)) alpha_i^v.  De ahi, leyendo por
# corraices,
#     alpha in R_P(q)  <=>  q | ( (alpha,alpha)/2 ) * ht(alpha^v)
# o sea: para las CORTAS de G (largas del dual) la condicion es la altura pelada, y para las LARGAS
# de G la altura MULTIPLICADA POR r.  Luego
#
#     R_P(q)^v = R(q) del dual   <=>   gcd(q, r) = 1,
#
# y si no, R_P(q)^v es ESTRICTAMENTE MAYOR, y lo que sobra son exactamente las corraices beta con
#     r*ht(beta) = 0 mod q   pero   ht(beta) != 0 mod q.
#
# QUE SE MIDE
#   C1  CONTROL: el generador de raices produce el numero correcto de raices positivas en cada tipo.
#   C2  CONTROL: en los tipos SIMPLEMENTE ENLAZADOS (r=1) los dos sistemas deben coincidir SIEMPRE.
#       Si no coincidieran, el instrumento estaria roto.
#   V1  la identidad  ht_P(alpha) = ((alpha,alpha)/2) * ht(alpha^v),  evaluada raiz a raiz.
#   V2  la equivalencia  [R_P(q)^v = R(q) dual]  <=>  gcd(q,r) = 1,  en B, C, F4, G2, todo q.
#   V3  la caracterizacion de lo que sobra.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python R_P_todos_los_tipos.py > R_P_todos_los_tipos_OUT.txt 2>&1

import sys
from fractions import Fraction
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 96


# ---------------------------------------------------------------- matrices de Cartan
# A[i][j] = <alpha_j, alpha_i^v>.  Convenio: la fila i dice como actua alpha_i^v.

def cartan(tipo, n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    if tipo in ("A", "B", "C", "F", "G"):
        for i in range(n - 1):
            A[i][i + 1] = A[i + 1][i] = -1
    if tipo == "B" and n >= 2:        # alpha_n corta;  <alpha_n, alpha_{n-1}^v> = -1, <alpha_{n-1}, alpha_n^v> = -2
        A[n - 1][n - 2] = -2
    if tipo == "C" and n >= 2:        # alpha_n larga
        A[n - 2][n - 1] = -2
    if tipo == "D":
        for i in range(n - 2):
            A[i][i + 1] = A[i + 1][i] = -1
        A[n - 3][n - 1] = A[n - 1][n - 3] = -1
    if tipo == "F":                   # F4: alpha_1,alpha_2 largas; alpha_3,alpha_4 cortas
        A = [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -2, 2, -1], [0, 0, -1, 2]]
    if tipo == "G":                   # G2: alpha_1 corta (d=1), alpha_2 larga (d=3).
        # Convenio de esta gate: A[i][j] = <alpha_j, alpha_i^v> = (alpha_i,alpha_j)/d_i.
        # Luego d_0*A[0][1] = d_1*A[1][0]  =>  1*A[0][1] = 3*(-1)  =>  A[0][1] = -3.
        # ⚠ La version del 25-ago tenia [[2,-1],[-3,2]], que es la TRANSPUESTA, y el control del
        #   numero de raices no podia verlo: la transpuesta tambien da 6 positivas.  La caza C2.
        A = [[2, -3], [-1, 2]]
    return A


def longitudes(tipo, n):
    """d_i = (alpha_i, alpha_i)/2,  normalizado a 1 en las CORTAS."""
    if tipo in ("A", "D"):
        return [1] * n
    if tipo == "B":                   # alpha_1..alpha_{n-1} largas, alpha_n corta
        return [2] * (n - 1) + [1]
    if tipo == "C":                   # alpha_1..alpha_{n-1} cortas, alpha_n larga
        return [1] * (n - 1) + [2]
    if tipo == "F":
        return [2, 2, 1, 1]
    if tipo == "G":
        return [1, 3]
    raise ValueError(tipo)


def razon(tipo):
    return {"A": 1, "D": 1, "B": 2, "C": 2, "F": 2, "G": 3}[tipo]


def raices_positivas(A, n):
    """Genera las raices positivas por cuerdas, en coordenadas de raices simples."""
    simples = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
    R = set(simples)
    frontera = list(simples)
    while frontera:
        nueva = []
        for a in frontera:
            for i in range(n):
                # p = max{k : a - k alpha_i es raiz}
                p = 0
                b = list(a)
                while True:
                    b[i] -= 1
                    if min(b) < 0 or tuple(b) not in R:
                        break
                    p += 1
                pair = sum(a[j] * A[i][j] for j in range(n))     # <a, alpha_i^v>
                q = p - pair
                if q >= 1:
                    c = tuple(a[j] + (1 if j == i else 0) for j in range(n))
                    if c not in R:
                        R.add(c)
                        nueva.append(c)
        frontera = nueva
    return sorted(R, key=lambda c: (sum(c), c))


def norma(c, A, d):
    """(alpha, alpha)/2  con (alpha_i,alpha_j) = d_i A[i][j]... simetrizado: d_i A[i][j] = d_j A[j][i]."""
    n = len(c)
    s = Fraction(0)
    for i in range(n):
        for j in range(n):
            s += Fraction(c[i] * c[j] * d[i] * A[i][j], 2)
    return s          # = (alpha,alpha)/2


CASOS = [("A", 4, 10), ("D", 4, 12), ("D", 5, 20),
         ("B", 2, 4), ("B", 3, 9), ("B", 4, 16), ("B", 5, 25),
         ("C", 2, 4), ("C", 3, 9), ("C", 4, 16), ("C", 5, 25),
         ("F", 4, 24), ("G", 2, 6)]

print(SEP)
print("C1 -- CONTROL: numero de raices positivas por tipo")
sistemas = {}
malos = 0
for tipo, n, esperado in CASOS:
    A = cartan(tipo, n)
    d = longitudes(tipo, n)
    R = raices_positivas(A, n)
    sistemas[(tipo, n)] = (A, d, R)
    ok = len(R) == esperado
    malos += (not ok)
    print(f"     {tipo}{n}: {len(R):>3} positivas, esperadas {esperado:>3}   {'ok' if ok else 'FALLA'}")
print(f"     {len(CASOS)} tipos, {malos} fallos"
      + ("  <-- el generador es correcto" if malos == 0 else "  <-- GENERADOR ROTO, nada de abajo vale"))

print(SEP)
print("C2 -- CONTROL QUE FALTABA: la forma simetrizada  d_i * A[i][j] = d_j * A[j][i].")
print("      Es el que caza una matriz de Cartan TRANSPUESTA, cosa que el conteo de raices no ve.")
malos = 0
for (tipo, n), (A, d, R) in sistemas.items():
    for i in range(n):
        for j in range(n):
            if d[i] * A[i][j] != d[j] * A[j][i]:
                malos += 1
                print(f"      FALLA {tipo}{n}: d[{i}]A[{i}][{j}]={d[i]*A[i][j]} != d[{j}]A[{j}][{i}]={d[j]*A[j][i]}")
print(f"      {len(sistemas)} tipos, {malos} fallos"
      + ("  <-- las matrices son coherentes con las longitudes" if malos == 0 else "  <-- MATRIZ ROTA"))

print(SEP)
print("C3 -- CONTROL: cuantas raices positivas CORTAS y LARGAS hay de verdad, contra lo conocido")
print("      (B_n: n cortas;  C_n: n largas;  F4: 12 y 12;  G2: 3 y 3;  A,D: todas iguales)")
ESPERADO = {("B", 2): (2, 2), ("B", 3): (3, 6), ("B", 4): (4, 12), ("B", 5): (5, 20),
            ("C", 2): (2, 2), ("C", 3): (6, 3), ("C", 4): (12, 4), ("C", 5): (20, 5),
            ("F", 4): (12, 12), ("G", 2): (3, 3)}
malos = 0
for (tipo, n), (A, d, R) in sistemas.items():
    normas = [norma(c, A, d) for c in R]
    mn = min(normas)
    cortas = sum(1 for x in normas if x == mn)
    largas = len(R) - cortas
    if (tipo, n) in ESPERADO:
        e = ESPERADO[(tipo, n)]
        ok = (cortas, largas) == e
        malos += (not ok)
        print(f"      {tipo}{n}: cortas={cortas:>2} largas={largas:>2}   esperado {e}   {'ok' if ok else 'FALLA'}")
    else:
        print(f"      {tipo}{n}: cortas={cortas:>2} largas={largas:>2}   (simplemente enlazado)")
print(f"      {malos} fallos")

print(SEP)
print("V1 -- ht_P(alpha) = ((alpha,alpha)/2) * ht(alpha^v).")
print("      ⚠ NO es una comprobacion: ht(alpha^v) se CALCULA dividiendo, asi que la igualdad es")
print("      una tautologia.  Se deja escrita como DEFINICION, no como control.  Lo que si se")
print("      comprueba es que ht(alpha^v) sale ENTERO en las 180 raices --- eso no es gratis, y")
print("      es lo que dice que la formula c_i^v = c_i d_i / nrm es la correcta.")
tot = malos = 0
for (tipo, n), (A, d, R) in sistemas.items():
    for c in R:
        nrm = norma(c, A, d)
        htdual = sum(Fraction(c[i] * d[i], 1) for i in range(n)) / nrm
        tot += 1
        if htdual.denominator != 1:
            malos += 1
            if malos <= 3:
                print(f"      NO ENTERO {tipo}{n} c={c} ht(alpha^v)={htdual}")
print(f"      {tot} raices, {malos} no enteros")

print(SEP)
print("V2 -- LOS TRES ENUNCIADOS, medidos por separado.  (La primera version de esta gate afirmaba")
print("      'iguales <=> gcd(q,r)=1' y daba 32 discrepancias: 26 porque los dos conjuntos estaban")
print("      VACIOS --- iguales por vacuidad --- y 6 por la matriz de G2 transpuesta.  Lo correcto")
print("      es una INCLUSION mas una caracterizacion de lo que sobra, y la coprimalidad solo da")
print("      una implicacion.)")
print()
print("      E1  INCLUSION:  R_P(q)^v contiene siempre a R(q) del dual.")
print("      E2  LO QUE SOBRA es exactamente { alpha LARGA de G : q | r*ht(alpha^v), q no| ht(alpha^v) }.")
print("      E3  gcd(q,r) = 1  =>  iguales.   (El reciproco es FALSO: fallan por vacuidad.)")
print()
e1 = e2 = e3 = casos = 0
f1 = f2 = f3 = 0
por_tipo = {}
for (tipo, n), (A, d, R) in sistemas.items():
    r = razon(tipo)
    h = max(sum(c) for c in R) + 1
    datos = []
    for c in R:
        nrm = norma(c, A, d)
        htP = sum(c[i] * d[i] for i in range(n))
        htdual = int(Fraction(htP) / nrm)
        datos.append((c, nrm, htP, htdual))
    mn = min(x[1] for x in datos)
    for q in range(2, 2 * h + 2):
        RP = {c for (c, nrm, htP, htd) in datos if htP % q == 0}
        RQ = {c for (c, nrm, htP, htd) in datos if htd % q == 0}
        sobra_obs = RP - RQ
        sobra_pred = {c for (c, nrm, htP, htd) in datos
                      if nrm != mn and (r * htd) % q == 0 and htd % q != 0}
        casos += 1
        if RQ <= RP:
            e1 += 1
        else:
            f1 += 1
            if f1 <= 4:
                print(f"      E1 FALLA {tipo}{n} q={q}: R(q) no esta dentro de R_P(q)")
        if sobra_obs == sobra_pred:
            e2 += 1
        else:
            f2 += 1
            if f2 <= 4:
                print(f"      E2 FALLA {tipo}{n} q={q}: sobra {sorted(sobra_obs)} vs predicho {sorted(sobra_pred)}")
        if gcd(q, r) != 1 or RP == RQ:
            e3 += 1
        else:
            f3 += 1
            if f3 <= 4:
                print(f"      E3 FALLA {tipo}{n} q={q}: coprimos pero distintos")
        if RP != RQ:
            por_tipo.setdefault(f"{tipo}{n}", []).append((q, len(sobra_obs)))
print(f"      {casos} pares (tipo, q):   E1 {e1}/{casos}   E2 {e2}/{casos}   E3 {e3}/{casos}")
print()
print("      Donde los dos sistemas DIFIEREN de verdad, por tipo (q: cuantas raices sobran):")
for k in sorted(por_tipo):
    print(f"        {k:>4}: {', '.join(f'{q}:{s}' for q, s in por_tipo[k])}")

print(SEP)
print("LECTURA.")
print("  1. En los simplemente enlazados los dos sistemas son el MISMO: rho = rho^v.  La pregunta")
print("     solo tiene contenido en B, C, F4 y G2.")
print("  2. R_P(q)^v CONTIENE siempre a R(q) del dual, y lo que sobra son siempre raices LARGAS de")
print("     G, con una caracterizacion explicita.  Coincidencia garantizada si gcd(q,r)=1; el")
print("     reciproco no vale, porque para q grande los dos conjuntos se vacian.")
print("  3. Luego la clasificacion de Polo DA la del lado a_P: R_P(q)^v es su R(q) mas un conjunto")
print("     explicito de raices largas.  La pregunta que Polo deja abierta tiene gemela, y la")
print("     gemela se contesta con la suya mas un termino que se escribe en una linea.")
print(SEP)
