# -*- coding: utf-8 -*-
# EL VERTICE DOMINANTE DEL DENOMINADOR, POR SUMAS DE MINKOWSKI.   16 de agosto de 2026.
#
# DE DONDE SALE.  mumax_odd.sage midio que el desplazamiento de la ley del peso superior es el mismo
# en las dos paridades, luego NO es el rho del factor libre (que cambia con la paridad) sino
# top Newt(N_delta).  Y eso se puede DEMOSTRAR contando, que es el teorema que la reseña marcaba
# como el mas barato de los que quedan.
#
# EL ARGUMENTO, y este guion ES el argumento, no solo su comprobacion.  El denominador es el
# Vandermonde del alfabeto,
#
#     N_delta = prod_{i<j} (x_i - x_j),   x = (1, zeta, ..., zeta^{t-1}, z_1^{+-}, ..., z_r^{+-}),
#
# y el politopo de Newton de un producto es la suma de Minkowski de los politopos de los factores
# (Ostrowski).  El maximo de un funcional lineal sobre una suma de Minkowski es la suma de los
# maximos, asi que basta elegir una direccion DOMINANTE GENERICA w = (r, r-1, ..., 1) y sumar, factor
# a factor, el exponente que la maximiza.  Cada factor tiene a lo sumo cuatro monomios, asi que la
# cuenta es O(N^2) y exacta.  Prediccion:
#
#     grado en z_1:  t (los factores zeta^a - z_1)  +  1 (el z_1 - z_1^{-1})
#                    + 2(r-1) (los cuatro factores con cada z_l, que aportan 1+1+0+0)  =  N-1,
#     y al gastar los mixtos cada variable siguiente pierde 2:
#
#            top Newt(N_delta) = (N-1, N-3, ..., N-2r+1).
#
# Ninguno de los tres sumandos es un rho.  Que se pueda ESCRIBIR como 2 rho_{C_r} + (t-1)(1..1) es
# aritmetica, y por eso vale igual en el impar, donde el factor libre es D_r.
#
# LO QUE SE MIDE
#   D1  el vertice por Minkowski contra la formula (N-1, N-3, ..., N-2r+1),  t = 2..12, r = 1..4.
#   D2  el reparto  t / 1 / 2(r-1)  del grado en z_1, contado por bloques de factores.
#
# CONTROLES
#   C0  FATAL, en los tamaños donde se puede: el vertice por Minkowski contra el vertice del
#       polinomio EXPANDIDO de verdad (determinante en el anillo de Laurent).  Solo N <= 8, porque
#       expandir N = 13 no termina -- y esa es la razon de existir de la ruta de Minkowski.
#   C1  el argmax de cada factor tiene que ser UNICO con esta w; si hubiera empate el metodo no
#       valdria y hay que avisar, no promediar.
#   C2  SEÑUELO: la formula (N-1, N-2, ..., N-r), el otro vector natural.  Tiene que fallar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage newt_denominator.sage

import json
import sys


def exponentes_factor(i, j, t, r):
    """los exponentes en (z_1..z_r) de los monomios del factor (x_i - x_j).
       las posiciones 0..t-1 son la orbita (exponente cero); despues van
       z_1, z_1^{-1}, z_2, z_2^{-1}, ..."""
    def ex(k):
        v = [0] * r
        if k >= t:
            idx = (k - t) // 2
            v[idx] = 1 if (k - t) % 2 == 0 else -1
        return tuple(v)
    return [ex(i), ex(j)]


def top_minkowski(t, r):
    """suma de los argmax por factor, en la direccion dominante generica w."""
    N = t + 2 * r
    w = [r - k for k in range(r)]          # (r, r-1, ..., 1), estrictamente decreciente y positiva
    total = [0] * r
    empates = 0
    reparto = {"orbita": 0, "auto": 0, "mixtos": 0}
    for i in range(N):
        for j in range(i + 1, N):
            ms = exponentes_factor(i, j, t, r)
            val = [sum(w[k] * m[k] for k in range(r)) for m in ms]
            # OJO: un empate solo es ambiguo si los dos EXPONENTES difieren.  Los t(t-1)/2 factores
            # internos de la orbita empatan con el mismo exponente (cero), y ahi no hay eleccion que
            # hacer: contarlos como empates era una falsa alarma de la primera version.
            if val[0] == val[1]:
                if ms[0] != ms[1]:
                    empates += 1
                mejor = ms[0]
            else:
                mejor = ms[0] if val[0] > val[1] else ms[1]
            for k in range(r):
                total[k] += mejor[k]
            # reparto del grado en z_1
            if mejor[0] != 0:
                if i < t or j < t:
                    reparto["orbita"] += mejor[0]
                elif (i - t) // 2 == (j - t) // 2:
                    reparto["auto"] += mejor[0]
                else:
                    reparto["mixtos"] += mejor[0]
    return tuple(total), empates, reparto


def top_expandido(t, r):
    """el vertice del polinomio de verdad; solo para N pequeño."""
    N = t + 2 * r
    K = CyclotomicField(t) if t > 2 else QQ
    zeta = K.gen() if t > 2 else K(-1)
    L = LaurentPolynomialRing(K, r, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(t)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(N - 1, -1, -1))
    det = matrix(L, N, N, lambda i, j: x[i] ** delta[j]).determinant()
    w = [r - k for k in range(r)]
    sop = [tuple(int(v) for v in e) if hasattr(e, '__iter__') else (int(e),)
           for e, c in zip(det.exponents(), det.coefficients()) if c != 0]
    return max(sop, key=lambda e: sum(w[k] * e[k] for k in range(r)))


print("=" * 116)
print("EL VERTICE DOMINANTE DEL DENOMINADOR, POR SUMAS DE MINKOWSKI")
print("=" * 116)
print("   t | r |  N | Minkowski          | formula (N-1,N-3,..) | ok | empates | expandido (C0)   | señuelo")
print("   " + "-" * 110)

RES = []
for t in range(2, 13):
    for r in range(1, 5):
        N = t + 2 * r
        top, emp, rep = top_minkowski(t, r)
        formula = tuple(N - 1 - 2 * k for k in range(r))
        senuelo = tuple(N - 1 - k for k in range(r))
        exp = None
        if N <= 8:
            exp = top_expandido(t, r)
        print("   %2d | %d | %2d | %-18s | %-20s | %-2s | %7d | %-16s | %s"
              % (t, r, N, str(top), str(formula), "si" if top == formula else "NO", emp,
                 str(exp) if exp else "-- (N>8)",
                 "falla" if top != senuelo else "ACIERTA (mal)"))
        sys.stdout.flush()
        RES.append({"t": int(t), "r": int(r), "N": int(N),
                    "top_minkowski": [int(v) for v in top],
                    "formula": [int(v) for v in formula],
                    "coincide": bool(top == formula),
                    "empates": int(emp),
                    "top_expandido": [int(v) for v in exp] if exp else None,
                    "C0": bool(exp is None or tuple(exp) == top),
                    "senuelo_falla": bool(top != senuelo),
                    "reparto_z1": {k: int(v) for k, v in rep.items()}})

n = len(RES)
c0 = [d for d in RES if d["top_expandido"] is not None]
nr1 = [d for d in RES if d["r"] > 1]
print("")
print("   TOTAL %d configuraciones | la formula acierta %d | sin empates ambiguos %d de %d"
      % (n, sum(1 for d in RES if d["coincide"]), sum(1 for d in RES if d["empates"] == 0), n))
print("   SEÑUELO (N-1, N-2, ...): falla en %d de %d con r >= 2.  Con r = 1 las dos formulas son"
      % (sum(1 for d in nr1 if d["senuelo_falla"]), len(nr1)))
print("   el mismo vector (N-1), asi que esos %d casos no testan nada y no se cuentan."
      % (n - len(nr1)))
print("   C0 (Minkowski == polinomio expandido) : %d de %d con N <= 8"
      % (sum(1 for d in c0 if d["C0"]), len(c0)))
print("")
print("   D2  el reparto del grado en z_1, en varias configuraciones:")
for d in RES:
    if d["r"] in (2, 3) and d["t"] in (3, 4, 6, 7):
        rp = d["reparto_z1"]
        print("       t=%2d r=%d :  orbita %2d  +  auto %d  +  mixtos %2d  =  %2d   (N-1 = %d)"
              % (d["t"], d["r"], rp["orbita"], rp["auto"], rp["mixtos"],
                 rp["orbita"] + rp["auto"] + rp["mixtos"], d["N"] - 1))
json.dump(RES, open("newt_denominator_DUMP.json", "w"), indent=1)
print("=" * 116)
print("DONE")
