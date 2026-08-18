# -*- coding: utf-8 -*-
# .POR QUE SE CANCELAN?  LA ESTRUCTURA DE LAS PROGRESIONES CON VARIOS TERMINOS.
#
# divided_differences.py dejo esto: la division por Delta_t es
#
#       c(X) = sum_{k impar >= 1} nu~(X + t k),
#
# y de los 1349 X con algun termino no nulo, 1016 tienen UNO (y ahi c = +-1) y 333 tienen DOS o
# CUATRO, que suman exactamente 0.  Nunca tres.  Nunca una cancelacion parcial.
#
# 1, 2, 4 son potencias de 2.  La hipotesis obvia es que el conjunto de k que contribuyen es un
# PRODUCTO de elecciones binarias independientes -- una caja combinatoria -- y que cada eleccion
# invierte el signo de nu.  Si eso es cierto, la suma es 0 en cuanto hay una eleccion, y (L1) queda
# probada: no hace falta disjuncion, hace falta una involucion libre que invierta el signo.
#
# LO QUE SE FALSA, en este orden:
#   P1  el numero de terminos es siempre una potencia de 2.
#   P2  el conjunto de k contribuyentes ES el producto cartesiano de sus proyecciones por
#       coordenada.  (Test duro: un conjunto de tamano 4 puede no ser una caja.)
#   P3  cada proyeccion tiene a lo sumo 2 valores.
#   P4  el signo es MULTIPLICATIVO sobre la caja: nu~(X+tk) = s0 * prod_j sigma_j(k_j) con
#       sigma_j(.) en {+-1} tomando los dos valores.  Equivale a que cada eleccion binaria invierta
#       el signo, que es lo que hace que la suma sea 0.
#   P5  y los k que contribuyen difieren en coordenadas por multiplos PARES de 2 (o sea, k-k' par),
#       que es lo que dice que estan en la misma progresion y no en otra.
#
# SENUELO
#   D1  la misma pregunta sobre cajas ALEATORIAS de nu (barajando los signos de nu): P4 tiene que
#       hundirse.  Si P4 saliera igual con signos barajados, no estaria midiendo nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_cancelacion.py

import itertools
import json
import random
from collections import Counter

from divided_differences import (CASOS, nu_de, nu_extendida, enderezar_D)


def contribuyentes(nu, t, r):
    """Para cada X dominante regular, los k impares con nu~(X+tk) != 0, y sus valores."""
    if not nu:
        return {}
    M = max(max(abs(v) for v in k) for k in nu)
    out = {}
    for X in itertools.product(range(-M, M + 1), repeat=r):
        if any(X[j] <= X[j + 1] for j in range(r - 2)):
            continue
        if not (X[r - 2] > abs(X[r - 1])):
            continue
        ks = []
        for j in range(r):
            L, k = [], 1
            while X[j] + t * k <= M:
                if X[j] + t * k >= -M:
                    L.append(k)
                k += 2
            ks.append(L)
        if any(not L for L in ks):
            continue
        viv = {}
        for k in itertools.product(*ks):
            v = nu_extendida(nu, tuple(X[j] + t * k[j] for j in range(r)))
            if v:
                viv[k] = v
        if viv:
            out[X] = viv
    return out


def es_caja(ks):
    """.El conjunto de tuplas ks es el producto cartesiano de sus proyecciones?"""
    r = len(next(iter(ks)))
    proy = [sorted({k[j] for k in ks}) for j in range(r)]
    return set(itertools.product(*proy)) == set(ks), proy


def signo_multiplicativo(viv):
    """.Existe s0 y sigma_j con nu(k) = s0 prod_j sigma_j(k_j), y cada sigma_j no constante donde
    la proyeccion tiene 2 valores?  Se resuelve fijando el k mas pequeno como referencia."""
    ks = sorted(viv)
    r = len(ks[0])
    base = ks[0]
    s0 = viv[base]
    sigma = []
    for j in range(r):
        vals = sorted({k[j] for k in ks})
        if len(vals) == 1:
            sigma.append({vals[0]: 1})
            continue
        if len(vals) > 2:
            return False, "proyeccion con %d valores" % len(vals)
        # el signo de cambiar SOLO la coordenada j respecto de base
        otro = [k for k in ks if k[j] != base[j] and all(k[i] == base[i] for i in range(r) if i != j)]
        if not otro:
            return False, "la caja no tiene el vecino en j=%d" % j
        s = viv[otro[0]] // s0
        if s not in (1, -1):
            return False, "razon no unitaria"
        sigma.append({base[j]: 1, otro[0][j]: s})
    for k in ks:
        pred = s0
        for j in range(r):
            if k[j] not in sigma[j]:
                return False, "coordenada fuera de la caja"
            pred *= sigma[j][k[j]]
        if pred != viv[k]:
            return False, "no multiplicativo"
    # y cada eleccion binaria tiene que INVERTIR el signo, no repetirlo
    for j in range(r):
        if len(sigma[j]) == 2 and set(sigma[j].values()) != {1, -1}:
            return False, "una eleccion binaria NO invierte el signo"
    return True, "ok"


print("=" * 104)
print("POR QUE SE CANCELAN: .son cajas con signo multiplicativo?")
print("=" * 104)
print("")

P1 = P1n = P2 = P2n = P3 = P3n = P4 = P4n = P5 = P5n = 0
D1 = D1n = 0
tam = Counter()
fallos = []
rnd = random.Random(20260816)

for (t, r, cota) in CASOS:
    Rp = (t - 1) // 2 + r
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu = nu_de(list(Lam), t, r)
        if not nu:
            continue
        con = contribuyentes(nu, t, r)
        # el senuelo: los MISMOS soportes con los signos barajados
        nu_falso = dict(nu)
        claves = list(nu_falso)
        for kk in claves:
            nu_falso[kk] = rnd.choice([1, -1])
        con_falso = contribuyentes(nu_falso, t, r)

        for X, viv in con.items():
            n = len(viv)
            tam[n] += 1
            P1n += 1
            if n & (n - 1) == 0:
                P1 += 1
            if n == 1:
                continue
            ok, proy = es_caja(set(viv))
            P2n += 1
            P2 += 1 if ok else 0
            P3n += 1
            P3 += 1 if all(len(p) <= 2 for p in proy) else 0
            P5n += 1
            P5 += 1 if all((a - b) % 2 == 0 for k1 in viv for k2 in viv
                           for a, b in zip(k1, k2)) else 0
            if not ok:
                if len(fallos) < 3:
                    fallos.append({"que": "no es caja", "t": t, "r": r, "Lambda": list(Lam),
                                   "X": list(X), "ks": [list(k) for k in viv]})
                continue
            good, motivo = signo_multiplicativo(viv)
            P4n += 1
            P4 += 1 if good else 0
            if not good and len(fallos) < 3:
                fallos.append({"que": "signo: " + motivo, "t": t, "r": r, "Lambda": list(Lam),
                               "X": list(X), "viv": {str(k): v for k, v in viv.items()}})

        for X, viv in con_falso.items():
            if len(viv) < 2:
                continue
            ok, _ = es_caja(set(viv))
            if not ok:
                continue
            good, _ = signo_multiplicativo(viv)
            D1n += 1
            D1 += 1 if good else 0

print("  reparto del numero de terminos por progresion : %s" % dict(sorted(tam.items())))
print("")
print("  P1  el numero de terminos es potencia de 2      : %d de %d" % (P1, P1n))
print("  P2  el conjunto de k es una CAJA (producto)     : %d de %d" % (P2, P2n))
print("  P3  cada proyeccion tiene a lo sumo 2 valores   : %d de %d" % (P3, P3n))
print("  P5  los k difieren en pares (misma progresion)  : %d de %d" % (P5, P5n))
print("  P4  el signo es MULTIPLICATIVO y cada eleccion")
print("      binaria lo INVIERTE  (esto es la cancelacion): %d de %d" % (P4, P4n))
print("")
print("  D1  SENUELO: lo mismo con los signos de nu barajados : %d de %d" % (D1, D1n))
if fallos:
    print("")
    print("  !! primeros fallos:")
    for f in fallos:
        print("    " + json.dumps(f)[:260])
print("")
print("  LECTURA: P2+P4 juntos SON la cancelacion.  Si salen enteros, cada progresion con mas de un")
print("  termino es un producto de elecciones binarias que invierten el signo, luego suma 0, y c(X)")
print("  vale +-1 exactamente cuando la progresion corta el soporte una sola vez.  Eso es (L1).")

json.dump({"tam": {str(k): v for k, v in sorted(tam.items())},
           "P1": [P1, P1n], "P2": [P2, P2n], "P3": [P3, P3n], "P4": [P4, P4n], "P5": [P5, P5n],
           "D1_senuelo": [D1, D1n], "fallos": fallos},
          open("_probe_cancelacion_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
