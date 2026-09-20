# -*- coding: utf-8 -*-
r"""reticulo_orbitas.py -- el reticulo de orbitas de los semisistemas de suma nula.

DE DONDE VIENE.  La conjetura "una sola orbita" murio en la auditoria del 20-sep (m = 19 tiene tres,
y N(m) > phi(m) para todo m impar >= 19).  Lo que queda vivo, y no lo da la formula, es QUE orbitas
hay: cuantas, de que tamanos y con que estabilizadores.

LA OBSERVACION QUE ORGANIZA EL CALCULO.  Si H = {a : aT = T} es el estabilizador de un semisistema,
entonces -1 NO esta en H --- si lo estuviera, T contendria a la vez t y -t --- luego **H tiene orden
IMPAR**.  Por tanto el estabilizador divide la parte impar de phi(m), y las orbitas tienen tamano
phi(m)/h con h impar.  Y ese h es exactamente el mismo h del censo de P3a: el estabilizador
multiplicativo del espectro.  Las dos preguntas son la misma.

QUE MIDE ESTE GUION.  Para cada m impar del rango: el reparto de estabilizadores, el numero de
orbitas de cada tamano, y si TODO divisor impar de phi(m) aparece como estabilizador (la conjetura
que sustituye a la muerta).

    python reticulo_orbitas.py [MMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import sys
from itertools import product
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 31


def unidades(m):
    return [a for a in range(1, m) if gcd(a, m) == 1]


def parte_impar(n):
    while n % 2 == 0:
        n //= 2
    return n


print("%-5s %-8s %-8s %-10s %-26s %s"
      % ("m", "N(m)", "phi(m)", "orbitas", "estabilizadores h: cuantas", "divisores impares de phi"))
for m in range(7, MMAX + 1, 2):
    n = (m - 1) // 2
    U = unidades(m)
    sols = []
    for eleccion in product(*[(t, m - t) for t in range(1, n + 1)]):
        T = frozenset(eleccion)
        if sum(T) % m == 0:
            sols.append(T)
    vistos, orbitas = set(), []
    for T in sols:
        if T in vistos:
            continue
        orb = {frozenset((a * t) % m for t in T) for a in U}
        orb = {O for O in orb if O in set(sols)} | {T}
        vistos |= orb
        h = sum(1 for a in U if frozenset((a * t) % m for t in T) == T)
        orbitas.append((len(orb), h))
    rep = {}
    for tam, h in orbitas:
        rep[h] = rep.get(h, 0) + 1
    phi_m = len(U)
    div_impares = [d for d in range(1, phi_m + 1) if phi_m % d == 0 and d % 2 == 1]
    faltan = [d for d in div_impares if d not in rep]
    print("%-5d %-8d %-8d %-10d %-26s %s%s"
          % (m, len(sols), phi_m, len(orbitas),
             " ".join("h=%d:%d" % (h, rep[h]) for h in sorted(rep)),
             div_impares, "   FALTAN %s" % faltan if faltan else "   (todos aparecen)"))
    sys.stdout.flush()

print("")
print("Si algun h par apareciera, la observacion de arriba seria falsa; no deberia ocurrir.")
