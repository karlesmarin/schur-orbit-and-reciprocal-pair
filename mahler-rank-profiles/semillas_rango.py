# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: control de dos corolarios del paper III frente a las filas medidas de fig_filas_datos.json
#      (chi = (-75619/.), p = 5, a = 4; filas dependientes calculadas por eliminacion, 624 filas positivas):
#        (1) cor:rank   rank B = [mu(Z_p) != 0] + sum_b sum_{t<a} max(p^t - r_b, 0);
#        (2) cor:D      semillas j_{b,i} = p i + [b-i]_{p-1} (j_{0,0} = p-1): fila dependiente <=> i < r_b,
#                       para b par y j_{b,i} < p^a (a nivel de filas de B; el corolario del paper las usa en j < P - p^{a-1}).
#      Control negativo: con r_0 cambiado a 5 o 7 la prediccion de semillas debe fallar.
# POR QUE: las dos formulas son corolarios del Teorema A en el paper (cor:rank, cor:D); este es su control.
# FECHA: 30-sep-2026
import json
import os

d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_filas_datos.json'), encoding='utf-8'))
p, a = d['p'], d['a']
P = p ** a
dep = set(d['ausentes'])
INF = 10 ** 9
r = [x if x is not None else INF for x in d['r_b']]

# La fila 0 no esta en 'ausentes' (el fichero lista filas positivas). Se decide aparte: fila 0 nueva <=>
# mu_chi(Z_p) = 2(1 - chi(2)) B_{1,chi} != 0 mod p, con B_{1,chi} = (1/n) sum a chi(a), chi = (-n/.) (Kronecker).
from sympy.functions.combinatorial.numbers import jacobi_symbol  # noqa: E402

n = d['n']
assert n % 4 == 3  # (-n/.) = (./n) para n primo = 3 mod 4
chi = lambda x: jacobi_symbol(x % n, n) if x % n else 0  # noqa: E731
suma = sum(x * chi(x) for x in range(1, n))
masa = 2 * (1 - chi(2)) * suma * pow(n, -1, p) % p  # 2(1-chi(2)) B_{1,chi} mod p
fila0_nueva = masa != 0
rango = (P - 1 - len(dep)) + fila0_nueva
formula = fila0_nueva + sum(max(p ** t - r[b], 0) for b in range(p - 1) for t in range(a))
print(f"(p, n, a) = ({p}, {n}, {a}); r_b = {d['r_b']}")
print(f"(0) mu_chi(Z_p) = 2(1-chi(2))B_1,chi = {masa} mod {p}: fila 0 {'nueva' if fila0_nueva else 'dependiente'}")
print(f"(1) rango medido {rango}, formula {formula}: {'OK' if rango == formula else 'FALLA'}")
for r0 in (5, 7):
    f2 = fila0_nueva + sum(max(p ** t - ([r0] + r[1:])[b], 0) for b in range(p - 1) for t in range(a))
    print(f"    control negativo r_0 = {r0}: formula {f2} (debe ser != {rango}): {'OK' if f2 != rango else 'FALLA'}")


def semillas(rr):
    malas = tot = 0
    for b in range(0, p - 1, 2):
        i = 0
        while True:
            j = (p - 1) if (b == 0 and i == 0) else p * i + ((b - i) % (p - 1))
            if j >= P:
                break
            tot += 1
            malas += ((j in dep) != (i < rr[b]))
            i += 1
    return tot, malas


tot, malas = semillas(r)
print(f"(2) semillas b par, j < p^a: {tot} comparadas, {malas} discrepancias: {'OK' if malas == 0 else 'FALLA'}")
for r0 in (5, 7):
    t2, m2 = semillas([r0] + r[1:])
    print(f"    control negativo r_0 = {r0}: {m2} discrepancias (debe ser > 0): {'OK' if m2 > 0 else 'FALLA'}")
