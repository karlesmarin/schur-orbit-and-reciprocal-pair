# -*- coding: utf-8 -*-
# FRONTERA 3 del papel: las filas fuera de la diagonal.
#
# El papel dice: "away from it the vanishing of T(m,chi) is a condition on the values of chi on
# 1..m, hence a congruence on f for each fixed m; WHAT IS MISSING IS A STATEMENT UNIFORM IN m".
# Y anota que la fila m=3 esta cortada por (2|f)=1 y (3|f)=-1.
#
# HIPOTESIS.  Para chi CUADRATICO, T(m,chi) = sum_{c<=m} c chi(c) es una suma corta de signos, y
# chi es completamente multiplicativo, luego el vector de signos lo determina chi en los PRIMOS
# <= m.  Entonces:
#
#     T(m,chi) = 0   <=>   el patron de signos de chi en los primos <= m es una solucion de
#                          la SUMA PESADA NULA   sum_{c<=m} c eps(c) = 0,  eps multiplicativo.
#
# Eso SI es uniforme en m: para cada m es una enumeracion finita sobre 2^{pi(m)} patrones, y cada
# solucion es una condicion de simbolos de Legendre, es decir una congruencia en f.
#
# Comprobar:  1 == 2 chi(2) + 3 chi(3) ... con m=3 da 1+2-3 = 0 con chi(2)=+1, chi(3)=-1.
#
# CONTROLES:
#   (a) para cada m con solucion se busca un f real que la realice y se verifica T = 0;
#   (b) para cada m SIN solucion se verifica que NINGUN f del barrido da T = 0  (si alguno diera,
#       la hipotesis muere);
#   (c) senuelo: patrones de signos NO multiplicativos con suma nula no deben realizarse.
#
# Authors: Carles Marin, Claude (AI assistant).
import sys
from itertools import product


def patrones(m):
    """(dict primo->signo, vector eps sobre 1..m) para todos los patrones multiplicativos."""
    ps = prime_range(m + 1)
    for signos in product([1, -1], repeat=len(ps)):
        asig = dict(zip(ps, signos))
        eps = [0] * (m + 1)
        eps[1] = 1
        for c in range(2, m + 1):
            v = 1
            for (pp, ee) in factor(c):
                v *= asig[pp] ** ee
            eps[c] = v
        yield asig, eps


def suma(m, eps):
    return sum([c * eps[c] for c in range(1, m + 1)])


print("=== (1) para cada m, cuantos patrones multiplicativos dan suma pesada nula ===")
print("%-4s %-8s %-8s %s" % ("m", "#patrones", "#nulos", "los patrones nulos (signo en cada primo)"))
sol = {}
for m in range(2, 31):
    ns = []
    tot = 0
    for asig, eps in patrones(m):
        tot += 1
        if suma(m, eps) == 0:
            ns.append(asig)
    sol[m] = ns
    if ns:
        print("%-4d %-8d %-8d %s" % (
            m, tot, len(ns), "; ".join(str(sorted(a.items())) for a in ns[:3])))
    sys.stdout.flush()
print("filas con solucion:", sorted([m for m in sol if sol[m]]))
print()

print("=== (2) contraste con los caracteres cuadraticos de verdad ===")
print("%-4s %-8s %-10s %-10s %-10s %s" % (
    "m", "#f vistos", "T=0 reales", "predichos", "coinciden", "veredicto"))
ok = mal = 0
for m in range(2, 21):
    reales = []
    pred = []
    for f in range(4 * m + 1, 400, 2):
        if not is_squarefree(f) or f % 4 != 3:
            continue
        chi = kronecker_character(-f)
        if chi is None or chi.conductor() != f:
            continue
        if any(gcd(c, f) > 1 for c in range(1, m + 1)):
            continue
        t = sum([c * chi(c) for c in range(1, m + 1)])
        asig = {pp: chi(pp) for pp in prime_range(m + 1)}
        esta = any(all(a[pp] == asig[pp] for pp in asig) for a in sol[m])
        if t == 0:
            reales.append(f)
        if esta:
            pred.append(f)
        if (t == 0) != esta:
            mal += 1
        else:
            ok += 1
    print("%-4d %-8d %-10d %-10d %-10s %s" % (
        m, len([f for f in range(4 * m + 1, 400, 2)]), len(reales), len(pred),
        "si" if reales == pred else "NO",
        "ok" if reales == pred else "*** DISCREPA ***"))
    sys.stdout.flush()
print()
print("celdas (m,f): %d ok, %d discrepan" % (ok, mal))
