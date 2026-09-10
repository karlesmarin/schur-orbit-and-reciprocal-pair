# -*- coding: utf-8 -*-
"""La clasificacion leida en el alfabeto, bien esta vez.

A = {zeta^{+-j} : j=1..m} como multiconjunto de mu_q.  Para cada multiplicidad c se mide
el coste sum_r |A[r] - c|.  Se toma el c que lo minimiza.

REGLA PROPUESTA:  (A) se cumple  <=>  A es un multiplo de mu_q MENOS a lo sumo dos letras,
                  y las letras que faltan estan en los PUNTOS FIJOS de w -> w^{-1}
                  (r=0, que es la letra 1;  y r=q/2 con q par, que es la letra -1),
                  sin ninguna letra SOBRANTE.       [ mas la esporadica q=6, v=1 ]

Eso es exactamente la Proposicion 3.1 del articulo, que hoy solo se enuncia para las potencias.
"""
import sys
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def alfabeto(m, q):
    c = Counter()
    for j in range(1, m + 1):
        c[j % q] += 1
        c[(-j) % q] += 1
    return [c.get(r, 0) for r in range(q)]


def mejor_c(A, q):
    """el c de coste minimo; EN CASO DE EMPATE, el MAYOR, para que la descripcion sea
    'multiplo de mu_q MENOS letras' y no 'menos letras MAS otras'.  Con q=2 y m impar los
    tres c consecutivos empatan, y solo el mayor da la descripcion pura."""
    best, bc = None, None
    for c in range(0, max(A) + 2):
        cost = sum(abs(x - c) for x in A)
        if best is None or cost < best or (cost == best and c > bc):
            best, bc = cost, c
    return bc, best


def regla(m, q):
    A = alfabeto(m, q)
    c, cost = mejor_c(A, q)
    fijos = {0} | ({q // 2} if q % 2 == 0 else set())
    sobra = [r for r in range(q) if A[r] > c]
    falta = [r for r in range(q) for _ in range(c - A[r]) if A[r] < c]
    ok = (not sobra) and len(falta) <= 2 and all(r in fijos for r in falta)
    return ok, c, falta, sobra


def prediccion(m, q):
    a = (q - 1) // 2
    u, v = divmod(m, q)
    return (v in (0, a, q - 1)) or (q == 6 and v == 1)


print("=" * 100)
print("REGLA DEL ALFABETO  contra  LA CLASIFICACION")
print("=" * 100)
ac = des = 0
malos = []
for q in range(2, 15):
    for m in range(1, 40):
        r, c, falta, sobra = regla(m, q)
        r = r or (q == 6 and m % 6 == 1)
        p = prediccion(m, q)
        if r == p:
            ac += 1
        else:
            des += 1
            malos.append((m, q, c, falta, sobra, p, r))
print("  ACUERDOS %d   DESACUERDOS %d" % (ac, des))
for x in malos[:12]:
    print("    m=%-3d q=%-3d  c=%-3d falta=%-14s sobra=%-14s pred=%-5s regla=%-5s" % x)

print()
print("=" * 100)
print("LAS CUATRO FAMILIAS, EN LENGUAJE DE ALFABETO")
print("=" * 100)
vistos = {}
for q in range(2, 13):
    for m in range(1, 60):
        if not prediccion(m, q):
            continue
        u, v = divmod(m, q)
        a = (q - 1) // 2
        A = alfabeto(m, q)
        c, cost = mejor_c(A, q)
        fijos = {0} | ({q // 2} if q % 2 == 0 else set())
        falta = tuple(sorted(r for r in range(q) for _ in range(c - A[r]) if A[r] < c))
        if v == 0:
            fam = "v=0    : A = c.mu_q EXACTO           (alfabeto completo)"
        elif v == q - 1:
            fam = "v=q-1  : A = c.mu_q - {1, 1}         (dos veces el punto fijo r=0)"
        elif v == a and q % 2:
            fam = "v=a, q impar : A = c.mu_q - {1}      (UNA letra, el punto fijo r=0)"
        elif v == a:
            fam = "v=a, q par   : A = c.mu_q - {1, -1}  (los dos puntos fijos)"
        else:
            fam = "esporadica q=6, v=1"
        vistos.setdefault(fam, []).append((m, q, falta))
for fam in sorted(vistos):
    ej = vistos[fam][:3]
    print("  %s" % fam)
    print("       ej: " + " | ".join("m=%d,q=%d faltan r=%s" % (m, q, list(f) if f else "nada")
                                     for m, q, f in ej))
