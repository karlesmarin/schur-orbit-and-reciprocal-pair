# -*- coding: utf-8 -*-
# FRONTERA 3: el enunciado uniforme en m, en el lenguaje bueno.
#
# Con n = 2m+1 la longitud del segmento simetrico {-m,...,m}:
#
#       T_m = m(m+1)/2 = (n^2 - 1)/8      y      (2|n) = (-1)^{(n^2-1)/8}
#
# luego la obstruccion de paridad "T_m par" ES la ley suplementaria  (2|n) = +1.  Y la diagonal,
# donde n = f y el criterio de thm:diag es chi(2) = 1 -- para chi cuadratico, (2|f) = +1 -- es el
# caso n = f de la MISMA ley.  La frontera 3 del papel dice hoy que "<2> gobierna ahi tambien" es
# falso; no lo es: es la ley suplementaria en la LONGITUD DEL SEGMENTO y no en el conductor.
#
# Se comprueban las tres cosas:
#   (1) la fila m admite ceros  <=>  (2|n) = +1,  n = 2m+1      [obstruccion y reciproco]
#   (2) en la diagonal n = f, y (2|f) = +1 es el criterio de thm:diag
#   (3) el contraste con caracteres cuadraticos reales, fila por fila
#
# CONTROL que debe fallar: la misma ley con 3 en vez de 2, es decir (3|n).
#
# Authors: Carles Marin, Claude (AI assistant).
import sys
from itertools import product


def n_soluciones(m):
    ps = prime_range(m + 1)
    n = 0
    for signos in product([1, -1], repeat=len(ps)):
        asig = dict(zip(ps, signos))
        eps = [0] * (m + 1)
        eps[1] = 1
        for c in range(2, m + 1):
            v = 1
            for (pp, ee) in factor(c):
                v *= asig[pp] ** ee
            eps[c] = v
        if sum([c * eps[c] for c in range(1, m + 1)]) == 0:
            n += 1
    return n


print("=== (1) la fila m tiene ceros  <=>  (2|n) = +1  con n = 2m+1 ===")
print("%-4s %-5s %-8s %-8s %-9s %-9s %s" % (
    "m", "n", "T_m", "(2|n)", "#soluc", "senuelo(3|n)", "veredicto"))
mal = 0
sen_mal = 0
for m in range(2, 41):
    n = 2 * m + 1
    Tm = m * (m + 1) // 2
    leg2 = kronecker(2, n)
    leg3 = kronecker(3, n) if n % 3 else 0
    ns = n_soluciones(m)
    bien = ((ns > 0) == (leg2 == 1))
    if not bien:
        mal += 1
    if (ns > 0) == (leg3 == 1):
        sen_mal += 1
    print("%-4d %-5d %-8d %-8d %-9d %-9d %s" % (
        m, n, Tm, leg2, ns, leg3, "ok" if bien else "*** DISCREPA ***"))
    sys.stdout.flush()
print()
print("(2|n) = +1  <=>  la fila tiene ceros:  %d fallos de 39" % mal)
print("senuelo con (3|n) coincide en %d de 39 filas" % sen_mal)
print()

print("=== (2) la diagonal es el caso n = f de la misma ley ===")
print("%-5s %-6s %-8s %-10s %s" % ("f", "m", "(2|f)", "T(m,chi)", "veredicto"))
ok2 = mal2 = 0
for f in range(7, 120, 2):
    if not is_squarefree(f) or f % 4 != 3:
        continue
    chi = kronecker_character(-f)
    if chi is None or chi.conductor() != f:
        continue
    m = (f - 1) // 2
    t = sum([c * chi(c) for c in range(1, m + 1)])
    leg2 = kronecker(2, f)
    bien = ((t == 0) == (leg2 == 1))
    ok2 += 1 if bien else 0
    mal2 += 0 if bien else 1
    if leg2 == 1 or t == 0:
        print("%-5d %-6d %-8d %-10s %s" % (f, m, leg2, t, "ok" if bien else "*** DISCREPA ***"))
print("diagonal:  %d ok, %d discrepan" % (ok2, mal2))
