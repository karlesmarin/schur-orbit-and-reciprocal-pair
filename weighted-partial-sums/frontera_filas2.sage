# -*- coding: utf-8 -*-
# FRONTERA 3, el detalle: la LEY MOD 4.
#
# El barrido anterior dio "filas con solucion: 3,4,7,8,11,12,15,16,19,20,23,24,27,28" -- todas
# m = 0 o 3 modulo 4.  La razon es una obstruccion de PARIDAD:
#
#     sum_{c<=m} c eps(c)  =  sum_{c<=m} c  =  m(m+1)/2   (mod 2),
#
# porque cambiar un signo cambia la suma en 2c.  Para que la suma se anule el NUMERO TRIANGULAR
# T_m = m(m+1)/2 tiene que ser PAR, y eso ocurre exactamente cuando m = 0 o 3 modulo 4.
#
# Dicho de otro modo: una fila puede tener ceros solo si {1,...,m} admite una particion en dos
# partes de igual suma; y las que ocurren son exactamente las MULTIPLICATIVAS.
#
# Aqui se comprueban las dos direcciones hasta m = 44:
#   (a) m = 1,2 mod 4  =>  CERO soluciones  (la obstruccion);
#   (b) m = 0,3 mod 4  =>  ?siempre hay al menos una?  (el reciproco, que no es obvio: la
#       particion generica existe, pero se exige que sea multiplicativa).
#
# Authors: Carles Marin, Claude (AI assistant).
import sys
from itertools import product


def cuenta(m):
    ps = prime_range(m + 1)
    n = 0
    ejemplos = []
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
            if len(ejemplos) < 1:
                ejemplos.append(asig)
    return n, ejemplos


print("%-4s %-6s %-8s %-8s %-9s %s" % ("m", "m mod 4", "T_m par", "#primos", "#soluc", "veredicto"))
mal_a = mal_b = 0
sec = []
for m in range(2, 45):
    Tm = m * (m + 1) // 2
    par = (Tm % 2 == 0)
    n, ej = cuenta(m)
    sec.append(n)
    espera_cero = not par
    if espera_cero and n > 0:
        ver = "*** rompe la obstruccion ***"
        mal_a += 1
    elif (not espera_cero) and n == 0:
        ver = "*** reciproco falla ***"
        mal_b += 1
    else:
        ver = "ok"
    print("%-4d %-6d %-8s %-8d %-9d %s" % (
        m, m % 4, "si" if par else "no", len(prime_range(m + 1)), n, ver))
    sys.stdout.flush()
print()
print("obstruccion (T_m impar => 0 soluciones): %d fallos" % mal_a)
print("reciproco (T_m par => alguna solucion): %d fallos" % mal_b)
print("sucesion de #soluciones para m = 2..44:")
print("  " + ", ".join(str(x) for x in sec))
