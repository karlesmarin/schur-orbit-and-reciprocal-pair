# -*- coding: utf-8 -*-
# COSTE del bialternante contra .expand().   19 de agosto de 2026.
#
# Continuacion de _probe_bialternante.sage, que ya calibro (B1: 45 de 45) y se corto ahi.
#
# POR QUE SE CORTO, y es una trampa que vale por si sola: `sage < fichero` alimenta el fichero al
# REPL, y en el REPL un `else:` a nivel superior --- despues de que el bloque `if` ya se ejecuto al
# ver una linea sin indentar --- no sobrevive.  El guion murio con EXIT 0 y sin decir nada.  Regla:
# en un `.sage` que se vaya a ejecutar por stdin, NADA de if/else, try/except ni elif a nivel
# superior.  Todo plano, o dentro de una funcion.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run -d --name pIII-bialt2 -v "E:/.../gates:/work" -w /work \
#         sage-normaliz:local bash -c "sage < _probe_bialt2.sage > _probe_bialt2_OUT.txt 2>&1"

import time
import sys

T    = 3
NLET = T + 2
K    = CyclotomicField(6 * T)
g    = K.gen()
w    = g ** 6
LzK  = LaurentPolynomialRing(K, 'zz')
zz   = LzK.gen()
sch  = SymmetricFunctions(QQ).schur()
RN   = PolynomialRing(K, NLET, 'x')
alphA = [LzK(w ** j) for j in range(T)] + [zz, zz ** -1]


def schur_expand(lam, alphabet):
    parts = [l for l in lam if l > 0]
    if not parts:
        return alphabet[0].parent().one()
    e = sch(Partition(sorted(parts, reverse=True))).expand(NLET, alphabet=list(RN.gens()))
    return e(*alphabet)


def schur_bialt(lam, alphabet):
    n = len(alphabet)
    L = sorted(list(lam) + [0] * (n - len(lam)), reverse=True)
    num = matrix(LzK, n, n, lambda i, j: alphabet[i] ** (L[j] + n - 1 - j)).det()
    den = LzK.one()
    for i in range(n):
        for j in range(i + 1, n):
            den *= (alphabet[i] - alphabet[j])
    q, r = num.quo_rem(den)
    return q


def crono(f, lam):
    t0 = time.time()
    f(lam, alphA)
    return time.time() - t0


print("=" * 96)
print("COSTE: .expand() contra el bialternante")
print("=" * 96)
print("")
print("  B2  las dos rutas, cronometradas")
print("      %-10s %14s %14s   %s" % ("lambda", ".expand()", "bialternante", "ganancia"))
sys.stdout.flush()

for lam in [(4, 2), (6, 3), (8, 4), (10, 5)]:
    d1 = crono(schur_expand, lam)
    d2 = crono(schur_bialt, lam)
    print("      %-10s %12.3f s %12.3f s   x%.0f" % (str(lam), d1, d2, d1 / max(d2, 1e-9)))
    sys.stdout.flush()

print("")
print("  B3  el bialternante donde .expand() no llega (con LMAX=12 el bucle T0 salia a ~39 h)")
sys.stdout.flush()

for lam in [(12, 6), (12, 12), (16, 8), (24, 12), (60, 30), (200, 100)]:
    d = crono(schur_bialt, lam)
    print("      lambda=%-12s  %8.3f s" % (str(lam), d))
    sys.stdout.flush()

print("")
print("  B4  el bucle T0 ENTERO por bialternante, con LMAX=12: las 91 lambdas")
sys.stdout.flush()
t0 = time.time()
n = 0
for a in range(13):
    for b in range(a + 1):
        schur_bialt((a, b), alphA)
        n += 1
print("      %d lambdas en %.2f s" % (n, time.time() - t0))
print("")
print("=" * 96)
print("DONE")
