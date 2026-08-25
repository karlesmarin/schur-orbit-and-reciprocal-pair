# -*- coding: utf-8 -*-
# EL BIALTERNANTE CONTRA .expand():  MISMO VALOR, OTRO COSTE.   19 de agosto de 2026.
#
# POR QUE.  `schur_at` de pIII_universal_t3.sage evalua el caracter asi:
#
#     sch(Partition(lam)).expand(NLET, alphabet=...)  y luego sustituye
#
# o sea que construye el POLINOMIO SIMETRICO entero en 5 variables y despues especializa.  El numero
# de monomios explota con lambda_1, y medido en _probe_t3_coste.sage el coste crece x10 cada 3 de
# |lambda|: 6 s en (8,4), 59 s en (10,5), y (12,6) no termino en 12 minutos.  El bucle T0 entero
# salia a ~39 horas.
#
# LA RUTA DEL GRUPO.  s_lambda es un caracter de GL_5 y la formula de caracteres de Weyl lo da como
# BIALTERNANTE --- un cociente de dos determinantes 5x5:
#
#     s_lambda(x_1..x_n) = det( x_i^{lambda_j + n - j} ) / det( x_i^{n - j} )
#
# con el denominador igual al Vandermonde prod_{i<j}(x_i - x_j).  El coste es un determinante 5x5 de
# monomios: NO depende de lambda mas que por elevar a una potencia.  Los cinco puntos del alfabeto
# --- 1, w, w^2, z, z^{-1} --- son DISTINTOS en el anillo, que es la unica hipotesis que pide.
#
# CONTROLES
#   B1  FATAL: el bialternante = .expand() en todas las lambdas donde .expand() es asequible.  Si no
#       coincide, no vale nada de lo de abajo.  Es la calibracion del instrumento, y va primero.
#   B2  el coste de las dos rutas, cronometrado, en lambdas crecientes.
#   B3  el bialternante en lambdas que .expand() NO alcanza --- (12,6), (16,8), (24,12) --- para ver
#       si de verdad es plano en lambda.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local bash -c "sage < _probe_bialternante.sage > _probe_bialternante_OUT.txt 2>&1"

import time
import sys

T    = 3
NLET = T + 2

K   = CyclotomicField(6 * T)
g   = K.gen()
w   = g ** 6
zK  = g ** 2
LzK = LaurentPolynomialRing(K, 'zz')
zz  = LzK.gen()

sch = SymmetricFunctions(QQ).schur()
RN  = PolynomialRing(K, NLET, 'x')

orbit = [LzK(w ** j) for j in range(T)]
alphA = orbit + [zz, zz ** -1]


def schur_expand(lam, alphabet):
    """La ruta vieja: construir el polinomio simetrico y sustituir."""
    parts = [l for l in lam if l > 0]
    if not parts:
        return alphabet[0].parent().one()
    e = sch(Partition(sorted(parts, reverse=True))).expand(NLET, alphabet=list(RN.gens()))
    return e(*alphabet)


def schur_bialt(lam, alphabet):
    """La ruta del grupo: formula de caracteres de Weyl para GL_n, un cociente de dos 5x5."""
    n = len(alphabet)
    lam = list(lam) + [0] * (n - len(lam))
    lam = sorted(lam, reverse=True)
    num = matrix(LzK, n, n, lambda i, j: alphabet[i] ** (lam[j] + n - 1 - j)).det()
    den = LzK.one()
    for i in range(n):
        for j in range(i + 1, n):
            den *= (alphabet[i] - alphabet[j])
    q, r = num.quo_rem(den)
    if r != 0:
        raise ArithmeticError("el bialternante no divide exacto en lambda=%s" % (lam,))
    return q


print("=" * 96)
print("EL BIALTERNANTE CONTRA .expand()")
print("=" * 96)
print("")

print("  B1  FATAL: calibracion --- el bialternante = .expand() donde .expand() es asequible")
ok = bad = 0
malas = []
for a in range(0, 9):
    for b in range(0, a + 1):
        lam = (a, b)
        try:
            v1 = LzK(schur_expand(lam, alphA))
            v2 = LzK(schur_bialt(lam, alphA))
            if v1 == v2:
                ok += 1
            else:
                bad += 1
                malas.append(lam)
        except Exception as e:
            bad += 1
            malas.append((lam, str(e)[:60]))
print("      %d de %d   %s" % (ok, ok + bad, "PASA" if bad == 0 else "FALLA en %s" % (malas[:5],)))
sys.stdout.flush()
if bad:
    print("      no se sigue adelante con un instrumento sin calibrar")
    print("DONE")
else:
    print("")
    print("  B2  coste de las dos rutas")
    print("      %-10s %14s %14s   %s" % ("lambda", ".expand()", "bialternante", "ganancia"))
    for lam in [(4, 2), (6, 3), (8, 4), (10, 5)]:
        t0 = time.time(); schur_expand(lam, alphA); d1 = time.time() - t0
        t0 = time.time(); schur_bialt(lam, alphA);  d2 = time.time() - t0
        print("      %-10s %12.3f s %12.3f s   %s"
              % (str(lam), d1, d2, ("x%.0f" % (d1 / d2)) if d2 > 0 else "-"))
        sys.stdout.flush()
    print("")
    print("  B3  el bialternante donde .expand() no llega")
    for lam in [(12, 6), (12, 12), (16, 8), (24, 12), (60, 30)]:
        t0 = time.time(); schur_bialt(lam, alphA); d = time.time() - t0
        print("      lambda=%-10s  %8.3f s" % (str(lam), d))
        sys.stdout.flush()
    print("")
    print("  Lectura: si B1 pasa y B3 es plano, el bucle T0 con LMAX=12 deja de costar 39 horas.")
    print("")
    print("=" * 96)
    print("DONE")
