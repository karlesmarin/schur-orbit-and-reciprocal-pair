# -*- coding: utf-8 -*-
# SONDA DE COSTE del bucle T0 de pIII_universal_t3.sage.   19 de agosto de 2026.
#
# POR QUE.  El bucle T0 recorre 91 lambdas y NO imprime nada por iteracion --- `say()` solo se llama
# despues del bucle ---, asi que "44 minutos en silencio" es indistinguible de "colgado".  Esto no
# razona sobre el coste: lo mide, en un contenedor APARTE, sin tocar la corrida viva.
#
# Cronometra schur_at para lambdas de tamano creciente y extrapola las 91.  Si el total estimado es
# de horas, la corrida esta bien y hay que esperar (o partirla); si es de minutos, esta colgada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local bash -c "sage < _probe_t3_coste.sage"

import time
import sys

T    = 3
NLET = T + 2
LMAX = 12

K   = CyclotomicField(6 * T)
g   = K.gen()
w   = g ** 6
zK  = g ** 2
Lz  = LaurentPolynomialRing(QQ, 'z')
LzK = LaurentPolynomialRing(K, 'zz')

sch = SymmetricFunctions(QQ).schur()
RN  = PolynomialRing(K, NLET, 'x')

orbit = [LzK(w ** j) for j in range(T)]
alphA = orbit + [LzK.gen(), LzK.gen() ** -1]


def schur_at(lam, alphabet):
    parts = [l for l in lam if l > 0]
    if not parts:
        return alphabet[0].parent().one()
    e = sch(Partition(sorted(parts, reverse=True))).expand(NLET, alphabet=list(RN.gens()))
    return e(*alphabet)


print("=" * 90)
print("SONDA DE COSTE: cuanto cuesta UNA lambda del bucle T0")
print("=" * 90)
print("")

lams = [(a, b) for a in range(LMAX + 1) for b in range(a + 1)]
print("  lambdas totales en T0: %d" % len(lams))
print("")

muestra = [(0, 0), (2, 1), (4, 2), (6, 3), (8, 4), (10, 5), (12, 6), (12, 12)]
tiempos = {}
for lam in muestra:
    t0 = time.time()
    try:
        schur_at(lam, alphA)
        dt = time.time() - t0
        tiempos[lam] = dt
        print("    lambda=%-8s  %8.2f s" % (str(lam), dt))
    except Exception as e:
        print("    lambda=%-8s  ERROR %s" % (str(lam), e))
    sys.stdout.flush()

print("")
if tiempos:
    # extrapolacion grosera: el coste crece con |lambda|, asi que se usa el tiempo medido mas
    # cercano por tamano para cada una de las 91.
    conocidos = sorted(tiempos.items(), key=lambda kv: sum(kv[0]))
    total = 0.0
    for lam in lams:
        n = sum(lam)
        mejor = min(conocidos, key=lambda kv: abs(sum(kv[0]) - n))
        total += mejor[1]
    print("  ESTIMACION del bucle T0 entero: %.0f s  =  %.1f min  =  %.2f h" % (total, total / 60, total / 3600))
    print("")
    print("  Lectura: si esto son horas, la corrida de las 05:50 esta SANA y solo es lenta.")
    print("           si son minutos, esta colgada y hay que matarla.")
print("")
print("=" * 90)
print("DONE")
