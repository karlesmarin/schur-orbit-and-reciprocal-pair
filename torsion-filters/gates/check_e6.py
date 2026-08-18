# -*- coding: utf-8 -*-
# SU CONTRAEJEMPLO A "e <= 4", CON NUESTRA PROPIA ANATOMIA.   16 de agosto de 2026.
#
# El dice: t=6, r=3, beta = (12,11,10,9,8,7,5,4,3,2,1,0) tiene lambda = (1^6), e = 6, |g_com| = 4 y
# Phi = 0.  Se comprueba con anatomia() y phi_zero(), que son NUESTRAS y no saben de esta discusion.
#
# Y se comprueba tambien el ALCANCE de nuestro hallazgo de ayer: si |g_com| = e - 2 y e <= 2r,
# entonces con r = 2 forzosamente |g_com| <= 2, o sea nuestro "g_com es un par" era un artefacto del
# rango.  Se mide |g_com| sobre los ceros con r = 3 para verlo directamente.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python check_e6.py
import itertools
from collections import Counter

from peel_gcom import anatomia, betas
from peel_zero import phi_zero

print("=" * 100)
print("SU CONTRAEJEMPLO:  t=6, r=3, beta = (12,11,10,9,8,7,5,4,3,2,1,0)")
print("=" * 100)
beta = (12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 0)
t, r = 6, 3
N = t + 2 * r
delta = list(range(N - 1, -1, -1))
lam = tuple(beta[i] - delta[i] for i in range(N))
print("  N = %d,  lambda = %s" % (N, str(tuple(v for v in lam if v > 0))))
print("  clases mod 6: %s" % dict(sorted(Counter(v % t for v in beta).items())))
a = anatomia(beta, t, r)
if a is None:
    print("  anatomia(): None -- falta una clase o no hay exceso")
else:
    S, C, gcom, nG = a
    cl = Counter(v % t for v in beta)
    e = len([i for i in cl if cl[i] >= 2])
    print("  e (clases con >= 2) = %d" % e)
    print("  |G| = %d   (g_com solo esta definida si |G| = 2)" % nG)
    print("  S  = %s" % str(sorted(S)))
    print("  C  = %d" % C)
    print("  g_com = %s   |g_com| = %d   (la relacion |g_com| = e-2 daria %d)"
          % (str(sorted(gcom)), len(gcom), e - 2))
    print("  C - g_com == g_com ?  %s" % (set(C - v for v in gcom) == set(gcom)))
print("  Phi == 0 ?  %s" % (phi_zero(beta, t, r) is True))

print("")
print("=" * 100)
print("EL ALCANCE DE NUESTRO HALLAZGO: |g_com| sobre los ceros, r = 2 contra r = 3")
print("=" * 100)
for (t, r, W) in [(6, 2, 13), (6, 3, 13), (4, 3, 13)]:
    tam = Counter()
    ceros = 0
    for b in betas(t, r, W):
        an = anatomia(b, t, r)
        if an is None:
            continue
        S, C, gcom, nG = an
        if phi_zero(b, t, r) is not True:
            continue
        ceros += 1
        if nG == 2:
            tam[len(gcom)] += 1
    print("  t=%d r=%d W=%d :  %4d ceros | reparto de |g_com| entre los de |G|=2 : %s"
          % (t, r, W, ceros, dict(sorted(tam.items()))))
print("")
print("  LECTURA: si con r=3 aparece |g_com| > 2, nuestro 'g_com es un par' era del rango y no del")
print("  fenomeno, exactamente como el dice.")
print("=" * 100)
print("DONE")
