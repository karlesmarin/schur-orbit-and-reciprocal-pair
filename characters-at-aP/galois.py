# -*- coding: utf-8 -*-
"""LA CLASIFICACION (I): en que puntos g_{1/q} de la recta son ENTEROS todos los caracteres.

Afirma:   todos enteros  <=>  q|2m  o  q|(2m+1)  o  q|(2m+2)  o  (q=6 y m = 1 mod 3).

La via de prueba es: todos los caracteres enteros <=> el espectro de la natural es estable
bajo Galois, xi_q -> xi_q^p para todo p coprimo con q.  Eso es finito y barato: el espectro es
el multiconjunto A = {+-j mod q : j=1..m}.

Se contrastan DOS cosas:
  (1) estabilidad de Galois del espectro   contra   la condicion (I);
  (2) y, en rangos pequenos, la integridad REAL de los caracteres, calculados como limite
      de la formula de senos, contra la estabilidad de Galois.
"""
import sys
from collections import Counter
from math import gcd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from mpmath import mp, mpf, sin as msin, pi as mpi, nstr
mp.dps = 40


def espectro(m, q):
    c = Counter()
    for j in range(1, m + 1):
        c[j % q] += 1
        c[(-j) % q] += 1
    return c


def galois_estable(m, q):
    A = espectro(m, q)
    for p in range(2, q):
        if gcd(p, q) != 1:
            continue
        B = Counter()
        for r, n in A.items():
            B[(p * r) % q] += n
        if B != A:
            return False
    return True


def condicion_I(m, q):
    return (2 * m) % q == 0 or (2 * m + 1) % q == 0 or (2 * m + 2) % q == 0 \
        or (q == 6 and m % 3 == 1)


print("=" * 96)
print("(1) ESTABILIDAD DE GALOIS DEL ESPECTRO   contra   LA CONDICION (I)")
print("=" * 96)
ac = des = 0
casos = []
for q in range(2, 61):
    for m in range(1, 201):
        g = galois_estable(m, q)
        c = condicion_I(m, q)
        if g == c:
            ac += 1
        else:
            des += 1
            if len(casos) < 12:
                casos.append((m, q, g, c))
print("  pares (m,q) contrastados: %d    ACUERDOS %d    DESACUERDOS %d" % (ac + des, ac, des))
for m, q, g, c in casos:
    print("     m=%-4d q=%-3d  galois=%-5s  condicion=%-5s" % (m, q, g, c))

print()
print("=" * 96)
print("(2) INTEGRIDAD REAL DE LOS CARACTERES  contra  ESTABILIDAD DE GALOIS")
print("    chi como limite de  prod_{alpha>0} sin(pi p (l,a)/q) / sin(pi p (rho,a)/q)")
print("=" * 96)


def raices_C(m):
    R = []
    for i in range(m):
        for j in range(i + 1, m):
            v = [0] * m; v[i] = 1; v[j] = -1; R.append(tuple(v))
            v = [0] * m; v[i] = 1; v[j] = 1; R.append(tuple(v))
    for i in range(m):
        v = [0] * m; v[i] = 2; R.append(tuple(v))
    return R


def chi_lim(lam, q, p=1):
    m = len(lam)
    R = raices_C(m)
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    eps = mpf(10) ** (-25)
    th = mpf(p) / q + eps
    num = mpf(1)
    for a in R:
        vl = sum(ell[i] * a[i] for i in range(m))
        vr = sum(rho[i] * a[i] for i in range(m))
        num *= msin(mpi * th * vl) / msin(mpi * th * vr)
    return num


def dominantes(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for r in rec(k - 1, v):
                yield (v,) + r
    return rec(m, top)


tol = mpf(10) ** (-12)
ac2 = des2 = 0
for m in range(2, 6):
    for q in range(2, 13):
        todos = True
        peor = None
        for lam in dominantes(m, 5):
            v = chi_lim(lam, q)
            r = round(float(v))
            if abs(v - r) > tol:
                todos = False
                if peor is None:
                    peor = (tuple(lam), nstr(v, 10))
                break
        g = galois_estable(m, q)
        if todos == g:
            ac2 += 1
        else:
            des2 += 1
            print("     DESACUERDO m=%d q=%-2d  enteros=%-5s galois=%-5s  %s"
                  % (m, q, todos, g, peor or ""))
print("  pares (m,q): %d   ACUERDOS %d   DESACUERDOS %d" % (ac2 + des2, ac2, des2))
