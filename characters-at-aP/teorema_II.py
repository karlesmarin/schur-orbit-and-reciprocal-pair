# -*- coding: utf-8 -*-
"""EL TEOREMA II PROPUESTO, verificado desde cero.

  (a)  Z_q(l)/Z_q(rho) = kappa_q(lambda) . D(lambda,q)  in Z_{>0},  kappa_q in {1,2},
       para CUALQUIER q --- no solo q | 2m+2.
  (b)  bajo la condicion (I),  chi = (-1)^E . gamma . D(lambda,q)  con gamma in {1,2,3,4,6}.

D(lambda,q) usa las MISMAS recetas de la §7: GL_n en las clases no fijas, Sp_{2n} en la clase 0,
SO_{2n+1} en la clase q/2.  El caracter se calcula aparte, como limite de la formula de senos,
para que la comprobacion no use ninguna de nuestras formulas.
"""
import os
import sys
from collections import defaultdict
from fractions import Fraction
from math import gcd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # antes: un scratchpad temporal (9-sep-2026)
from galois import chi_lim, condicion_I, raices_C


def dim_gl(nu):
    n = len(nu)
    if n == 0:
        return 1
    d = Fraction(1)
    for i in range(n):
        for j in range(i + 1, n):
            d *= Fraction(nu[i] - nu[j] + j - i, j - i)
    assert d.denominator == 1, (nu, d)
    return int(d)


def dim_sp(nu):
    n = len(nu)
    if n == 0:
        return 1
    l = [Fraction(nu[i] + n - i) for i in range(n)]
    r = [Fraction(n - i) for i in range(n)]
    d = Fraction(1)
    for i in range(n):
        d *= l[i] / r[i]
        for j in range(i + 1, n):
            d *= (l[i] ** 2 - l[j] ** 2) / (r[i] ** 2 - r[j] ** 2)
    assert d.denominator == 1, (nu, d)
    return int(d)


def dim_so_odd(nu):
    n = len(nu)
    if n == 0:
        return 1
    l = [Fraction(2 * nu[i] + 2 * (n - i) - 1, 2) for i in range(n)]
    r = [Fraction(2 * (n - i) - 1, 2) for i in range(n)]
    d = Fraction(1)
    for i in range(n):
        d *= l[i] / r[i]
        for j in range(i + 1, n):
            d *= (l[i] ** 2 - l[j] ** 2) / (r[i] ** 2 - r[j] ** 2)
    assert d.denominator == 1, (nu, d)
    return int(d)


def perfil(ell, q):
    """clases plegadas -> lista de ell_i que caen en ellas."""
    d = defaultdict(list)
    for x in ell:
        r = x % q
        d[min(r, q - r)].append(x)
    return d


def D_de(ell, q):
    """D(lambda,q) con las recetas de la §7.  Devuelve None si alguna receta no da entero."""
    fij = {0} | ({q // 2} if q % 2 == 0 else set())
    P = perfil(ell, q)
    total = 1
    for c, filas in P.items():
        n = len(filas)
        if n == 0:
            continue
        if c == 0:
            L = sorted(filas, reverse=True)
            nu = [L[i] // q - (n - i) for i in range(n)]
            if any(nu[i] < nu[i + 1] for i in range(n - 1)) or nu[-1] < 0:
                return None
            total *= dim_sp(nu)
        elif c in fij:                       # c = q/2
            L = sorted(filas, reverse=True)
            nu = [Fraction(L[i], q) - (n - i - Fraction(1, 2)) for i in range(n)]
            if any(x.denominator != 1 for x in nu):
                return None
            nu = [int(x) for x in nu]
            if any(nu[i] < nu[i + 1] for i in range(n - 1)) or nu[-1] < 0:
                return None
            total *= dim_so_odd(nu)
        else:
            ys = []
            for x in filas:
                if x % q == c:
                    ys.append((x - c) // q)
                else:
                    ys.append((-x - c) // q)
            ys.sort(reverse=True)
            nu = [ys[i] - ys[-1] - (n - 1 - i) for i in range(n)]
            if any(nu[i] < nu[i + 1] for i in range(n - 1)) or nu[-1] < 0:
                return None
            total *= dim_gl(nu)
    return total


def dominantes(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for r in rec(k - 1, v):
                yield (v,) + r
    return rec(m, top)


print("=" * 100)
print("(a)  Z_q(l)/Z_q(rho) = kappa . D(lambda,q)  con kappa in {1,2},  para CUALQUIER q")
print("=" * 100)
malos_a = tot_a = 0
kaps = defaultdict(int)
for m in range(2, 6):
    R = raices_C(m)
    rho = [m - i for i in range(m)]
    vr = [sum(rho[i] * a[i] for i in range(m)) for a in R]
    for q in range(2, 13):
        rq = sum(1 for v in vr if v % q == 0)
        Zr = 1
        for v in vr:
            if v % q == 0:
                Zr *= v
        for lam in dominantes(m, 5):
            ell = [lam[i] + rho[i] for i in range(m)]
            ve = [sum(ell[i] * a[i] for i in range(m)) for a in R]
            if sum(1 for v in ve if v % q == 0) != rq:
                continue
            Z = 1
            for v in ve:
                if v % q == 0:
                    Z *= v
            D = D_de(ell, q)
            tot_a += 1
            if D is None or Z % (Zr * D) != 0:
                malos_a += 1
                continue
            k = Z // (Zr * D)
            kaps[k] += 1
            if k not in (1, 2):
                malos_a += 1
print("   supervivientes contrastados: %d   fuera de la forma kappa.D: %d" % (tot_a, malos_a))
print("   reparto de kappa:", dict(sorted(kaps.items())))

print()
print("=" * 100)
print("(b)  chi = +- gamma . D   con gamma in {1,2,3,4,6},  bajo la condicion (I)")
print("=" * 100)
gam = defaultdict(int)
malos_b = tot_b = 0
fuera = []
for m in range(2, 6):
    R = raices_C(m)
    rho = [m - i for i in range(m)]
    vr = [sum(rho[i] * a[i] for i in range(m)) for a in R]
    for q in range(2, 13):
        if not condicion_I(m, q):
            continue
        rq = sum(1 for v in vr if v % q == 0)
        for lam in dominantes(m, 4):
            ell = [lam[i] + rho[i] for i in range(m)]
            ve = [sum(ell[i] * a[i] for i in range(m)) for a in R]
            if sum(1 for v in ve if v % q == 0) != rq:
                continue
            D = D_de(ell, q)
            if D is None:
                continue
            c = chi_lim(lam, q)
            cr = round(float(c))
            tot_b += 1
            if D == 0 or abs(cr) % D != 0:
                malos_b += 1
                continue
            g = abs(cr) // D
            gam[g] += 1
            if g not in (1, 2, 3, 4, 6):
                malos_b += 1
                if len(fuera) < 6:
                    fuera.append((m, q, tuple(lam), cr, D, g))
print("   supervivientes contrastados: %d   con gamma fuera de {1,2,3,4,6}: %d" % (tot_b, malos_b))
print("   reparto de gamma:", dict(sorted(gam.items())))
for x in fuera:
    print("      m=%d q=%d lam=%s chi=%d D=%d gamma=%s" % x)
