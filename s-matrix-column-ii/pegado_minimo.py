# -*- coding: utf-8 -*-
r"""P3_pegado_minimo.py -- el pegado reducido a su minima expresion, y comprobado.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Forma minima (k = p - n):
    gl(n,p)/C(p-1,n-1) = [(n-2)(k-2) - 2] / (2(p-1))  +  p/(2 C(p,n)) * sum_{[T]} (1 - 2 delta(T)),
la suma sobre las orbitas multiplicativas de n-conjuntos centrados de F_p.  Cada termino es simetrico
en n <-> k (C(p,n) = C(p,k); el complemento empareja orbitas y conserva delta).
Se comprueba contra la enumeracion (P3_pegado_su3.pegado) para 3 <= n <= p-3, p <= 19, y se
reescriben SU(3) y SU(4) con sympy en su forma mas corta.
"""
import sys
from fractions import Fraction
from math import comb
from itertools import combinations
import sympy as sp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
_src = open(__file__.replace("P3_pegado_minimo.py", "P3_pegado_su3.py"), encoding="utf-8").read()
_src = _src.split("for p in primos(5, PMAX):")[0].replace("PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 43", "")
_src = _src.replace('__file__.replace("P3_pegado_su3.py", "P3_socratico.py")',
                    repr(__file__.replace("P3_pegado_minimo.py", "P3_socratico.py")))
_src = _src.replace('repr(__file__.replace("P3_pegado_su3.py", "P3_delta.py"))',
                    repr(repr(__file__.replace("P3_pegado_minimo.py", "P3_delta.py"))))
exec(_src)

malos = tot = 0
for p in primos(5, 19):
    for n in range(2, p - 1):
        k = p - n
        vistos = set(); suma = 0
        for T in combinations(range(p), n):
            if sum(T) % p:
                continue
            key = min(tuple(sorted(a * t % p for t in T)) for a in range(1, p))
            if key in vistos:
                continue
            vistos.add(key)
            h, e, gG, gS, cS = datos(list(T), p)
            suma += 1 - 2 * gS
        minima = Fraction((n - 2) * (k - 2) - 2, 2 * (p - 1)) + Fraction(p, 2 * comb(p, n)) * suma
        real = pegado(n, p) / comb(p - 1, n - 1)
        tot += 1
        if minima != real:
            malos += 1
            print("  FALLA n=%d p=%d: %s vs %s" % (n, p, minima, real))
print("forma minima: %d pares (n,p), fallos %d" % (tot, malos))

p = sp.symbols('p')
g3 = sp.factor(((p - 5) / 2) ** 2)
print("SU(3): gl = ((p-5)/2)^2 + 2*[p=1 mod 3]   (", g3, ")")
q4 = sp.Rational(1, 12) * (2 * p ** 3 - 25 * p ** 2 + 92 * p - 97)
x = sp.symbols('x')
print("SU(4) parte polinomica en x = p-5:", sp.expand(q4.subs(p, x + 5)))
print("SU(4) parte polinomica factorizada:", sp.factor(q4 - sp.Rational(-1)), " ... y en (p-3)(p-5):",
      sp.expand(q4 - (p - 3) * (p - 5) * (2 * p - 7) / 12))
