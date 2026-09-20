# -*- coding: utf-8 -*-
r"""altura_compuesta.py -- columns at composite height q = n+k: Gorenstein test and a non-monomial order.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

(1) Gorenstein test by integer lattices (P3a.orden_columna): A = Z[e_1..e_{n-1}] of the eigenvalues
    of g_mu, O = its integral closure (the H-fixed part of Z[zeta_N], N = nq), conductor
    f = {x in O : xO in A}; always [A:f] <= [O:A], with equality iff A is Gorenstein.
    Controls at prime height 17: SU(5)_12, mu = (0,0,2,1) (semigroup <3,4,5>: [O:A] = 17^2,
    [A:f] = 17, not Gorenstein) and mu = 0 (maximal).
    Composite height 14: SU(3)_11, mu = (0,1), and every column of SU(3)_11.
(2) SU(3)_6, q = 9, mu = (1,0): in Z_3[zeta_27] modulo (3, lambda^6), lambda = zeta - 1, the two
    fundamental characters, and that no cube of a uniformiser lies in the image of the order.
Public domain (CC0).
"""
import sys
from math import comb, gcd

from P3a import orden_columna, exponentes, factor_str

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MAL = []


def check(ok, msg):
    print(("  ok     " if ok else "  FAILED ") + msg)
    if not ok:
        MAL.append(msg)
    sys.stdout.flush()


def show(n, k, mu):
    r = orden_columna(n, k, mu)
    print("   SU(%d)_%d q=%d mu=%s: exponents mod %d = %s, degree %d, [O:A] = %s, [A:f] = %s, [O:f] = %s -> %s"
          % (n, k, n + k, tuple(mu), r["N"], r["exponentes"], r["grado"], factor_str(r["OA"]),
             factor_str(r["Af"]), factor_str(r["Of"]), "Gorenstein" if r["gorenstein"] else "NOT Gorenstein"))
    return r


print("(1) Gorenstein test")
a = show(5, 12, (0, 0, 2, 1))
b = show(5, 12, (0, 0, 0, 0))
check((a["OA"], a["Af"], a["gorenstein"]) == (17 ** 2, 17, False) and (b["OA"], b["Af"]) == (1, 1),
      "controls at prime height as predicted by the semigroups <3,4,5> and N_0")
t = show(3, 11, (0, 1))
check((t["OA"], t["Af"], t["gorenstein"]) == (2 ** 8, 2 ** 4, False),
      "SU(3)_11 (q = 14), mu = (0,1): [O:A] = 2^8, [A:f] = 2^4, not Gorenstein")
malos = []
for x in range(12):
    for y in range(12 - x):
        r = orden_columna(3, 11, (x, y))
        if r["Af"] > r["OA"]:
            MAL.append("[A:f] > [O:A] at %s" % ((x, y),))
        if not r["gorenstein"]:
            malos.append(((x, y), factor_str(r["OA"]), factor_str(r["Af"])))
print("   every column of SU(3)_11 (78): %d not Gorenstein:" % len(malos))
for m in malos:
    print("      mu=%s [O:A]=%s [A:f]=%s" % m)

print("(2) SU(3)_6, q = 9, mu = (1,0): the order is not monomial")
X, N = exponentes(3, 6, (1, 0))
check(sorted(X) == [5, 23, 26] and N == 27, "eigenvalues zeta_27^%s" % sorted(X))
H = [c for c in range(1, 27) if gcd(c, 3) == 1 and sorted(c * x % 27 for x in X) == sorted(X)]
check(H == [1], "trivial stabiliser in (Z/27)^x: the field is Q(zeta_27), totally ramified at 3")
R = 6


def pot(e):
    """zeta^e = (1 + lambda)^e modulo (3, lambda^6)."""
    return [comb(e, i) % 3 for i in range(R)]


def suma(*vs):
    return [sum(x) % 3 for x in zip(*vs)]


e1 = suma(*[pot(x) for x in X])
e2 = suma(*[pot((x + y) % 27) for x, y in ((X[0], X[1]), (X[0], X[2]), (X[1], X[2]))])
check(e1 == [0, 0, 0, 1, 2, 1] and e2 == [0, 0, 0, 2, 2, 0],
      "e_1 = lambda^3 + 2 lambda^4 + lambda^5, e_2 = 2 lambda^3 + 2 lambda^4  mod (3, lambda^6)")
check(3 + 3 >= R, "products of e_1, e_2 have valuation >= 6: the image of the order is span{1, e_1, e_2}")
imagen = {tuple(suma([u * c for c in e1], [v * c for c in e2], [w] + [0] * (R - 1)))
          for u in range(3) for v in range(3) for w in range(3)}
check(all((0, 0, 0, c, 0, 0) not in imagen for c in (1, 2)),
      "no u^3 lambda^3 (u a unit) lies in the image: no cube of a uniformiser is in the order")
check(any(x[3] and not any(x[:3]) for x in imagen),
      "yet the image has elements of valuation exactly 3 (e_1 itself)")
print()
print("FAILED: %d %s" % (len(MAL), MAL))
