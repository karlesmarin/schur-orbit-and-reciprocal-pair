# -*- coding: utf-8 -*-
r"""maximo.py -- the largest defect (n-1)(n-2)/2 and the split trinomials X^n + aX + b.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Proposition 'the largest defect': delta(T) <= (n-1)(n-2)/2, with equality exactly when
prod_{t in T}(X - t) = X^n + aX + b with ab != 0; then Gamma_T = S_T = <n-1, n>.
(1) For n = 4, 5 and every prime p <= 71: every trinomial X^n + aX + b, ab != 0, with n distinct
    roots in F_p gives a centred T with delta = (n-1)(n-2)/2 (semigroup of moments, and the exact
    semigroup); number of such trinomials against p^2/n!.
(2) Conversely, at the primes of (1), every centred T with h = 1 whose Gamma_T has genus
    (n-1)(n-2)/2 comes from such a trinomial (n = 4 all p <= 71; n = 5, p <= 41, by enumeration).
(3) The example of the note: p = 67, T = {9,27,50,52,63}, X^5 + 2X + 27, delta = 6 by lattices.
Public domain (CC0).
"""
import sys
from itertools import combinations
from math import factorial

from P3a import es_primo, gamma_T, s_exacto, semigrupo, delta_red, estabilizador

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def poly(T, p):
    pol = [1]
    for t in T:
        pol = [(b - t * a) % p for a, b in zip(pol + [0], [0] + pol)]   # low to high
    return pol


def es_trinomio(T, p):
    c = poly(T, p)      # c[0] + c[1] X + ... + X^n
    n = len(T)
    return all(c[i] == 0 for i in range(2, n)) and c[0] and c[1]


mal = 0
print("(1) split trinomials give the largest defect")
for n in (4, 5):
    top = (n - 1) * (n - 2) // 2
    for p in range(n + 2, 72):
        if not es_primo(p):
            continue
        cuenta = malos = 0
        for a in range(1, p):
            for b in range(1, p):
                raices = [x for x in range(p) if (pow(x, n, p) + a * x + b) % p == 0]
                if len(raices) == n:
                    cuenta += 1
                    G = semigrupo(gamma_T(raices, p))
                    S = semigrupo(s_exacto(raices, p))
                    if not (G["genero"] == S["genero"] == top and G["minimales"] == [n - 1, n]):
                        malos += 1
        mal += malos
        if cuenta:
            print("   n=%d p=%2d: %4d split trinomials (p^2/n! = %6.1f), delta = %d in all: %s"
                  % (n, p, cuenta, p * p / factorial(n), top, malos == 0))
    sys.stdout.flush()

print("(2) the converse, by enumeration of centred spectra with h = 1")
for n, pmax in ((4, 71), (5, 41)):
    top = (n - 1) * (n - 2) // 2
    for p in range(n + 2, pmax + 1):
        if not es_primo(p):
            continue
        k = m = 0
        for T in combinations(range(p), n):
            if sum(T) % p or len(estabilizador(T, p)) != 1:
                continue
            if semigrupo(gamma_T(T, p))["genero"] == top:
                k += 1
                if not es_trinomio(T, p):
                    m += 1
        mal += m
        if k or n == 5:
            print("   n=%d p=%2d: %4d spectra of genus %d, not from a trinomial: %d" % (n, p, k, top, m))
    sys.stdout.flush()

print("(3) the example")
T = [9, 27, 50, 52, 63]
c = poly(T, 67)
print("   p=67 T=%s: coefficients (X^0..X^5) %s ; Gamma_T = <%s> ; delta by lattices = %d"
      % (T, c, ",".join(map(str, semigrupo(gamma_T(T, 67))["minimales"])), delta_red(T, 67)))
print("failures: %d" % mal)
