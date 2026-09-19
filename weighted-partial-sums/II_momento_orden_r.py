# -*- coding: utf-8 -*-
r"""II_momento_orden_r.py -- el segundo momento para valores en mu_r y la constante C_4, comprobados.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

f(p) independiente y uniforme en mu_r, completamente multiplicativa, S_r(m) = sum_{n<=m} n f(n).
(11)  E|S_r(m)|^2 = sum_{(a,b)=1} a^r b^r sum_{d <= m/max(a,b)^r} d^2
      (n = d a^r, k = d b^r: n/k es potencia r-esima exactamente cuando f(n) conj f(k) tiene media 1).
Comprobacion 1: media EXACTA sobre las r^pi(m) asignaciones (fuerza bruta) contra (11), r = 2, 3, 4.
Comprobacion 2: C_4 = (1/3) sum_{(a,b)=1} a^4 b^4 / max^12 contra la forma cerrada
      (2/(3 zeta(4))) (zeta(3)/5 + zeta(5)/3 - zeta(7)/30) = 0.340281785...
Comprobacion 3: E|S_4(m)|^2 / m^3 por (11) a m grande, tiende a C_4.
"""
import sys, math, cmath
from itertools import product
from fractions import Fraction
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def primos(m):
    return [p for p in range(2, m + 1) if all(p % q for q in range(2, int(p ** .5) + 1))]


def formula11(m, r):
    tot = 0
    b = 1
    while b ** r <= m:
        for a in range(1, b + 1):
            if math.gcd(a, b) != 1:
                continue
            M = m // (b ** r)
            s = M * (M + 1) * (2 * M + 1) // 6
            tot += (a ** r) * (b ** r) * s * (1 if a == b else 2)
        b += 1
    return tot


def bruta(m, r):
    P = primos(m)
    F = []
    for n in range(1, m + 1):
        v = []; x = n
        for i, p in enumerate(P):
            while x % p == 0:
                v.append(i); x //= p
        F.append(v)
    w = [cmath.exp(2j * math.pi * k / r) for k in range(r)]
    tot = 0.0
    for k in product(range(r), repeat=len(P)):
        s = sum((n + 1) * w[sum(k[i] for i in F[n]) % r] for n in range(m))
        tot += abs(s) ** 2
    return tot / r ** len(P)


for r, ms in ((2, (10, 16, 20)), (3, (10, 14)), (4, (10, 13))):
    for m in ms:
        fb, f11 = bruta(m, r), formula11(m, r)
        print("r=%d m=%2d  fuerza bruta=%.6f  (11)=%d  %s" % (r, m, fb, f11, "OK" if abs(fb - f11) < 1e-6 * f11 else "DISTINTO"))

z = lambda s: sum(1.0 / n ** s for n in range(1, 200000)) + 1 / ((s - 1) * 200000 ** (s - 1))
cerrada = 2 / (3 * z(4)) * (z(3) / 5 + z(5) / 3 - z(7) / 30)
serie = 0.0
for b in range(1, 4001):
    for a in range(1, b + 1):
        if math.gcd(a, b) == 1:
            serie += (a ** 4) * (b ** 4) / b ** 12 * (1 if a == b else 2)
serie /= 3
print("C_4 forma cerrada=%.9f  serie (b<=4000)=%.9f" % (cerrada, serie))
for m in (10 ** 4, 10 ** 6, 10 ** 8, 10 ** 10):
    print("m=%.0e  E|S_4|^2/m^3 = %.9f" % (m, formula11(m, 4) / m ** 3))
