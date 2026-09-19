# -*- coding: utf-8 -*-
r"""II_exceso_mod4.py -- el exceso de a(m) sobre la gaussiana en m = 0 (mod 4): se sostiene hasta 120?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

a(m) exacto: enumeracion de los primos <= m/2 y programacion dinamica sobre los de (m/2, m].
Gaussiana con los momentos exactos del modelo de signos aleatorios en los primos:
  mu = sum_{j^2<=m} j^2,   sigma^2 = sum_{1<d<=m libre de cuadrados} W_d^2,  W_d = d sum_{j<=sqrt(m/d)} j^2,
  a_G = 2^{pi(m)} * 2/sqrt(2 pi sigma^2) * exp(-mu^2/(2 sigma^2))   (reticulo de paso 2).
Se imprime a/a_G separado por m mod 4 y por si m es primo, y la media geometrica de cada grupo.
"""
import math
import sys
from itertools import product

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120


def primos(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p ** .5) + 1))]


def fact(c):
    out, d = {}, 2
    while d * d <= c:
        while c % d == 0:
            out[d] = out.get(d, 0) + 1
            c //= d
        d += 1
    if c > 1:
        out[c] = out.get(c, 0) + 1
    return out


def a(m):
    ps = primos(m)
    L = [p for p in ps if 2 * p > m]
    peq = [p for p in ps if 2 * p <= m]
    idx = {p: i for i, p in enumerate(peq)}
    fijos = []
    for c in range(1, m + 1):
        f = fact(c)
        if any(p in L for p in f):
            continue
        mask = 0
        for p, k in f.items():
            if k % 2:
                mask |= 1 << idx[p]
        fijos.append((c, mask))
    dp = {0: 1}
    for p in L:
        nd = {}
        for s, k in dp.items():
            nd[s + p] = nd.get(s + p, 0) + k
            nd[s - p] = nd.get(s - p, 0) + k
        dp = nd
    tot = 0
    for bits in range(1 << len(peq)):
        base = 0
        for c, mask in fijos:
            base += -c if bin(bits & mask).count("1") % 2 else c
        tot += dp.get(-base, 0)
    return tot


def libre(d):
    return all(v == 1 for v in fact(d).values())


filas = []
for m in range(3, MMAX + 1):
    if m % 4 not in (0, 3):
        continue
    am = a(m)
    mu = sum(j * j for j in range(1, math.isqrt(m) + 1))
    s2 = sum((d * sum(j * j for j in range(1, math.isqrt(m // d) + 1))) ** 2
             for d in range(2, m + 1) if libre(d))
    aG = 2 ** len(primos(m)) * 2 / math.sqrt(2 * math.pi * s2) * math.exp(-mu * mu / (2 * s2))
    esp = all(m % d for d in range(2, m)) and m > 1
    filas.append((m, am, aG, am / aG, m % 4, esp))
    print("m=%3d a=%7d aG=%9.1f  a/aG=%5.2f  m%%4=%d  %s" % (m, am, aG, am / aG, m % 4, "primo" if esp else ""))


def gm(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float("nan")


for nombre, sel in (("m=0 mod 4", lambda r: r[4] == 0), ("m=3 mod 4, m primo", lambda r: r[4] == 3 and r[5]),
                    ("m=3 mod 4, m compuesto", lambda r: r[4] == 3 and not r[5])):
    xs = [r[3] for r in filas if sel(r) and r[0] >= 20]
    print("%-24s  media geometrica de a/aG (m>=20): %.3f  en %d filas" % (nombre, gm(xs), len(xs)))
