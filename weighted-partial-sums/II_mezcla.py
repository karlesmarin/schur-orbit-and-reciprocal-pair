# -*- coding: utf-8 -*-
r"""II_mezcla.py -- el factor ~1.41 entre a(m) y la gaussiana: es una mezcla de escalas?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Fijados los signos de los primos <= sqrt(m), cada n <= m tiene a lo sumo un primo > sqrt(m) y con
exponente 1, asi que
    S = R + sum_{sqrt m < p <= m} A_p eps(p),   A_p = p * sum_{k <= m/p} k eps(k)
con R y los A_p funcion SOLO de los signos pequenos, y los eps(p) grandes independientes.
Prediccion de mezcla (TCL condicional, reticulo de paso 2):
    a_mix = 2^{pi(m)} * E_small[ 2/sqrt(2 pi V) * exp(-R^2/(2V)) ],   V = sum A_p^2,
con la esperanza EXACTA sobre las 2^{pi(sqrt m)} asignaciones pequenas.  Se compara a/a_mix con
a/a_G (gaussiana global).  Si a/a_mix ~ 1 donde a/a_G ~ 1.4, el exceso es la mezcla.
"""
import math
import sys
from itertools import product

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


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


def libre(d):
    return all(v == 1 for v in fact(d).values())


def eps_val(c, e):
    v = 1
    for p, k in fact(c).items():
        if k % 2 and e.get(p, 1) == -1:
            v = -v
    return v


reales = {}
for line in open("II_exceso_mod4_OUT.txt", encoding="utf-8"):
    if line.startswith("m="):
        import re
        mm = re.match(r"m=\s*(\d+)\s+a=\s*(\d+)", line)
        if mm:
            reales[int(mm.group(1))] = int(mm.group(2))

print("   m     a      a/a_G   a/a_mix")
for m in sorted(reales):
    if m < 20:
        continue
    r = math.isqrt(m)
    peq = primos(r)
    grandes = [p for p in primos(m) if p > r]
    acc = 0.0
    for sg in product((1, -1), repeat=len(peq)):
        e = dict(zip(peq, sg))
        R = sum(n * eps_val(n, e) for n in range(1, m + 1)
                if all(p <= r for p in fact(n)))
        V = 0
        for p in grandes:
            Ap = p * sum(k * eps_val(k, e) for k in range(1, m // p + 1))
            V += Ap * Ap
        acc += 2 / math.sqrt(2 * math.pi * V) * math.exp(-R * R / (2 * V))
    amix = 2 ** len(grandes) * acc          # 2^{pi(m)} * (1/2^{#peq}) * suma
    mu = sum(j * j for j in range(1, r + 1))
    s2 = sum((d * sum(j * j for j in range(1, math.isqrt(m // d) + 1))) ** 2
             for d in range(2, m + 1) if libre(d))
    aG = 2 ** len(primos(m)) * 2 / math.sqrt(2 * math.pi * s2) * math.exp(-mu * mu / (2 * s2))
    print("%4d %8d   %5.2f    %5.2f" % (m, reales[m], reales[m] / aG, reales[m] / amix))
