# -*- coding: utf-8 -*-
r"""II_s_menos_y_131.py -- dos afirmaciones de la nota II, por suma directa y por la recursion (D).

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

(e) s = -1  =>  S(q^2 m) = q^2 S(m)   (q = 11, 19, 23; m = 1, 38, 75, ... hasta 3000).
(131) q = 131, s = +1: 27 + 43 (131 + ... + 131^t) es cero para t = 0..5 (directo hasta 10^6).
Uso: python II_s_menos_y_131.py
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def tablas(q):
    L = [0] + [-1] * (q - 1)
    for r in range(1, q):
        L[r * r % q] = 1
    C = [0] * q; T = [0] * q; c = t = 0
    for r in range(q):
        c += L[r]; t += r * L[r]; C[r] = c; T[r] = t
    return L, C, T


def S_dir(q, s, m):
    L, C, T = tablas(q); tot = 0
    for n in range(1, m + 1):
        x, sg = n, 1
        while x % q == 0:
            x //= q; sg *= s
        tot += n * sg * L[x % q]
    return tot


def S_rec(q, s, m, C, T):
    if m == 0:
        return 0
    k, r = divmod(m, q)
    return s * q * S_rec(q, s, k, C, T) + k * (T[q - 1] + q * C[r]) + T[r]


ok = tot = 0
for q in (11, 19, 23):
    L, C, T = tablas(q)
    for m in range(1, 3001, 37):
        tot += 1; ok += (S_rec(q, -1, q * q * m, C, T) == q * q * S_dir(q, -1, m))
print("(e) s=-1: S(q^2 m) = q^2 S(m) en %d de %d casos (q = 11, 19, 23)" % (ok, tot))
q = 131
L, C, T = tablas(q)
res = []
for t in range(6):
    m = 27 + 43 * sum(q ** j for j in range(1, t + 1))
    res.append((m, S_dir(q, 1, m) if m < 10 ** 6 else S_rec(q, 1, m, C, T)))
print("(131) m_t y S(m_t):", res)
