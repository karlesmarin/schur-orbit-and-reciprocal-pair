# -*- coding: utf-8 -*-
r"""II_ceros_heegner.py -- la cadena de ceros de numero de clases uno, y los q = 3 mod 8 con h > 1.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Afirmacion (de (D) con r = 0):  q = 3 mod 4, q > 3, s = chi(q) = +1, h(-q) = 1  =>
   S(q) = T_q + q s = q (s - h) = 0,  y  S(z) = 0, chi(z+1) = 1  =>  S(q(z+1)) = q (z+1)(chi(z+1) - h) = 0;
   como z_i = 0 mod q, chi(z_i + 1) = 1 siempre: z_0 = q, z_{i+1} = q (z_i + 1) son todos ceros.
Comprobacion DIRECTA (suma explicita, sin usar (D)) de los z_i <= ZMAX para q en {7, 11, 19, 43, 67, 163}.
Control: con s = -1 la misma cadena NO es de ceros (S(q) = -2q).
Pregunta abierta (medida): q = 3 mod 8, h(-q) > 1, s = +1: numero de ceros <= MMAX.
Uso: python II_ceros_heegner.py [ZMAX] [MMAX]
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ZMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30_000_000
MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 10_000_000


def legendre(q):
    t = np.full(q, -1, dtype=np.int64); t[0] = 0
    for r in range(1, q):
        t[r * r % q] = 1
    return t


def S_de(q, s, N):
    L = legendre(q)
    x = np.arange(N + 1, dtype=np.int64)
    eps = np.ones(N + 1, dtype=np.int64)
    while True:
        div = (x % q == 0) & (x > 0)
        if not div.any():
            break
        x[div] //= q; eps[div] *= s
    eps = eps * L[x % q]; eps[0] = 0
    return np.cumsum(np.arange(N + 1, dtype=np.int64) * eps)


def h_menos(q):
    L = legendre(q)
    return -int(sum(r * int(L[r]) for r in range(1, q))) // q


for q in (7, 11, 19, 43, 67, 163):
    assert h_menos(q) == 1
    z = [q]
    while q * (z[-1] + 1) <= ZMAX:
        z.append(q * (z[-1] + 1))
    S = S_de(q, 1, z[-1])
    Sm = S_de(q, -1, z[-1])
    print("q=%3d h=1  cadena %s  S=0 (s=+1): %s   control s=-1, S(z_i): %s"
          % (q, z, [int(S[x]) == 0 for x in z], [int(Sm[x]) for x in z[:3]]))
    sys.stdout.flush()


def es_primo(n):
    return n > 1 and all(n % d for d in range(2, int(n ** .5) + 1))


print("q = 3 mod 8, h > 1, s = +1: ceros <= %d" % MMAX)
for q in [p for p in range(11, 400) if es_primo(p) and p % 8 == 3]:
    h = h_menos(q)
    if h == 1:
        continue
    S = S_de(q, 1, MMAX)
    z = np.nonzero(S[1:] == 0)[0] + 1
    print("  q=%3d h=%2d  ceros=%d  %s" % (q, h, len(z), list(map(int, z[:6]))))
    sys.stdout.flush()
