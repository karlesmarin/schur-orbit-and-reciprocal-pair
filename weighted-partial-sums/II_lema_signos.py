# -*- coding: utf-8 -*-
r"""II_lema_signos.py -- el lema de sumas con signo que falta en la linea A de II, medido.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

L = primos de (m/2, m], P = sum L.  Las sumas sum_{p in L} +-p tienen todas la paridad de P y viven en
[-P, P].  K(m) = el menor K tal que TODOS los enteros de esa paridad en [-P+K, P-K] se alcanzan.
Equivalente: las sumas de subconjuntos de L cubren [K/2, P-K/2] (s = (P - suma con signo)/2).
Se mide K(m) y K/P; el lema util para A es K pequeno frente a P (la ventana [-P+K, P-K] casi entera).
Bitset con enteros de Python.
Uso: python II_lema_signos.py MMAX
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000


def primos_hasta(n):
    cr = bytearray([1]) * (n + 1); cr[0:2] = b"\x00\x00"
    for i in range(2, int(n ** .5) + 1):
        if cr[i]:
            cr[i * i::i] = bytearray(len(cr[i * i::i]))
    return [i for i in range(n + 1) if cr[i]]


PR = primos_hasta(MMAX)
peor = 0
for m in list(range(10, 200, 10)) + list(range(200, MMAX + 1, 100)):
    L = [p for p in PR if m / 2 < p <= m]
    P = sum(L)
    bits = 1
    for p in L:
        bits |= bits << p
    # menor s0 tal que todas las sumas de subconjunto s0..P-s0 existen
    s0 = 0
    while s0 <= P // 2:
        ok = True
        # comprobar el rango completo es caro; buscamos el primer hueco desde abajo
        break
    faltan = [s for s in range(0, P // 2 + 1) if not (bits >> s) & 1]
    s0 = (max(faltan) + 1) if faltan else 0
    K = 2 * s0
    peor = max(peor, K / P if P else 0)
    if m % 500 == 0 or m < 100:
        print("m=%5d |L|=%4d P=%9d  K=%6d  K/P=%.4f  (huecos de subconjunto: %d, el mayor en %d)"
              % (m, len(L), P, K, K / P if P else 0, len(faltan), max(faltan) if faltan else -1))
print("max K/P con m >= 200:", round(peor, 4))
