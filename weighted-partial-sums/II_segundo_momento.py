# -*- coding: utf-8 -*-
r"""II_segundo_momento.py -- la cota del segundo momento del Lema 2 de la nota II, medida exacta.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Signos primos independientes y equiprobables, extendidos completamente multiplicativos:
E eps(n) eps(k) = 1 si nk es un cuadrado, 0 si no.  Con n = d a^2, k = d b^2 (d libre de cuadrados),
  E S(m)^2 = sum_{d libre de cuadrados <= m} d^2 (sum_{a <= sqrt(m/d)} a^2)^2.
Afirmacion a juzgar:  E S(m)^2 <= m^3 (log m + 13) / 9  para todo m >= 5000 (y de hecho m >= 1?).
E R_Q^2 = E S^2 - sum_{q in Q} q^2 <= E S^2, asi que basta S.
Calculo: criba de libres de cuadrados hasta 10^8, sumas prefijas de d^2 (enteros exactos en int64
por tramos -> se usa object/float128 no: se acumula en Python int por bloques de k).
Se agrupa por k = floor(sqrt(m/d)): d en (m/(k+1)^2, m/k^2].
Se mide el cociente E S^2 / (m^3 (log m + 13)/9) en: todo m <= 20000, y 4000 puntos log-espaciados
hasta 10^8 mas los extremos de cada bloque del certificado.  Tambien la constante asintotica
E S^2 / (m^3 log m) frente a 2/(3 pi^2).
"""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
N = 10 ** 8
sf = np.ones(N + 1, dtype=bool)
sf[0] = False
for p in range(2, int(N ** .5) + 1):
    sf[p * p::p * p] = False
d = np.arange(N + 1, dtype=np.float64)
# suma prefija de d^2 sobre libres de cuadrados; float64 con error relativo ~1e-16*N: suficiente
# para un cociente (se informa el margen minimo, que debe ser grande frente a 1e-8)
pref = np.cumsum(np.where(sf, d * d, 0.0))
del d


def ES2(m):
    tot = 0.0
    k = 1
    while k * k <= m:
        hi = m // (k * k)
        lo = m // ((k + 1) * (k + 1))
        s2 = k * (k + 1) * (2 * k + 1) // 6
        tot += (pref[hi] - pref[lo]) * s2 * s2
        k += 1
    return tot


def cota(m):
    return m ** 3 * (math.log(m) + 13) / 9


puntos = set(range(3, 20001))
puntos |= {int(round(10 ** (math.log10(20000) + i * (8 - math.log10(20000)) / 4000))) for i in range(4001)}
a = 5000
while a <= N - 1:
    b = min(11 * a // 10, N - 1)
    puntos |= {a, b}
    a = b + 1
puntos = sorted(x for x in puntos if x <= N - 1)
peor = (0, None)
peor5000 = (0, None)
for m in puntos:
    r = ES2(m) / cota(m)
    if r > peor[0]:
        peor = (r, m)
    if m >= 5000 and r > peor5000[0]:
        peor5000 = (r, m)
for m in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, N - 1):
    e = ES2(m)
    print("m=%9d  E S^2/cota=%.4f  E S^2/(m^3 log m)=%.5f  (2/(3pi^2)=%.5f)"
          % (m, e / cota(m), e / (m ** 3 * math.log(m)), 2 / (3 * math.pi ** 2)))
print("puntos medidos: %d; max cociente (m>=3): %.4f en m=%d; max (m>=5000): %.4f en m=%d  (debe ser < 1)"
      % (len(puntos), peor[0], peor[1], peor5000[0], peor5000[1]))
