# -*- coding: utf-8 -*-
r"""II_cola_analitica.py -- las desigualdades de la prueba de A para M >= 10^8, en nuestra version.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Prueba:  L = log M, Q = primos de (M/2, M], N = |Q|.
 (a) Dusart (6.6): pi(x) >= x/(ln x - 1) (x >= 5393), pi(x) <= x/(ln x - 1.1) (x >= 60184)
     => N >= M/(L-1) - (M/2)/(L - log 2 - 1.1) >= M/(2L)   (se comprueba g(L) > 1/2 para L >= 8.66).
 (b) X = n primos de Q, n = floor(N/8) >= M/(17L).  BW Thm 4.23 con l = 4M (4 max X <= 4M = l) pide
     l <= n^2/(12 log(4 max X/n)), implicado por n^2 >= 48 M log(4M/n), implicado por
     M >= 13872 L^2 log(68 L).  Da una PA de l terminos, paso s <= 4M/n, terminos <= 6 l M/n = 24M^2/n.
 (c) Y = s-1 primos de Q \ X (s <= 4M/n <= 68L): sus sumas cubren Z/s.  S_{X u Y} contiene el intervalo
     [alpha, alpha + w], alpha <= 24M^2/n + 4M^2/n = 28M^2/n, w = (l-1)s - Sigma_Y >= M.
 (d) Resto Z: cada primo <= M <= w; Sigma_Z >= P_Q/2 => S_Q contiene [alpha, P_Q - alpha];
     sumas con signo: toda la paridad de P_Q en [-D, D], D = P_Q - 2 alpha, K := 2 alpha <= 56M^2/n.
 (e) Segundo momento: |R| <= B = (M^{3/2}/3) sqrt(L + 13).  P_Q >= N M/2 >= M^2/(4L).
     B/P_Q <= (4/3) L sqrt(L+13)/sqrt(M),   K/P_Q <= 224 L/n <= 3808 L^2/M.
Se comprueba: g(L) > 1/2; la hipotesis de (b); las dos cotas y su suma < 1 en M = 10^8 y en una malla
logaritmica hasta 10^40 (y que son decrecientes en la malla); ademas, con pi EXACTO en M = 10^8
(criba), N, n y las cotas sin redondeo por Dusart.
"""
import math, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def g(L):
    return L / (L - 1) - L / (2 * (L - math.log(2) - 1.1))


def cotas(M):
    L = math.log(M)
    hip = M >= 13872 * L * L * math.log(68 * L)
    b = (4 / 3) * L * math.sqrt(L + 13) / math.sqrt(M)
    k = 3808 * L * L / M
    return L, hip, b, k


print("g(L) - 1/2 en L = 8.7, 18.42, 50, 200: %s" % [round(g(L) - .5, 5) for L in (8.7, 18.42, 50, 200)])
malla = [10 ** (8 + i * 0.05) for i in range(641)]
prev = None; ok = True
for M in malla:
    L, hip, b, k = cotas(M)
    ok &= hip and b + k < 1 and g(L) > .5
    if prev is not None:
        ok &= b <= prev[0] + 1e-15 and k <= prev[1] + 1e-15
    prev = (b, k)
L, hip, b, k = cotas(1e8)
print("M=1e8: L=%.3f  hipotesis BW (M >= 13872 L^2 log 68L = %.3e): %s  B/P <= %.5f  K/P <= %.5f  suma %.5f"
      % (L, 13872 * L * L * math.log(68 * L), hip, b, k, b + k))
print("malla 1e8..1e40 (641 puntos): hipotesis, suma < 1, g > 1/2 y cotas decrecientes: %s" % ok)
# umbral de la hipotesis
lo, hi = 1e6, 1e8
for _ in range(100):
    mid = math.sqrt(lo * hi)
    if cotas(mid)[1]:
        hi = mid
    else:
        lo = mid
print("la hipotesis de (b) vale desde M ~ %.3e" % hi)
# comprobacion con pi exacto en 10^8
N8 = 10 ** 8
es = np.ones(N8 + 1, dtype=bool); es[:2] = False
for p in range(2, int(N8 ** .5) + 1):
    if es[p]:
        es[p * p::p] = False
idx = np.nonzero(es[N8 // 2 + 1:])[0] + N8 // 2 + 1
N = len(idx); PQ = int(idx.sum()); n = N // 8
M = N8; L = math.log(M)
B = M ** 1.5 / 3 * math.sqrt(L + 13); K = 56 * M * M / n
print("M=1e8 exacto: N=%d (cota M/2L=%.0f)  n=%d (cota M/17L=%.0f)  P_Q=%.4e (cota M^2/4L=%.4e)"
      % (N, M / (2 * L), n, M / (17 * L), PQ, M * M / (4 * L)))
print("  n^2=%.3e >= 48 M log(4M/n)=%.3e: %s   (B+K)/P_Q = %.5f" %
      (n * n, 48 * M * math.log(4 * M / n), n * n >= 48 * M * math.log(4 * M / n), (B + K) / PQ))
