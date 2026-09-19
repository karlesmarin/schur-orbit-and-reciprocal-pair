# -*- coding: utf-8 -*-
r"""II_ruta_A.py -- la ruta de prueba de la linea A, ejecutada en m grande.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Ruta: (i) signos en los primos <= sqrt(m) fijos (base: eps(2)=eps(3)=-1, resto +1); primos de
(sqrt(m), m/2] con peso c_p = p*T(floor(m/p)), signos por equilibrado voraz (de mayor a menor |c_p|,
el signo que acerca la suma parcial a 0), partiendo de S_liso; (ii) los primos de (m/2, m] (peso 1)
tienen que absorber el resto R: basta |R| <= P - K con K ~ 4m (II_lema_signos_OUT.txt) y la paridad.
Se imprime |R|, P - 4m y el cociente; si el cociente es < 1 con holgura, la ruta cierra en ese m.
Tambien max|c_p| y S_liso/m^2, que es lo que una prueba tiene que acotar.
Uso: python II_ruta_A.py
"""
import sys, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def criba_menor_factor(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** .5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


MS = [10 ** 4, 3 * 10 ** 4, 10 ** 5, 3 * 10 ** 5, 10 ** 6]
N = max(MS)
spf = criba_menor_factor(N)
for m in MS:
    raiz = math.isqrt(m)
    signo = {2: -1, 3: -1}
    # eps en todo c <= m con los primos grandes a +1 (para la parte lisa y los T(k))
    eps = [0] * (m + 1); eps[1] = 1
    mayor = [1] * (m + 1)
    for c in range(2, m + 1):
        p = spf[c]
        eps[c] = (signo.get(p, 1) if p <= raiz else 1) * eps[c // p]
        mayor[c] = max(p, mayor[c // p])
    T = [0] * (raiz + 2)
    for k in range(1, raiz + 2):
        T[k] = T[k - 1] + k * eps[k]
    S_liso = sum(c * eps[c] for c in range(1, m + 1) if mayor[c] <= raiz)
    medios = [p for p in range(raiz + 1, m // 2 + 1) if spf[p] == p]
    pesos = sorted((p * T[m // p] for p in medios), key=abs, reverse=True)
    R = S_liso
    for w in pesos:
        R = R - w if abs(R - w) < abs(R + w) else R + w
    grandes = [p for p in range(m // 2 + 1, m + 1) if spf[p] == p]
    P = sum(grandes)
    ventana = P - 4 * m
    print("m=%8d  S_liso/m^2=%+.4f  max|c_p|/m^1.5=%.3f  |R|=%d  P-4m=%d  |R|/(P-4m)=%.2e  paridad R+P par: %s"
          % (m, S_liso / m ** 2, max(abs(w) for w in pesos) / m ** 1.5, abs(R), ventana, abs(R) / ventana,
             (R + P) % 2 == 0))
    sys.stdout.flush()
