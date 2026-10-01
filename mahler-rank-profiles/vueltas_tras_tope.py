# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: las vueltas DESPUES DEL TOPE p-1. Punto de partida: (5,75619), lambda = 6, primera vuelta
#      presente del piso 3 en k = 8). Caracteres cuadraticos chi = (-n/.), -n fundamental, p no divide n phi(n):
#        - lambda exacta por la serie de Iwasawa (filas_negativas.serie_iwasawa, rama <x>, N = 4, i < p^3);
#        - piso 3: vueltas extendidas J_k = p^2 - 1 + (k-1) p (p-1), k = 1..KMAX (J_k <= L), cadena P/A completa;
#        - piso 2 (tope): primera vuelta presente.
#      Solo se imprimen los casos con lambda >= p - 1 (los demas ya estan medidos: vueltas = lambda).
#      Pistas a contrastar (literatura): Childress 1989 tiene un tope lambda <= p; Saikia-Barman lo suben a 2p; Satoh lo
#      quita. Digitos base p de lambda frente a la primera vuelta presente.
# FECHA: 2-oct-2026
# Uso (docker, /work, con -e OPENBLAS_NUM_THREADS=4): sage -python vueltas_tras_tope.py p nmax [KMAX] [COLS]
import sys
import time

sys.path.insert(0, '/work')
import filas_negativas as FN  # noqa: E402
import mahler_masas as MM  # noqa: E402
from sage.all import ZZ, euler_phi, kronecker_character  # noqa: E402

p, nmax = int(sys.argv[1]), int(sys.argv[2])
KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 3 * p
COLS = int(sys.argv[4]) if len(sys.argv) > 4 else 3  # piso de columnas para leer las vueltas del piso 3
N = 4
Jk3 = [p * p - 1 + (k - 1) * p * (p - 1) for k in range(1, KMAX + 1)]
L3 = Jk3[-1]
Jk2 = [p - 1 + (k - 1) * (p - 1) for k in range(1, p)]
t0 = time.time()
tot = altos = 0
casos = [75619] + list(range(3, nmax + 1, 4))
for n in casos:
    if not ZZ(-n).is_fundamental_discriminant() or n % p == 0 or euler_phi(n) % p == 0:
        continue
    chi = kronecker_character(-n)
    ex = [-1 if v == 0 else (0 if v == 1 else 1) for v in chi.values()]
    mu, g, f = FN.medida(ex, 2, n, p, N)
    ls = FN.lam(FN.serie_iwasawa(mu, p, N, 0, 4 * p))
    tot += 1
    if ls is None or ls < p - 1:
        continue
    altos += 1
    _, piv2, c2, cp = MM.cuadratico(n, p, 2, Jk2[-1])
    v2 = ''.join('P' if J in piv2 else 'A' for J in Jk2)
    # Filas J >= p^a no son intrinsecas en el piso a (Teorema F solo vale para J < p^a). Con COLS = 4 las vueltas extendidas del piso 3
    # (J_k = p^2-1 + (k-1)p(p-1), J <= L3) se leen en la matriz de p^4 columnas, donde J < p^4 si lo son.
    _, piv3, _, _ = MM.cuadratico(n, p, COLS, L3)
    v3 = ''.join('P' if J in piv3 else 'A' for J in Jk3)
    k3 = v3.find('P') + 1 if 'P' in v3 else None
    print("n=%d: lambda %d (base %d: %s) | chi(2)=%s chi(p)=%s | piso2 %s | piso3 %s | primera P piso3 k=%s (k-1=%s)" % (
        n, ls, p, ZZ(ls).digits(p), c2, cp, v2, v3, k3, (k3 - 1) if k3 else None), flush=True)
assert tot > 0
print("p=%d n<=%d: %d caracteres, %d con lambda >= p-1 (%.0f s)" % (p, nmax, tot, altos, time.time() - t0), flush=True)
