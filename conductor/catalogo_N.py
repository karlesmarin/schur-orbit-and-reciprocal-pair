# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# CATALOGO DIRIGIDO: e contra N, el cociente del criterio de soporte (C2):
#   tipo "m"    : m = N n        tipo "m+1"  : m + 1 = N n        tipo "2m+1" : 2m + 1 = N n  (n impar)
# Conjetura en observacion (NO pre-registrada todavia): e = e0 . p^{v_p(N)}, e0 = e con p no | N.
# Perfil hasta K = min(E, KMAX); e visible si hay al menos 2 niveles totales por encima.
#
# Uso: python catalogo_N.py [KMAX QMAX]
#
# Authors: Carles Marin, Claude (AI assistant).
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perfil_W import dims_B  # noqa: E402
from reparto_p2 import phi  # noqa: E402


def vp(x, p):
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    QMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 420
    res = defaultdict(list)
    for p in (2, 3, 5):
        for n in (5, 7, 9, 11, 13):
            if n % p == 0:
                continue
            f, x = 1, p % n
            while x not in (1, n - 1):
                x = (x * p) % n
                f += 1
            full = phi(n) // 2
            for v in range(1, 9):
                q = p ** v * n
                if q > QMAX:
                    break
                E = phi(p ** v)
                K = min(E, KMAX)
                if K < 3 or K * phi(n) > 150:
                    continue
                for N in range(1, 60):
                    for tipo, m in (("m", N * n), ("m+1", N * n - 1), ("2m+1", (N * n - 1) // 2 if (N * n) % 2 else None)):
                        if m is None or m < 1 or not (2 * m < q) or m == q // 2 or (q % 2 == 0 and m == q // 2 - 1):
                            continue
                        dB = dims_B(q, p, m, n, K)
                        W = [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, K)]
                        e = None
                        for kk in range(1, K - 2):
                            if all(W[i] == full for i in range(kk, K)):
                                e = kk
                                break
                        res[(p, tipo, n, f)].append((vp(N, p), N, v, E, m, q, e, W))
    for key in sorted(res):
        p, tipo, n, f = key
        porv = defaultdict(set)
        for (vN, N, v, E, m, q, e, W) in res[key]:
            porv[vN].add(e)
        print("p=%d %-5s q'=%2d f=%d :  e por v_p(N): %s" % (p, tipo, n, f, {k: sorted(s, key=str) for k, s in sorted(porv.items())}))
    print()
    print("detalle de las casillas con v_p(N) >= 1:")
    for key in sorted(res):
        for (vN, N, v, E, m, q, e, W) in res[key]:
            if vN >= 1:
                print("   p=%d %-5s q'=%2d  N=%2d v_p(N)=%d  v=%d E=%2d m=%3d q=%3d  e=%s  W=%s" % (key[0], key[1], key[2], N, vN, v, E, m, q, e, W))


if __name__ == "__main__":
    main()
