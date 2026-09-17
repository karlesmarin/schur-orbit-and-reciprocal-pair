# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# CATALOGO SINTETICO del exponente e y del perfil W (instrumento: perfil_W.dims_B, 184/184 con el
# reticulo).  Para m en Est(n), A != Z, 2m < q:
#   W_0..W_{K-1} con K = min(E, KMAX);  e visible = menor k <= K-2 con W_i total para i en [k, K-1]
#   (se exige al menos UN nivel total por encima de k para no confundir e con el borde).
# Informes:
#   E2   1 <= W_1 < total  =>  e = 2 ?   (cuenta confirmaciones, contraejemplos, y no decidibles)
#   tabla de e por (p, tipo de r, v_p(escala), f) -- la escala de S: 2k+1 (r=(n-1)/2), k (r=0),
#   k+1 (r=n-1), y para n par la que toque (se marca aparte).
#
# Uso: python catalogo_e.py [NMAX QMAX KMAX]
#
# Authors: Carles Marin, Claude (AI assistant).
import os
import sys
from math import gcd
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perfil_W import dims_B  # noqa: E402
from reparto_p2 import phi  # noqa: E402


def vp(x, p):
    if x == 0:
        return 99
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    QMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 260
    KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    e2 = Counter()
    e2_contra = []
    tabla = defaultdict(Counter)
    filas = []
    for p in (2, 3, 5, 7):
        for n in range(5, NMAX + 1):
            if n % p == 0 or n == 6:
                continue
            f, x = 1, p % n
            while x not in (1, n - 1):
                x = (x * p) % n
                f += 1
            full = phi(n) // 2
            g = full // f
            for v in range(1, 8):
                q = p ** v * n
                if q > QMAX:
                    break
                E = phi(p ** v)
                K = min(E, KMAX)
                if K * phi(n) > 170 or K < 2:
                    continue
                if n % 2:
                    tipos = {0: "r=0", (n - 1) // 2: "r=(n-1)/2", n - 1: "r=n-1"}
                else:
                    tipos = {0: "r=0", n // 2 - 1: "r=n/2-1", n // 2: "r=n/2", n - 1: "r=n-1"}
                for k in range(0, 40):
                    for r, tipo in tipos.items():
                        m = k * n + r
                        if m < 1 or not (2 * m < q) or m == q // 2 or (q % 2 == 0 and m == q // 2 - 1):
                            continue
                        escala = {"r=0": k, "r=(n-1)/2": 2 * k + 1, "r=n-1": k + 1}.get(tipo, None)
                        dB = dims_B(q, p, m, n, K)
                        W = [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, K)]
                        e = None
                        for kk in range(1, K - 1):
                            if all(W[i] == full for i in range(kk, K)):
                                e = kk
                                break
                        if e is None and W[K - 1] == full and K >= 2 and all(W[i] == full for i in range(K - 1, K)):
                            e = "borde"
                        filas.append((p, n, f, g, v, E, tipo, k, m, q, W, e))
                        if 1 <= W[1] < full:
                            if K < 4:
                                e2["no decidible (K<4)"] += 1
                            elif e == 2:
                                e2["confirma"] += 1
                            else:
                                e2["CONTRA"] += 1
                                e2_contra.append((m, q, p, v, n, f, g, E, W, e))
                        vs = vp(escala, p) if escala is not None else "n par"
                        tabla[(p, tipo, vs)][e if e is not None else "?"] += 1
    print("E2  1 <= W_1 < total => e = 2 :", dict(e2))
    for t in e2_contra[:25]:
        print("   E2 CONTRA m=%d q=%d p=%d v=%d q'=%d f=%d g=%d E=%d W=%s e=%s" % t)
    print("e por (p, tipo de r, v_p(escala)):")
    for key in sorted(tabla, key=str):
        print("   %-28s %s" % (key, dict(tabla[key])))
    print("perfiles con W_1 = 0 (primeros 40):")
    c = 0
    for t in filas:
        if t[10][1] == 0 and c < 40:
            print("   p=%d q'=%d f=%d g=%d v=%d E=%d %-10s k=%d m=%d q=%d W=%s e=%s" % t)
            c += 1


if __name__ == "__main__":
    main()
