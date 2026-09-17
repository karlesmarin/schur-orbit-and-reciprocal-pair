# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# TRES LEMAS LOCALES SOBRE LAS CAPAS, contra los datos (hipotesis del Teorema 5: q' no en {1,2,3,4,6}, E >= 2,
# q' | h, h+1, h+2; d = phi(q')/2):
#   (M)  monotonia general:  W_s > 0  =>  W_{i+s} >= W_i  (cualquier numero de primos sobre p)
#   (B)  desigualdad dual:   W_{i+j} < d  =>  W_i + W_j <= d
#   (G)  simetria Gorenstein: con e el exponente del reticulo,
#          [ W_i + W_{e-1-i} = d para todo i < e ]  <=>  a == b   (a = v_p[O:A], b = v_p[A:c])
#        y ademas la desigualdad 2 sum_{i<e} W_i <= e d.
# Y la CERTIFICACION del catalogo sintetico (catalogo_e.py) con (M): si W_1 > 0, e = primera capa llena, sin
# pedir un solo primo sobre p.
#
# Fuentes: perfil_W_OUT.txt (filas con reticulo: W, a, b, forma (f,e)) y perfil_W.dims_B (catalogo).
# Uso: python capas_locales.py [NMAX QMAX KMAX]
import os
import re
import sys
from collections import Counter

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from perfil_W import dims_B  # noqa: E402
from reparto_p2 import phi  # noqa: E402

FILA = re.compile(r"^\s+ok\s+m=\s*(\d+) q=\s*(\d+) p=(\d+) v=(\d+) q'=\s*(\d+) E=\s*(\d+) forma=(\[.*?\]) a=(\d+) b=(\d+)\s+W=(\[.*?\])")


def lemas(W, d, K=None):
    """Devuelve (medidas_M, fallos_M, medidas_B, fallos_B) sobre el perfil W (capas 0..len-1)."""
    K = len(W)
    mM = fM = mB = fB = 0
    for s in range(1, K):
        if W[s] > 0:
            for i in range(0, K - s):
                mM += 1
                fM += W[i + s] < W[i]
    for i in range(0, K):
        for j in range(i, K - i):
            if i + j < K and W[i + j] < d:
                mB += 1
                fB += W[i] + W[j] > d
    return mM, fM, mB, fB


def reticulo():
    c = Counter()
    ejemplos = []
    for ln in open(os.path.join(AQUI, "perfil_W_OUT.txt"), encoding="utf-8"):
        mt = FILA.match(ln)
        if not mt:
            continue
        m, q, p, v, n, E = (int(mt.group(k)) for k in range(1, 7))
        forma = eval(mt.group(7))
        a, b = int(mt.group(8)), int(mt.group(9))
        W = eval(mt.group(10))
        if n in (1, 2, 3, 4, 6) or E < 2:
            c["fuera de las hipotesis del Teorema 5"] += 1
            continue
        d = phi(n) // 2
        e = forma[0][1]
        if any(fe[1] != e for fe in forma):
            c["exponentes distintos (contradice el Lema 3)"] += 1
        mM, fM, mB, fB = lemas(W, d)
        c["M medidas"] += mM
        c["M fallos"] += fM
        c["B medidas"] += mB
        c["B fallos"] += fB
        if len(W) >= e:  # capas 0..e-1 calculadas
            c["G filas con capas < e completas"] += 1
            sim = all(W[i] + W[e - 1 - i] == d for i in range(e))
            c["G desigualdad 2 sum W <= e d falla"] += 2 * sum(W[:e]) > e * d
            if sim == (a == b):
                c["G coincide"] += 1
            else:
                c["G NO coincide"] += 1
                ejemplos.append((m, q, p, W, e, a, b))
            c["G Gorenstein (a == b)"] += a == b
        if W[1:] and W[1] > 0:
            llenas = [i for i in range(len(W)) if W[i] == d]
            if llenas:
                c["e == primera capa llena (W_1>0)"] += llenas[0] == e
                c["e != primera capa llena (W_1>0)"] += llenas[0] != e
            c["W_1>0 con g>1"] += len(forma) > 1
    return c, ejemplos


def catalogo(NMAX, QMAX, KMAX):
    c = Counter()
    e_cert = Counter()
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
                tipos = [0, (n - 1) // 2, n - 1] if n % 2 else [0, n // 2 - 1, n // 2, n - 1]
                for k in range(0, 40):
                    for r in tipos:
                        m = k * n + r
                        if m < 1 or not (2 * m < q) or m == q // 2 or (q % 2 == 0 and m == q // 2 - 1):
                            continue
                        dB = dims_B(q, p, m, n, K)
                        W = [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, K)]
                        mM, fM, mB, fB = lemas(W, full)
                        c["M medidas"] += mM
                        c["M fallos"] += fM
                        c["B medidas"] += mB
                        c["B fallos"] += fB
                        if 2 * W[1] > full and K >= 3:
                            c["2W_1 > d medidas"] += 1
                            c["2W_1 > d y W_2 < d (falla)"] += W[2] < full
                        if not (1 <= W[1] < full) or K < 4:
                            continue
                        c["decidibles (1<=W_1<d, K>=4)"] += 1
                        llenas = [i for i in range(1, K) if W[i] == full]
                        if llenas:
                            c["certificados por (M), cualquier g"] += 1
                            c["  de ellos con g>1"] += g > 1
                            e_cert[llenas[0]] += 1
                        else:
                            c["sin capa llena calculada"] += 1
    return c, e_cert


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    QMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 260
    KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    c, ej = reticulo()
    print("RETICULO (perfil_W_OUT.txt):")
    for k in sorted(c):
        print("   %-50s %d" % (k, c[k]))
    print("   ejemplos G no coincide:", ej[:5])
    c2, e_cert = catalogo(NMAX, QMAX, KMAX)
    print("CATALOGO SINTETICO (dims_B, mismo rango que catalogo_e.py):")
    for k in sorted(c2):
        print("   %-50s %d" % (k, c2[k]))
    print("   exponente certificado:", dict(e_cert))


if __name__ == "__main__":
    main()
