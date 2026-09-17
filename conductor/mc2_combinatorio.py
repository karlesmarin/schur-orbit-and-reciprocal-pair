# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# MC2 combinatorio (prediccion escrita antes de los datos, paso 1).
#   Para m <= M, q <= Q, 2m < q, p | q con f >= 2:  N_m D-estable  =>  N_m G-estable ?
# Imprime cuantas ternas se miraron, el reparto (D-estable / G-estable) para que se vea que el test no
# es vacio, y TODOS los contraejemplos (casillas asesinas de H2', con g y el numero de multiconjuntos
# distintos que D3 predice como dimension de la imagen).
#
# Uso: python mc2_combinatorio.py [M Q]
#
# Authors: Carles Marin, Claude (AI assistant).
import sys
from math import gcd
from collections import Counter


def factor_p(q, p):
    v = 0
    while q % p == 0:
        q //= p
        v += 1
    return v, q


def primos(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def clase(j, qp):
    j %= qp
    return min(j, (qp - j) % qp)


def multiconjunto(m, qp):
    c = Counter(clase(j, qp) for j in range(1, m + 1))
    return tuple(c.get(k, 0) for k in range(qp // 2 + 1))


def actuar(N, u, qp):
    out = [0] * len(N)
    for k, n in enumerate(N):
        if n:
            out[clase(k * u, qp)] += n
    return tuple(out)


def main():
    M, Q = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (12, 400)
    ternas = 0
    cnt = Counter()
    asesinas = []
    for q in range(3, Q + 1):
        for p in primos(q):
            v, qp = factor_p(q, p)
            if qp <= 2:
                continue
            unidades = [u for u in range(1, qp) if gcd(u, qp) == 1]
            # f = orden de p en (Z/q')*/+-1
            f, x = 1, p % qp
            while x not in (1, qp - 1):
                x = (x * p) % qp
                f += 1
            if f < 2:
                continue
            g = (len(unidades) // (2 if qp > 2 else 1)) // f
            for m in range(1, M + 1):
                if not (2 * m < q):
                    continue
                ternas += 1
                N = multiconjunto(m, qp)
                d_est = actuar(N, p, qp) == N
                orbita = set(actuar(N, u, qp) for u in unidades)
                g_est = len(orbita) == 1
                cnt[(d_est, g_est)] += 1
                if d_est and not g_est:
                    asesinas.append((m, q, p, v, qp, f, g, len(orbita)))
    print("ternas (m,q,p) con f>=2 miradas: %d   (M=%d, Q=%d)" % (ternas, M, Q))
    print("  D-estable y G-estable     %d" % cnt[(True, True)])
    print("  D-estable, NO G-estable   %d   <- contraejemplos de MC2" % cnt[(True, False)])
    print("  no D-estable              %d" % (cnt[(False, True)] + cnt[(False, False)]))
    print("  (G-estable sin D-estable: %d -- tiene que ser 0, D esta dentro de G)" % cnt[(False, True)])
    print("MC2: %s" % ("SIN contraejemplos en el rango" if not asesinas else "*** MUERE: %d ternas ***" % len(asesinas)))
    for (m, q, p, v, qp, f, g, k) in sorted(asesinas, key=lambda t: (t[1] if False else 0, t[0], t[1]))[:60]:
        print("  m=%2d q=%3d p=%d v=%d q'=%3d f=%d g=%d  multiconjuntos distintos=%d  grado=%d"
              % (m, q, p, v, qp, f, g, k, sum(1 for u in range(1, q) if gcd(u, q) == 1) // 2))


if __name__ == "__main__":
    main()
