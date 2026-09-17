# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# Identidades de columnas del lema de duplicacion (L_{2r} = 2 L_r) y de la subfamilia m = Br, Br-1:
#   s_{2r}(2c)      = 2 s_r(c)
#   s_{2r}(nu(c))   = 2 (s_r(c) - s_r(2^{-1} c))      nu(c) = representante impar mod 2r de c mod r
#   h_{2r}(nu(c))   = s_r(2^{-1} c)
#   h_{2r}(2c)      = s_r(2c) - s_r(c)
#   primer momento de m = B r o B r - 1 (B impar) sobre Z/2r  =  B h_{2r}
# con s_n(0) = 0, s_n(x) = 2x - n (1 <= x < n), h_{2r}(c) = c (0<c<r), 0 (c=0,r), c-2r (r<c<2r).
# Uso: python duplicacion.py [RMAX]
import sys


def s(n, x):
    x %= n
    return 0 if x == 0 else 2 * x - n


def h2(r, x):
    x %= 2 * r
    if x == 0 or x == r:
        return 0
    return x if x < r else x - 2 * r


def S(m, n, c):
    c %= n
    return sum(j for j in range(1, m + 1) if j % n == c) - sum(j for j in range(1, m + 1) if (j + c) % n == 0)


def main():
    RMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 45
    fallos = 0
    medidas = 0
    for r in range(3, RMAX + 1, 2):
        inv2 = pow(2, -1, r)
        for c in range(r):
            nu = c if c % 2 == 1 else c + r  # representante impar mod 2r
            nu %= 2 * r
            pruebas = [
                (s(2 * r, 2 * c), 2 * s(r, c)),
                (s(2 * r, nu), 2 * (s(r, c) - s(r, inv2 * c))),
                (h2(r, nu), s(r, inv2 * c)),
                (h2(r, 2 * c), s(r, 2 * c) - s(r, c)),
            ]
            for a, b in pruebas:
                medidas += 1
                fallos += a != b
        for B in (1, 3, 5):
            for m in (B * r, B * r - 1):
                if m < 1:
                    continue
                for c in range(2 * r):
                    medidas += 1
                    fallos += S(m, 2 * r, c) != B * h2(r, c)
    print("r impar en [3, %d]: medidas %d  fallos %d" % (RMAX, medidas, fallos))


if __name__ == "__main__":
    main()
