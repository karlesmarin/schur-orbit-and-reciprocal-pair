# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# SONDA (no pre-registrada; exploratoria): W_1 no total  <=>  p | h^-(Q(zeta_n)) ?
# Para m = n - 1 (o m = kn con p no | k), S(c) = N.(2c - n) y su transformada en un caracter impar chi
# es 2nN.B_{1,chi}.  W_1 = dim F_p[G].S (teorema del reparto, G = (Z/n)*/+-1).  Aqui n primo.
# Tabla de h^- de Washington (Introduction to Cyclotomic Fields, tabla §3): valores citados de memoria
# por Claude -- SE COMPRUEBAN abajo calculando h^- = Q.w.prod(-B_{1,chi}/2) numericamente.
#
# Authors: Carles Marin, Claude (AI assistant).
import cmath
from math import gcd

H_MENOS_MEMORIA = {23: 3, 29: 8, 31: 9, 37: 37, 41: 121, 43: 211, 47: 695, 53: 4889, 59: 41241, 61: 76301,
                   67: 853513, 71: 3882809, 73: 11957417, 79: 100146415, 83: 838216959}
for n in (3, 5, 7, 11, 13, 17, 19):
    H_MENOS_MEMORIA[n] = 1


def primitiva(n):
    for g in range(2, n):
        if all(pow(g, (n - 1) // r, n) != 1 for r in set(factores(n - 1))):
            return g


def factores(x):
    out, d = [], 2
    while d * d <= x:
        while x % d == 0:
            out.append(d)
            x //= d
        d += 1
    if x > 1:
        out.append(x)
    return out


def h_menos_numerico(n):
    # n primo impar: h^- = Q w prod_{chi impar} (-B_{1,chi}/2), Q = 1, w = 2n
    g = primitiva(n)
    log = {}
    x = 1
    for i in range(n - 1):
        log[x] = i
        x = x * g % n
    prod = 1
    for k in range(1, n - 1, 2):  # caracteres impares: chi(g) = e^{2 pi i k/(n-1)}, k impar
        B1 = sum(a * cmath.exp(2j * cmath.pi * k * log[a] / (n - 1)) for a in range(1, n)) / n
        prod *= -B1 / 2
    return round((2 * n * prod).real)


def rango_mod_p(filas, p):
    filas = [r[:] for r in filas]
    rk, ncol = 0, len(filas[0])
    for col in range(ncol):
        piv = next((i for i in range(rk, len(filas)) if filas[i][col] % p), None)
        if piv is None:
            continue
        filas[rk], filas[piv] = filas[piv], filas[rk]
        inv = pow(filas[rk][col], -1, p)
        filas[rk] = [(v * inv) % p for v in filas[rk]]
        for i in range(len(filas)):
            if i != rk and filas[i][col] % p:
                t = filas[i][col]
                filas[i] = [(a - t * b) % p for a, b in zip(filas[i], filas[rk])]
        rk += 1
    return rk


def W1(n, p):
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    filas = [[(2 * ((u * c) % n) - n) % p for c in range(1, n)] for u in units]
    return rango_mod_p(filas, p)


def es_primo(x):
    return x > 1 and all(x % d for d in range(2, int(x ** 0.5) + 1))


def main():
    print("comprobacion de la tabla de h^- de memoria contra el calculo numerico:")
    mal = 0
    for n, h in sorted(H_MENOS_MEMORIA.items()):
        hn = h_menos_numerico(n)
        if hn != h:
            mal += 1
        print("   n=%2d  memoria %10d  numerico %10d  %s" % (n, h, hn, "ok" if hn == h else "*** DISTINTO ***"))
    print("   tabla: %d discrepancias" % mal)
    print()
    print("W_1 no total  <=>  p | h^-   (n primo, p primo impar != n, p <= 13 y los primos de h^-)")
    tab = {}
    contra = []
    for n in sorted(H_MENOS_MEMORIA):
        h = h_menos_numerico(n)
        full = (n - 1) // 2
        ps = set(q for q in range(3, 14) if es_primo(q)) | set(q for q in factores(h) if q != 2 and q < 5000)
        for p in sorted(ps):
            if p == n:
                continue
            w = W1(n, p)
            no_total = w < full
            divide = (h % p == 0)
            semis = (full % p != 0)
            clave = (no_total, divide, "semisimple" if semis else "p | |G|")
            tab[clave] = tab.get(clave, 0) + 1
            if no_total != divide:
                contra.append((n, p, h, w, full, semis))
            if divide or no_total:
                print("   n=%2d p=%4d  h^-=%10d  W_1=%2d de %2d  %s" % (n, p, h, w, full, "semisimple" if semis else "p | |G|"))
    print("tabla (W_1 no total, p | h^-, estrato):", tab)
    print("contraejemplos:", contra)


if __name__ == "__main__":
    main()
