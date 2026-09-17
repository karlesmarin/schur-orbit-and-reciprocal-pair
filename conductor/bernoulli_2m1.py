# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# SONDA tipo 2m+1 (exploratoria).  2m+1 = N n, n primo, S(c) = N . chat(c) (residuo con signo).
# Derivacion:  S^(chi) = N n (chibar(2) - 1) B_{1,chi}  para chi impar
#              (usando sum_{a<n/2} chi(a) = (chibar(2) - 2) B_{1,chi}; la primera version decia chi(2) y el
#              control numerico la tumbo -- ver la salida vieja en el commit).
# Prediccion (semisimple p no | |G|, p impar no | n N, p no | h^-):
#     total - W_1 = #{chi impar : chi(2) = 1} = t/2 si -1 no esta en <2>, 0 si esta   (t = [(Z/n)^* : <2>])
# Controles: (a) la identidad de media suma se comprueba numericamente; (b) casos p | h^- se listan aparte.
#
# Authors: Carles Marin, Claude (AI assistant).
import cmath
import sys
import os
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bernoulli_W1 import rango_mod_p, primitiva, h_menos_numerico, factores, es_primo  # noqa: E402


def chat(c, n):
    c %= n
    return c if c <= n // 2 else c - n


def W1_media(n, p):
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    filas = [[chat(u * c, n) % p for c in range(1, n)] for u in units]
    return rango_mod_p(filas, p)


def identidad(n):
    g = primitiva(n)
    log = {}
    x = 1
    for i in range(n - 1):
        log[x] = i
        x = x * g % n
    peor = 0.0
    for k in range(1, n - 1, 2):
        chi = lambda a: cmath.exp(2j * cmath.pi * k * log[a % n] / (n - 1))
        B1 = sum(a * chi(a) for a in range(1, n)) / n
        izq = sum(chi(a) for a in range(1, (n + 1) // 2))
        # CORREGIDO: es chi(2) CONJUGADO (= chi(2)^-1).  La primera version con chi(2) daba error ~ 6..158:
        # el recuento de abajo acertaba igual porque chibar(2) = 1 <=> chi(2) = 1.
        der = (chi(2).conjugate() - 2) * B1
        izq2 = sum(chat(c, n) * chi(c) for c in range(1, n))
        der2 = n * (chi(2).conjugate() - 1) * B1
        peor = max(peor, abs(izq - der), abs(izq2 - der2))
    return peor


def main():
    print("identidades de media suma (error maximo numerico):")
    for n in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        print("   n=%2d  %.2e" % (n, identidad(n)))
    tab = {}
    contra = []
    for n in [x for x in range(5, 90) if es_primo(x)]:
        h = h_menos_numerico(n)
        full = (n - 1) // 2
        ord2 = 1
        y = 2 % n
        while y != 1:
            y = y * 2 % n
            ord2 += 1
        menos1_en_2 = pow(2, ord2 // 2, n) == n - 1 if ord2 % 2 == 0 else False
        t = (n - 1) // ord2
        pred_def = 0 if menos1_en_2 else t // 2
        for p in [x for x in range(3, 40) if es_primo(x)]:
            if p == n or full % p == 0:
                continue
            w = W1_media(n, p)
            if h % p == 0:
                tab["p | h^- (aparte)"] = tab.get("p | h^- (aparte)", 0) + 1
                print("   aparte n=%d p=%d h^-=%d W_1=%d de %d, deficit chi(2)=1 predicho %d" % (n, p, h, w, full, pred_def))
                continue
            ok = (full - w == pred_def)
            clave = ("ok" if ok else "FALLA", "pred_def>0" if pred_def else "pred_def=0")
            tab[clave] = tab.get(clave, 0) + 1
            if not ok:
                contra.append((n, p, w, full, pred_def, ord2, menos1_en_2))
    print("resultado:", tab)
    print("contraejemplos (n, p, W_1, total, deficit predicho, ord 2, -1 en <2>):", contra[:20])


if __name__ == "__main__":
    main()
