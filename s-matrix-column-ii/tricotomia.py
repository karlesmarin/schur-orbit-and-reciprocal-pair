# -*- coding: utf-8 -*-
r"""tricotomia.py -- la tricotomia, contra el censo.

LO QUE SE COMPRUEBA: para h = 1, FUERA de la union de los cuatro lugares (M2,M3) = 0, (M2,M4) = 0, (M2,M5) = 0, (M3,M5) = 0, se tiene

        delta - 1 = [M2 = 0] + [M3 = 0],

porque: M2 M3 != 0 da <2,3>; solo M2 = 0 con M3,M4,M5 != 0 da <3,4,5>; solo M3 = 0 con M2,M5 != 0
da <2,5>.  Aqui se comprueba en todas las orbitas del censo, separando las que caen en los lugares
excluidos (que son las unicas donde se permite otra cosa).

TAMBIEN: su version de la formula de columnas maximales con el factor h/(p-1) en vez de n,

        z(n,p) = sum_{h>1, h|p-1, n mod h in {0,1}} [h/(p-1)] [ C(e,d) - Z_{p,e}(d) ],

que cuenta ORBITAS con delta = 0 (la nuestra, con el factor n, cuenta columnas).

    python tricotomia.py [PMAX] [NMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from fractions import Fraction
from itertools import combinations
from math import comb

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import centrados, es_primo, estabilizador, momento, semigrupo, s_exacto   # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 23
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8


def Z(p, e, d):
    """subconjuntos de mu_e de tamano d con suma 0."""
    g = pow(min(x for x in range(2, p) if pow(x, e, p) == 1 and
               all(pow(x, e // q, p) != 1 for q in range(2, e + 1) if e % q == 0)), 1, p) if e > 1 else 1
    mu = [pow(g, i, p) for i in range(e)] if e > 1 else [1]
    return sum(1 for V in combinations(mu, d) if sum(V) % p == 0)


print("(1) TRICOTOMIA  delta - 1 = [M2=0] + [M3=0]  fuera de los cuatro lugares")
print("%-4s %-3s %-8s %-10s %-10s %s" % ("p", "n", "orbitas", "h=1 fuera", "fallos", "en los lugares excluidos"))
fallos_tot = fuera_tot = 0
for p in range(7, PMAX + 1):
    if not es_primo(p):
        continue
    for n in range(3, min(NMAX, p - 2) + 1):
        fuera = dentro = fallos = 0
        for T in centrados(p, n):
            if len(estabilizador(T, p)) != 1:
                continue
            M = {r: momento(T, r, p) % p for r in (2, 3, 4, 5)}
            lugares = [(2, 3), (2, 4), (2, 5), (3, 5)]
            if any(M[a] == 0 and M[b] == 0 for a, b in lugares):
                dentro += 1
                continue
            fuera += 1
            d = semigrupo(s_exacto(T, p))["genero"]
            if d - 1 != (1 if M[2] == 0 else 0) + (1 if M[3] == 0 else 0):
                fallos += 1
        fallos_tot += fallos
        fuera_tot += fuera
        print("%-4d %-3d %-8d %-10d %-10d %d" % (p, n, len(centrados(p, n)), fuera, fallos, dentro))
        sys.stdout.flush()
print("orbitas h=1 fuera de los lugares: %d ; fallos de la tricotomia: %d" % (fuera_tot, fallos_tot))

print("")
print("(2) z(n,p) = orbitas con delta = 0, con el factor h/(p-1)")
print("%-4s %-3s %-10s %-14s %s" % ("p", "n", "medido", "formula", "ok"))
fallos_z = casos_z = 0
for p in range(7, min(PMAX, 23) + 1):
    if not es_primo(p):
        continue
    for n in range(3, min(NMAX, p - 2) + 1):
        medido = sum(1 for T in centrados(p, n) if semigrupo(s_exacto(T, p))["genero"] == 0)
        tot = Fraction(0)
        for h in range(2, p):
            if (p - 1) % h == 0 and n % h in (0, 1):
                e, d = (p - 1) // h, n // h
                if d <= e:
                    tot += Fraction(h, p - 1) * (comb(e, d) - Z(p, e, d))
        ok = (tot == medido)
        casos_z += 1
        fallos_z += 0 if ok else 1
        print("%-4d %-3d %-10d %-14s %s" % (p, n, medido, str(tot), "si" if ok else "NO"))
        sys.stdout.flush()
print("casos: %d ; fallos: %d" % (casos_z, fallos_z))

print("")
print("VEREDICTO: tricotomia %s ; z(n,p) %s"
      % ("SOBREVIVE" if not fallos_tot else "FALLA", "SOBREVIVE" if not fallos_z else "FALLA"))
