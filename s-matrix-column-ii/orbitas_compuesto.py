# -*- coding: utf-8 -*-
r"""orbitas_compuesto.py -- el reticulo de orbitas para m COMPUESTO: que se rompe y que aguanta.

LO QUE TENEMOS PARA m PRIMO (periodos_gauss.py):
  (i) el estabilizador H de un semisistema tiene orden IMPAR;
  (ii) un semisistema H-invariante es una eleccion de signos sobre los periodos de Gauss;
  (iii) para h > 1 todos los periodos son CERO, porque la suma de un subgrupo de F_p^x de orden >= 2
       es cero, luego A(h) = 2^{(p-1)/(2h)};
  (iv) Mobius da el reticulo entero.

QUE PUEDE ROMPERSE CON m COMPUESTO, y este guion lo mide en vez de suponerlo:
  * (Z/m)^x NO es ciclico, asi que puede haber subgrupos de orden PAR que no contengan a -1; si los
    hay, (i) es falso fuera de los primos;
  * el conjunto {1,...,m-1} tiene elementos NO UNIDADES, que las unidades permutan pero cuyas
    orbitas tienen otro tamano;
  * la suma de un subgrupo de (Z/m)^x NO tiene por que ser 0 mod m, asi que (iii) puede fallar y
    entonces A(H) ya no es una potencia de dos.

    python orbitas_compuesto.py [MMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import sys
from itertools import product
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 45


def es_primo(x):
    return x > 1 and all(x % d for d in range(2, int(x ** 0.5) + 1))


print("%-5s %-5s %-8s %-8s %-30s %s"
      % ("m", "tipo", "N(m)", "orbitas", "estabilizadores (orden: cuantas orbitas)", "hay h par?"))
hay_par_global = []
for m in range(7, MMAX + 1, 2):
    n = (m - 1) // 2
    U = [a for a in range(1, m) if gcd(a, m) == 1]
    sols = []
    for elec in product(*[(t, m - t) for t in range(1, n + 1)]):
        T = frozenset(elec)
        if sum(T) % m == 0:
            sols.append(T)
    conj = set(sols)
    vistos, orbitas = set(), []
    for T in sols:
        if T in vistos:
            continue
        orb = {frozenset((a * t) % m for t in T) for a in U} & conj
        vistos |= orb
        H = [a for a in U if frozenset((a * t) % m for t in T) == T]
        orbitas.append((len(orb), len(H), sorted(H)))
    rep = {}
    for tam, h, Hs in orbitas:
        rep[h] = rep.get(h, 0) + 1
    pares = [h for h in rep if h % 2 == 0]
    if pares:
        hay_par_global.append((m, pares))
    print("%-5s %-5s %-8d %-8d %-30s %s"
          % (m, "primo" if es_primo(m) else "comp.", len(sols), len(orbitas),
             " ".join("h=%d:%d" % (h, rep[h]) for h in sorted(rep)),
             "SI %s" % pares if pares else "no"))
    sys.stdout.flush()

print("")
if hay_par_global:
    print("ROTO FUERA DE LOS PRIMOS: hay estabilizadores de orden PAR en %s" % hay_par_global)
    print("La demostracion de (i) usa que -1 genera el unico subgrupo de orden 2 --- cierto en")
    print("F_p^x, que es ciclico, y FALSO en (Z/m)^x cuando no lo es.")
else:
    print("ningun estabilizador par en el rango: (i) aguanta tambien para m compuesto aqui.")

print("")
print("SUMA DE LOS SUBGRUPOS: para m primo vale 0 en cuanto el orden es >= 2.  Aqui, medido:")
for m in range(9, min(MMAX, 45) + 1, 2):
    if es_primo(m):
        continue
    U = [a for a in range(1, m) if gcd(a, m) == 1]
    subs = []
    for a in U:
        S, x = set(), 1
        while x not in S:
            S.add(x)
            x = x * a % m
        if len(S) >= 2:
            subs.append((len(S), sum(S) % m))
    nz = sorted({(o, s) for o, s in subs if s})
    print("   m=%-4d subgrupos ciclicos con suma NO nula: %s" % (m, nz if nz else "ninguno"))
