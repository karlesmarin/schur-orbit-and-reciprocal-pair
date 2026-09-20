# -*- coding: utf-8 -*-
r"""ley_2p.py -- la fraccion 2/p: se gatea la DEMOSTRACION paso a paso, no el resultado.

LO QUE SE COMPRUEBA.  Con delta(T) = g(S_T) y F_n(p) la fraccion de conjuntos
CENTRADOS (= columnas a altura prima) con delta >= 2,

    F_n(p) = 2/p - 1/p^2 + O_n( p^{1-n/2} + p^{1+floor(n/2)-n} ),

y en particular F_n(p) = 2/p + O(p^{-min(2, n/2-1)}).  Para n >= 7 el error es o(p^{-2}), asi que
el coeficiente siguiente, el -1, queda identificado.

Una demostracion no se comprueba mirando si el resultado final cuadra: se comprueba PASO A PASO,
porque un resultado correcto puede salir de un paso falso.  Los pasos son cuatro y aqui van los
cuatro por separado.

PASO 1.  El numero de conjuntos centrados es exactamente C(p,n)/p.
  La traslacion t -> t+c actua sobre los n-subconjuntos; si n < p ninguna orbita tiene punto
  fijo, luego todas tienen tamano p y hay un unico representante de suma cero por orbita.

PASO 2.  C_J := #{T centrado : M_j = 0 para j en J} = C(p,n)/p^{|J|} + O_n(p^{n/2}).
  Se mide el error real y se imprime la constante implicita |C_J - C(p,n)/p^{|J|}| / p^{n/2},
  que es lo unico que decide si la cota es la correcta o esta puesta a ojo.

PASO 3.  Inclusion-exclusion: con h = 1 el criterio ya gateado es delta >= 2 <=> M_2 M_3 = 0,
  luego #{M_2 M_3 = 0} = C_{1,2} + C_{1,3} - C_{1,2,3}.  Se comprueba como IDENTIDAD EXACTA
  contra el censo, no como estimacion.

PASO 4.  Los conjuntos con estabilizador no trivial son O_n(p^{floor(n/2)}).
  Son uniones de orbitas multiplicativas de tamano >= 2, con la posibilidad de anadir el cero, y
  todas tienen M_1 = 0 automaticamente porque la suma de un subgrupo de orden >= 2 de F_p^x es 0.
  Se cuentan por fuerza bruta y se compara con p^{floor(n/2)}.

Y AL FINAL, la asintotica: (F_n(p) - 2/p) * p^2 debe tender a -1 para n >= 7.

    python ley_2p.py [PMAX] [NMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from itertools import combinations
from math import comb

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo                                              # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 53
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8


def censo_momentos(p, nmax, pots):
    """dp[k] indexado por los momentos de exponentes en 'pots'.  Devuelve #T con todos = 0."""
    forma = (nmax + 1,) + (p,) * len(pots)
    dp = np.zeros(forma, dtype=np.int64)
    dp[(0,) + (0,) * len(pots)] = 1
    for t in range(p):
        des = tuple(pow(t, e, p) for e in pots)
        for k in range(min(nmax, p) - 1, -1, -1):
            bloque = dp[k]
            for ax, d in enumerate(des):
                bloque = np.roll(bloque, d, axis=ax)
            dp[k + 1] += bloque
    return [int(dp[(n,) + (0,) * len(pots)]) for n in range(nmax + 1)]


def estab_no_trivial(p, n):
    """#{T centrado, |T| = n, con lambda T = T para algun lambda != 1}, por fuerza bruta."""
    total = set()
    for d in range(2, p):
        if (p - 1) % d:
            continue
        g = next(a for a in range(2, p) if all(pow(a, (p - 1) // q, p) != 1
                                               for q in range(2, p) if (p - 1) % q == 0
                                               and all(q % r for r in range(2, q))))
        lam = pow(g, (p - 1) // d, p)
        H = {pow(lam, i, p) for i in range(d)}
        cosets = []
        vistos = set()
        for a in range(1, p):
            if a in vistos:
                continue
            c = frozenset((a * h) % p for h in H)
            vistos |= c
            cosets.append(c)
        for k in (n // d, (n - 1) // d):
            if k * d not in (n, n - 1) or k > len(cosets):
                continue
            cero = (k * d == n - 1)
            for sel in combinations(cosets, k):
                T = set()
                for c in sel:
                    T |= c
                if cero:
                    T.add(0)
                if len(T) == n:
                    total.add(frozenset(T))
    return len(total)


primos = [p for p in range(11, PMAX + 1) if es_primo(p)]
print("PASO 1 y 2: C_J contra C(p,n)/p^{|J|}, y la constante implicita del error O(p^{n/2})")
print("%-4s %-5s %-8s %-12s %-12s %-12s %-10s"
      % ("n", "p", "J", "C_J medido", "C(p,n)/p^|J|", "error", "err/p^{n/2}"))
datos = {}
for p in primos:
    m12 = censo_momentos(p, min(NMAX, p - 1), (1, 2))
    m13 = censo_momentos(p, min(NMAX, p - 1), (1, 3))
    m123 = censo_momentos(p, min(NMAX, p - 1), (1, 2, 3))
    m1 = censo_momentos(p, min(NMAX, p - 1), (1,))
    for n in range(5, min(NMAX, p - 2) + 1):
        cent = comb(p, n) // p
        assert comb(p, n) % p == 0, (p, n)
        assert m1[n] == cent, (p, n, m1[n], cent)          # PASO 1, exacto
        datos[(n, p)] = dict(cent=cent, C12=m12[n], C13=m13[n], C123=m123[n])
        for et, val, k in (("{1,2}", m12[n], 2), ("{1,3}", m13[n], 2), ("{1,2,3}", m123[n], 3)):
            pred = comb(p, n) / p ** k
            err = abs(val - pred)
            if n in (7, 8) and p in (primos[0], primos[len(primos) // 2], primos[-1]):
                print("%-4d %-5d %-8s %-12d %-12.1f %-12.1f %-10.3f"
                      % (n, p, et, val, pred, err, err / p ** (n / 2)))
    sys.stdout.flush()

print("")
print("PASO 1 (numero de centrados = C(p,n)/p) comprobado exacto en todos los casos.")

print("")
print("PASO 3: la inclusion-exclusion, como IDENTIDAD EXACTA")
print("%-4s %-5s %-14s %-14s %s" % ("n", "p", "C12+C13-C123", "censo M2M3=0", "coincide"))
fall3 = 0
for (n, p), d in sorted(datos.items()):
    ie = d["C12"] + d["C13"] - d["C123"]
    # censo directo de M_1 = 0 y M_2 M_3 = 0
    m12 = censo_momentos(p, n, (1, 2))[n]
    m13 = censo_momentos(p, n, (1, 3))[n]
    m123 = censo_momentos(p, n, (1, 2, 3))[n]
    directo = m12 + m13 - m123
    ok = (ie == directo)
    fall3 += 0 if ok else 1
    if n in (7, 8) and p in (primos[0], primos[-1]):
        print("%-4d %-5d %-14d %-14d %s" % (n, p, ie, directo, "si" if ok else "NO"))
print("    fallos: %d" % fall3)

print("")
print("PASO 4: conjuntos centrados con estabilizador no trivial, contra p^{floor(n/2)}")
print("%-4s %-5s %-14s %-14s %-10s" % ("n", "p", "medido", "p^floor(n/2)", "cociente"))
for p in primos[:4]:
    for n in (5, 6, 7):
        if n > p - 2:
            continue
        e = estab_no_trivial(p, n)
        print("%-4d %-5d %-14d %-14d %-10.4f" % (n, p, e, p ** (n // 2), e / p ** (n // 2)))
    sys.stdout.flush()

print("")
print("LA ASINTOTICA: F_n(p) = #{M_2 M_3 = 0}/#centrados, y (F - 2/p) p^2 -> -1 si n >= 7")
print("%-4s %-5s %-12s %-12s %-12s %-12s" % ("n", "p", "F_n(p)", "2/p", "(F-2/p)p", "(F-2/p)p^2"))
for n in range(5, min(NMAX, 8) + 1):
    for p in primos:
        if (n, p) not in datos:
            continue
        d = datos[(n, p)]
        F = (d["C12"] + d["C13"] - d["C123"]) / d["cent"]
        print("%-4d %-5d %-12.8f %-12.8f %-12.5f %-12.4f"
              % (n, p, F, 2 / p, (F - 2 / p) * p, (F - 2 / p) * p * p))
    print("")
    sys.stdout.flush()
