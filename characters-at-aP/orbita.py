# -*- coding: utf-8 -*-
"""LA REFORMULACION POR ORBITAS.

R_x(X) = sum_{alpha in Phi} X^{(x,alpha)}  en Z[X]/(X^q-1)  es invariante bajo:
  (i)  W            -- porque permuta Phi;
  (ii) qQ^v         -- porque (x+q g, alpha) = (x,alpha) + q(g,alpha), y (g,alpha) es entero;
  (iii) el NUCLEO de P/qP -> Hom(Q, Z/q), que en tipo C es {0, (q/2)(1,...,1)} con q par,
        porque ese vector da 0 en e_i-e_j y q en e_i+e_j y en 2e_i.

Sea W_q^+ el grupo generado por los tres.  CONJETURA:

        chi_lambda(a_P^k) != 0   <=>   ell in W_q^+ . rho ,

y entonces R_ell = R_rho es INMEDIATO por invariancia, en cualquier tipo.
"""
import os
import sys
from itertools import permutations, product
from math import gcd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # antes: un scratchpad temporal (9-sep-2026)
from verifica_informe import chi, raices_C


def dominantes(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for r in rec(k - 1, v):
                yield (v,) + r
    return rec(m, top)


def en_orbita(ell, rho, q, m, con_nucleo):
    """existe w en W(C_m) (permutacion con signos) y s en el nucleo con ell = w(rho) + s (mod q)?"""
    desplaz = [[0] * m]
    if con_nucleo and q % 2 == 0:
        desplaz.append([q // 2] * m)
    for s in desplaz:
        objetivo = [(ell[i] - s[i]) % q for i in range(m)]
        for sigma in permutations(range(m)):
            for eps in product((1, -1), repeat=m):
                if all((eps[i] * rho[sigma[i]] - objetivo[i]) % q == 0 for i in range(m)):
                    return True
    return False


print("=" * 100)
print("chi != 0   <=>   ell in W_q . rho   (sin el nucleo)  /  in W_q^+ . rho  (con el nucleo)")
print("=" * 100)
tot = {"sin": [0, 0, 0], "con": [0, 0, 0]}   # [vivos capturados, vivos perdidos, muertos capturados]
for m in (2, 3, 4, 5):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    for k in range(1, t):
        d = gcd(k, t); q = t // d
        if q < 2:
            continue
        filas = {}
        for modo in ("sin", "con"):
            cap = perd = falso = 0
            for lam in dominantes(m, 6 if m <= 3 else 4):
                ell = [lam[i] + rho[i] for i in range(m)]
                dentro = en_orbita(ell, rho, q, m, modo == "con")
                vivo = chi(lam, m, k) != 0
                if vivo and dentro:
                    cap += 1
                elif vivo and not dentro:
                    perd += 1
                elif (not vivo) and dentro:
                    falso += 1
            filas[modo] = (cap, perd, falso)
            for i, v in enumerate((cap, perd, falso)):
                tot[modo][i] += v
        (c1, p1, f1), (c2, p2, f2) = filas["sin"], filas["con"]
        marca = ""
        if p2 == 0 and f2 == 0:
            marca = "  <-- con nucleo: EXACTO"
        print("  C_%d k=%-2d d=%-2d q=%-2d | sin nucleo: captura %4d, PIERDE %4d, falsos %4d"
              " | con nucleo: captura %4d, pierde %4d, falsos %4d%s"
              % (m, k, d, q, c1, p1, f1, c2, p2, f2, marca))
print()
print("TOTAL sin nucleo : captura %d, PIERDE %d, falsos %d" % tuple(tot["sin"]))
print("TOTAL con nucleo : captura %d, pierde %d, falsos %d" % tuple(tot["con"]))
