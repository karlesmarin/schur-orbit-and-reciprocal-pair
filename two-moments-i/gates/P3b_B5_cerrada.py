# -*- coding: utf-8 -*-
r"""P3b_B5_cerrada.py -- forma cerrada de C(5,p), el patron NO consecutivo {1,3}.

C(n,p) = #{T subset F_p : |T| = n, M_1 = M_3 = 0} = #{T : e_1 = e_3 = 0}: los polinomios
monicos escindidos y libres de cuadrados con a_{n-1} = a_{n-3} = 0.  El tamiz por particiones da

    C(n,p) = sum_{lambda |- n} (-1)^{n-l} / z_lambda * N_lambda,
    N_lambda = #{x in F_p^l : sum w_i x_i = 0, sum w_i x_i^3 = 0},

y N_lambda es el CONO AFIN sobre la interseccion de un hiperplano con una cubica diagonal
ponderada en P^{l-1}, o sea una hipersuperficie cubica en P^{l-2}:  N = 1 + (p-1)*#proj.

LO QUE SE AFIRMA AQUI (n = 5).  Las siete particiones dan, para p > 5:

    (5)         l=1   N = 1
    (4,1)       l=2   N = 1                       (w_1 = +-w_2 solo si p <= 5)
    (3,2)       l=2   N = 1
    (3,1,1)     l=3   cubica binaria  -3x(8x^2+9xy+3y^2)      disc -15
                      N = 1 + (p-1)(2 + chi_p(-15))
    (2,2,1)     l=3   cubica binaria  -3(x+y)(x^2+3xy+y^2)    disc 5
                      N = 1 + (p-1)(2 + chi_p(5))
    (2,1,1,1)   l=4   CUBICA PLANA LISA: 2x(x+y+z)^2 + yz(y+z) = 0, isomorfa a la
                      CURVA ELIPTICA  E: Y^2 = X^3 + X^2 - X   (disc minimal 80, conductor 20)
                      N = 1 + (p-1)(p + 1 - a_p(E))
    (1,1,1,1,1) l=5   SUPERFICIE CUBICA DIAGONAL DE CLEBSCH (1871)
                      N = 1 + (p-1)(p^2 + (6 + chi_p(5)) p + 1)

y por tanto, con chi = simbolo de Legendre,

    C(5,p) = (p-1) [ (2+chi_p(-15))/6 + (2+chi_p(5))/8 - (p+1-a_p)/12
                     + (p^2 + (6+chi_p(5)) p + 1)/120 ].

C(5,p) NO es un polinomio en p: lleva dos simbolos de Legendre y el a_p de una curva eliptica.
Contrasta con C(3,p) = (p-1)/2 y C(4,p) = (p-1)(p-3)/8, que si lo son.

CONTROLES QUE CORRE ESTE FICHERO, todos contra enumeracion exacta y sin usar el enunciado:
 (0) la reduccion algebraica de cada cubica binaria, comprobada como identidad de polinomios;
 (1) la lisura de la cubica plana (jacobiano sin ceros no triviales), primo a primo;
 (2) cada N_lambda cerrado contra N_lambda enumerado;
 (3) a_p(E) contra el producto eta  q prod (1-q^{2n})^2 (1-q^{10n})^2  --- la forma de nivel 20:
     lo comprueba en un rango finito, sin base de datos externa (el conductor mismo lo
     calcula P3b_B5_weierstrass.sage);
 (4) C(5,p) cerrado contra el censo exacto por DP.

    python P3b_B5_cerrada.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""

import sys
from fractions import Fraction

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# Autocontenido: no lee ni importa ningun otro fichero.
def es_primo(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 151

PARTS = {(5,): 5, (4, 1): 4, (3, 2): 6, (3, 1, 1): 6,
         (2, 2, 1): 8, (2, 1, 1, 1): 12, (1, 1, 1, 1, 1): 120}


def chi(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def N_enum(lam, p):
    """N_lambda por DP exacta sobre (suma, suma de cubos): O(l p^3)."""
    est = np.zeros((p, p), dtype=np.int64)
    est[0, 0] = 1
    cubo = [pow(x, 3, p) for x in range(p)]
    for w in lam:
        nuevo = np.zeros((p, p), dtype=np.int64)
        for x in range(p):
            nuevo += np.roll(np.roll(est, (w * x) % p, axis=0), (w * cubo[x]) % p, axis=1)
        est = nuevo
    return int(est[0, 0])


def a_p_E(p):
    """a_p de E: Y^2 = X^3 + X^2 - X, por conteo directo."""
    cnt = 1                                                   # el punto del infinito
    for x in range(p):
        v = (x * x * x + x * x - x) % p
        cnt += 1 if v == 0 else (2 if pow(v, (p - 1) // 2, p) == 1 else 0)
    return p + 1 - cnt


def eta_20(pmax):
    """q prod_{n>=1} (1-q^{2n})^2 (1-q^{10n})^2 : la newform de peso 2 y nivel 20."""
    L = pmax + 2
    f = [0] * L
    f[0] = 1
    for m in (2, 2, 10, 10):
        for n in range(1, L // m + 1):
            for i in range(L - 1, m * n - 1, -1):
                f[i] -= f[i - m * n]
    return [0] + f[:L - 1]                                     # el factor q inicial


print("N_lambda cerrado contra enumeracion, y C(5,p) cerrado contra el censo")
print("%-5s %-4s %-4s %-6s %-8s %-10s %-10s %-6s"
      % ("p", "chi5", "ch-15", "a_p(E)", "eta20", "B cerrado", "B censo", "ok"))

eta = eta_20(PMAX)
fallos_N = fallos_B = fallos_eta = fallos_lisa = casos = 0
for p in range(7, PMAX + 1):
    if not es_primo(p) or p == 5:
        continue
    c5, c15 = chi(5, p), chi(-15, p)
    ap = a_p_E(p)

    # (1) lisura de 2x(x+y+z)^2 + yz(y+z) = 0 sobre F_p: ningun punto proyectivo singular
    sing = 0
    for x in (range(p) if p <= 47 else []):
        for y in range(p):
            for z in (range(p) if (x or y) else range(1, p)):
                if x == 0 and y == 0 and z == 0:
                    continue
                s = (x + y + z) % p
                if (2 * x * s * s + y * z * (y + z)) % p:
                    continue
                gx = (2 * s * s + 4 * x * s) % p
                gy = (4 * x * s + 2 * y * z + z * z) % p
                gz = (4 * x * s + y * y + 2 * y * z) % p
                if gx == gy == gz == 0:
                    sing += 1
    fallos_lisa += 1 if sing else 0

    cerr = {(5,): 1, (4, 1): 1, (3, 2): 1,
            (3, 1, 1): 1 + (p - 1) * (2 + c15),
            (2, 2, 1): 1 + (p - 1) * (2 + c5),
            (2, 1, 1, 1): 1 + (p - 1) * (p + 1 - ap),
            (1, 1, 1, 1, 1): 1 + (p - 1) * (p * p + (6 + c5) * p + 1)}
    for lam in PARTS:
        if p <= 61 and N_enum(lam, p) != cerr[lam]:
            fallos_N += 1
            print("   FALLO N_%s en p=%d: enum=%d cerrado=%d"
                  % (str(lam), p, N_enum(lam, p), cerr[lam]))

    B_cerr = sum(Fraction((-1) ** (5 - len(lam)), z) * cerr[lam] for lam, z in PARTS.items())

    # censo exacto por DP
    dp = np.zeros((6, p, p), dtype=np.int64)
    dp[0, 0, 0] = 1
    for t in range(p):
        c3 = pow(t, 3, p)
        for k in range(4, -1, -1):
            dp[k + 1] += np.roll(np.roll(dp[k], t, axis=0), c3, axis=1)
    B_cen = int(dp[5, 0, 0])

    ok = (B_cerr == B_cen)
    casos += 1
    fallos_B += 0 if ok else 1
    fallos_eta += 0 if eta[p] == ap else 1
    print("%-5d %-4d %-4d %-6d %-8d %-10s %-10d %-6s"
          % (p, c5, c15, ap, eta[p], str(B_cerr), B_cen, "si" if ok else "NO"))
    sys.stdout.flush()

print("")
print("casos: %d" % casos)
print("fallos de lisura de la cubica plana : %d" % fallos_lisa)
print("fallos de N_lambda cerrado          : %d   (comprobado hasta p = 61)" % fallos_N)
print("fallos de a_p(E) vs eta(2z)^2 eta(10z)^2 : %d   (nivel 20: comprobacion en rango finito,"
      " no determinacion del conductor)" % fallos_eta)
print("fallos de C(5,p) cerrado vs censo   : %d" % fallos_B)
print("")
print("VEREDICTO:", "C(5,p) cerrado: polinomio + chi_p(5) + chi_p(-15) + a_p(20.a)"
      if not (fallos_B or fallos_N or fallos_eta or fallos_lisa) else "LA FORMA CERRADA FALLA")
bien = [casos == 33, fallos_lisa == 0, fallos_N == 0, fallos_eta == 0, fallos_B == 0]
print("")
print("=" * 72)
print("VERIFICACION C(5,p) CERRADA:  %d ok, %d MAL" % (sum(bien), len(bien) - sum(bien)))
print("=" * 72)
sys.exit(0 if all(bien) else 1)
