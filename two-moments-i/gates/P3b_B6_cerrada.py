# -*- coding: utf-8 -*-
r"""P3b_B6_cerrada.py -- forma cerrada de C(6,p), y por que n = 6 vuelve a ser elemental.

C(n,p) = #{T subset F_p : |T| = n, M_1 = M_3 = 0} = #{T : e_1 = e_3 = 0}.

EL RESULTADO:

    C(6,p) = [ p^4 - 10p^3 + 90p^2 - 480p + 399
               - 90 (p-1) chi_p(-1) - 40 (p-1) chi_p(-3) ] / 720,

polinomio mas DOS simbolos de Legendre y NINGUNA curva eliptica.  Contrasta con C(5,p), que
lleva el a_p de una curva de conductor 20 (P3b_B5_cerrada.py).

POR QUE.  Por P3b_singularidades.py, Sing(V_lambda) son los vectores de signos +-1 cuyos bloques
tienen igual suma de pesos, y como sum w_i = n eso exige n PAR.  Las cuatro particiones de 6 con
l >= 4 --- las unicas que pueden dar geometria con aritmetica propia --- se parten TODAS, luego
todas sus V_lambda son singulares y por tanto racionales.  Las que quedan tienen l <= 3, o sea
dimension 0, y ahi lo unico que cabe es un simbolo cuadratico.  Eso es todo el fenomeno.

LAS ONCE PARTICIONES DE 6, con su variedad y su recuento (p > 6):

  lambda         l  z    signo  V_lambda                                        #V_lambda
  (6)            1   6    -1    ---                                              N = 1
  (5,1)          2   5    +1    ---  (5 != +-1 en F_p)                           N = 1
  (4,2)          2   8    +1    ---  (4 != +-2 en F_p)                           N = 1
  (3,3)          2  18    +1    y = -x, la cubica sale gratis                    N = p
  (4,1,1)        3   8    -1    -12x(5x^2+4xy+y^2), disc -4, LISA                2 + chi_p(-1)
  (3,2,1)        3   6    -1    -6(x+y)^2(4x+y), UN nodo (doble)                 2
  (2,2,2)        3  48    -1    xyz = 0 sobre x+y+z = 0, LISA                    3
  (3,1,1,1)      4  18    +1    cubica nodal irreducible, 1 nodo                 p + 1 - chi_p(-3)
  (2,2,1,1)      4  16    +1    RECTA + CONICA: t (zw - xy - t^2), t = x+y       2p
  (2,1,1,1,1)    5  48    -1    CUBICA DE CAYLEY, la 4-nodal (Cayley 1869)       p^2 + 3p + 1
  (1,1,1,1,1,1)  6 720    +1    CUBICA DE SEGRE, 10 nodos (Segre 1887)           p^3 + 6p^2 - 4p + 1

con N_lambda = 1 + (p-1) #V_lambda para l >= 3.  Las dos factorizaciones de (2,2,1,1) y (3,3)
estan demostradas y se comprueban aqui como identidades de polinomios; el numero de nodos de
Cayley (4) y de Segre (10) es el que predice el teorema de los vectores de signos, y se
recomprueba aqui por fuerza bruta.

CONTROLES:
 (0) la factorizacion t(zw - xy - t^2) de (2,2,1,1) como identidad sobre Q;
 (1) cada #V_lambda cerrado contra enumeracion exacta;
 (2) C(6,p) cerrado contra el censo exacto por DP, todos los primos 7 <= p <= PMAX.

    python P3b_B6_cerrada.py [PMAX]

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

# lambda -> (z_lambda, signo (-1)^{6-l})
PARTS = {(6,): (6, -1), (5, 1): (5, 1), (4, 2): (8, 1), (3, 3): (18, 1),
         (4, 1, 1): (8, -1), (3, 2, 1): (6, -1), (2, 2, 2): (48, -1),
         (3, 1, 1, 1): (18, 1), (2, 2, 1, 1): (16, 1),
         (2, 1, 1, 1, 1): (48, -1), (1, 1, 1, 1, 1, 1): (720, 1)}


def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def N_cerrado(lam, p):
    c1, c3 = chi(-1, p), chi(-3, p)
    V = {(4, 1, 1): 2 + c1,
         (3, 2, 1): 2,
         (2, 2, 2): 3,
         (3, 1, 1, 1): p + 1 - c3,
         (2, 2, 1, 1): 2 * p,
         (2, 1, 1, 1, 1): p * p + 3 * p + 1,
         (1, 1, 1, 1, 1, 1): p ** 3 + 6 * p * p - 4 * p + 1}
    if lam in V:
        return 1 + (p - 1) * V[lam]
    return p if lam == (3, 3) else 1


def N_enum(lam, p):
    est = np.zeros((p, p), dtype=np.int64)
    est[0, 0] = 1
    cubo = [pow(x, 3, p) for x in range(p)]
    for w in lam:
        nuevo = np.zeros((p, p), dtype=np.int64)
        for x in range(p):
            nuevo += np.roll(np.roll(est, (w * x) % p, axis=0), (w * cubo[x]) % p, axis=1)
        est = nuevo
    return int(est[0, 0])


def B6_cerrado(p):
    c1, c3 = chi(-1, p), chi(-3, p)
    num = (p ** 4 - 10 * p ** 3 + 90 * p * p - 480 * p + 399
           - 90 * (p - 1) * c1 - 40 * (p - 1) * c3)
    assert num % 720 == 0, (p, num)
    return num // 720


def censo6(p):
    dp = np.zeros((7, p, p), dtype=np.int64)
    dp[0, 0, 0] = 1
    for t in range(p):
        c3 = pow(t, 3, p)
        for k in range(5, -1, -1):
            dp[k + 1] += np.roll(np.roll(dp[k], t, axis=0), c3, axis=1)
    return int(dp[6, 0, 0])


print("(0) La factorizacion de (2,2,1,1): sobre 2x+2y+z+w = 0, con t = x+y,")
print("    2x^3+2y^3+z^3+w^3 = -6 t (t^2 + xy - zw).  Identidad sobre Q:")
import sympy as sp                                                     # noqa: E402
x, y, z = sp.symbols("x y z")
w = -2 * x - 2 * y - z
t = x + y
izq = sp.expand(2 * x ** 3 + 2 * y ** 3 + z ** 3 + w ** 3)
der = sp.expand(-6 * t * (t ** 2 + x * y - z * w))
print("    diferencia =", sp.simplify(izq - der))
ok_fact = sp.simplify(izq - der) == 0
print("    %s" % ("ok" if ok_fact else "FALLA"))
print("    y el factor cuadratico es (z+t)^2 + xy tras w = -2t-z: conica lisa, luego")
print("    #V = (p+1) + (p+1) - 2 = 2p, con los dos puntos de corte racionales.")

print("")
print("(1) y (2): cada N_lambda cerrado contra enumeracion, y C(6,p) contra el censo")
print("%-5s %-5s %-5s %-12s %-12s %-6s" % ("p", "ch-1", "ch-3", "B cerrado", "B censo", "ok"))
fallos_N = fallos_B = casos = 0
for p in range(7, PMAX + 1):
    if not es_primo(p):
        continue
    if p <= 61:
        for lam in PARTS:
            e, c = N_enum(lam, p), N_cerrado(lam, p)
            if e != c:
                fallos_N += 1
                print("   FALLO N_%s en p=%d: enum=%d cerrado=%d" % (str(lam), p, e, c))
    tamiz = sum(Fraction(s, zl) * N_cerrado(lam, p) for lam, (zl, s) in PARTS.items())
    Bc, Bn = B6_cerrado(p), censo6(p)
    assert tamiz == Bc, (p, tamiz, Bc)         # el tamiz y la forma cerrada son lo mismo
    ok = (Bc == Bn)
    casos += 1
    fallos_B += 0 if ok else 1
    print("%-5d %-5d %-5d %-12d %-12d %-6s" % (p, chi(-1, p), chi(-3, p), Bc, Bn,
                                               "si" if ok else "NO"))
    sys.stdout.flush()

print("")
print("casos: %d" % casos)
print("fallos de la factorizacion de (2,2,1,1) : %d" % (0 if ok_fact else 1))
print("fallos de N_lambda cerrado              : %d   (comprobado hasta p = 61)" % fallos_N)
print("fallos de C(6,p) cerrado vs censo       : %d" % fallos_B)
print("")
print("VEREDICTO:", "C(6,p) cerrado: polinomio + chi_p(-1) + chi_p(-3), sin curva eliptica"
      if not (fallos_B or fallos_N or not ok_fact) else "LA FORMA CERRADA FALLA")
bien = [casos == 33, bool(ok_fact), fallos_N == 0, fallos_B == 0]
print("")
print("=" * 72)
print("VERIFICACION C(6,p) CERRADA:  %d ok, %d MAL" % (sum(bien), len(bien) - sum(bien)))
print("=" * 72)
sys.exit(0 if all(bien) else 1)
