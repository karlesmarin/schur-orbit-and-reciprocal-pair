# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: control de M1 carácter a carácter (no por dimensiones totales). Para (n, p, a):
#      (1) capas G_i del orden real (conjetura_C.Caso: cierre de A en F_p[x]/psi_n^L, x = zeta_Q + zeta_Q^{-1});
#      (2) cada vector de G_i (coordenadas en la base de potencias de psi) pasa a la coordenada t:
#          coeficiente de t^i = r(alpha + alpha^{-1}) * w^i, w = psi_n'(alpha + alpha^{-1}) (alpha - alpha^{-1});
#      (3) componente chi: pi_chi(y) = sum_c conj(chi)(c) y(omega^c), con omega raiz n-esima primitiva fija en K;
#          dim_chi G_i = rango de {pi_chi(y)} (0 o 1). Estabilidad de Galois: sum_chi dim_chi G_i == W_i;
#      (4) filas de Mahler con la MISMA omega y los MISMOS valores de chi:
#          A_{j,i} = pi_chi [t^j] ((E_w - 1)^i U)(alpha), U(Z) = 2Z/(Z^2 - 1), E_w: Z -> Z E(t);
#      (5) M1: para chi impar y 0 <= j <= P - p^{a-1} - 1, dim_chi G_{P+j} = 1  sii  fila j nueva.
#      Sin etiquetas de Conrey ni convención de primo: los dos lados usan la misma reducción.
# FECHA: 3-oct-2026
# Uso (docker, /work = ext_clasicos): sage -python M1_por_caracter.py n p a
import sys
import time

sys.path.insert(0, '/work')
from conjetura_C import Caso  # noqa: E402
from nucleo import psi  # noqa: E402
from sage.all import GF, DirichletGroup, PowerSeriesRing, PolynomialRing, Integers, carmichael_lambda, factorial, lcm, matrix  # noqa: E402

n, p, a = (int(x) for x in sys.argv[1:4])
P = p ** a
jmax = P - p ** (a - 1) - 1
L = P + jmax + 1
t0 = time.time()
caso = Caso(n, p, L, a)
B = caso.cierre()
G = caso.capas(B)
W = [len(g) for g in G]
d = caso.d
print("(n,p,a) = (%d,%d,%d): P=%d, j <= %d, capas hasta %d, d=%d, cierre %.1f s"
      % (n, p, a, P, jmax, L - 1, d, time.time() - t0), flush=True)

# cuerpo con raices n-esimas y valores de caracteres
ex = lcm(n, carmichael_lambda(n))
f = 1
while (p ** f - 1) % ex:
    f += 1
K = GF(p ** f, 'g')
om = K.multiplicative_generator() ** ((p ** f - 1) // n)
assert om ** n == 1 and all(om ** k != 1 for k in range(1, n))
unidades = [c for c in range(1, n) if Integers(n)(c).is_unit()]
chars = list(DirichletGroup(n, K))
print("K = GF(%d^%d); %d caracteres" % (p, f, len(chars)), flush=True)

Fp = PolynomialRing(GF(p), 'x')
ps = Fp([int(c) for c in psi(n)])
dps = ps.derivative()


def ev(poly, val):
    return sum(K(int(c)) * val ** k for k, c in enumerate(poly.list()))


# valores en cada rama c: x0 = omega^c + omega^-c, w = psi'(x0)(omega^c - omega^-c)
x0 = {c: om ** c + om ** (-c) for c in unidades}
w = {c: ev(dps, x0[c]) * (om ** c - om ** (-c)) for c in unidades}
assert all(ev(ps, x0[c]) == 0 for c in unidades)


def pi(chi, vals):
    return sum(chi(c) ** -1 * vals[c] for c in unidades)


# (3) componentes medidas
med = {}
for i in range(P, L):
    for k, chi in enumerate(chars):
        vs = [pi(chi, {c: ev(Fp(v), x0[c]) * w[c] ** i for c in unidades}) for v in G[i]]
        med[(i, k)] = 1 if any(x != 0 for x in vs) else 0
estable = all(sum(med[(i, k)] for k in range(len(chars))) == W[i] for i in range(P, L))
print("(3) sum_chi dim_chi G_i == W_i en i = %d..%d: %s" % (P, L - 1, estable), flush=True)

# (4) filas de Mahler
S = PowerSeriesRing(K, 't', default_prec=jmax + 1)
t = S.gen()
h = sum(t ** j / K(factorial(j)) for j in range(p))
E = (h / h(-t)).sqrt()
Ur = {}
for c in unidades:
    Z = om ** c
    Ur[c] = [2 * Z * E ** r / (Z ** 2 * E ** (2 * r) - 1) for r in range(P)]
dif = {}
for c in unidades:
    col = []
    cur = Ur[c]
    for i in range(P):
        col.append(cur[0])
        cur = [cur[r + 1] - cur[r] for r in range(len(cur) - 1)]
    dif[c] = col  # dif[c][i] = ((E_w - 1)^i U)(omega^c)

fallos = 0
tot = 0
for k, chi in enumerate(chars):
    if chi(-1) != -1 or not chi.is_primitive():
        continue
    A = matrix(K, jmax + 1, P, lambda j, i: pi(chi, {c: dif[c][i][j] for c in unidades}))
    rk_prev = 0
    linea = []
    for j in range(jmax + 1):
        rk = A.matrix_from_rows(range(j + 1)).rank()
        nueva = 1 if rk > rk_prev else 0
        rk_prev = rk
        m_ = med[(P + j, k)]
        tot += 1 - (j % 2)
        if nueva != m_ or (j % 2 and nueva):
            fallos += 1
        linea.append('.' if nueva == m_ == 1 else ('x' if nueva == m_ == 0 else '!'))
    orden = chi.multiplicative_order()
    aus = [j for j in range(0, jmax + 1, 2) if med[(P + j, k)] == 0]
    print("chi #%d (orden %d): filas pares ausentes %s" % (k, orden, aus), flush=True)
    print("   %s" % ''.join(linea), flush=True)
print("(5) chi impares primitivos: %d comparaciones con j par (las de j impar son nulas en los dos lados), %d discrepancias; %.1f s"
      % (tot, fallos, time.time() - t0), flush=True)
print("RESULTADO:", "OK" if fallos == 0 and estable else "FALLA", flush=True)
