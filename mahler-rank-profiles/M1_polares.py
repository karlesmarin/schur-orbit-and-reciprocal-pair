# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: control de los pasos 1-2 de M1 (version coherente, PRUEBA_VUELTAS.md) en T[t]/(t^e), sin formula de producto:
#      c_k = e_k(z^{2i} : -m <= i <= m), z = alpha*E(t), calculado por el producto prod_i (1 + z^{2i} y).
#      Comprueba, para k <= m:
#       (a) n | k, k = n l: v(c_k) = P - p^{v_p(l)}, y c_k mod t^D tiene todos los coeficientes en F_p (sin alpha);
#       (b) n no divide k: (-1)^{k-1} c_k * [k] == [N] mod t^D  (es decir, g_k == S U(z^k) mod t^D, con [k] unidad);
#       (c) el termino S^2 U_k Y_h con k = 2 n p^{a-1}, h = n p^{a-1} tiene valoracion 2P - 2p^{a-1} < D (V2 falso),
#           y aun asi (a) se cumple para ese k;
#       (d) precision: e >= D + B, B = a(p-1)p^{a-1}/2.
#      D = 2P - p^{a-1}, N = Pn, [h] = z^h - z^{-h}.
# FECHA: 3-oct-2026
# Uso (docker, /work = ext_clasicos): sage -python M1_polares.py p n a v
import sys

from sage.all import GF, PolynomialRing, cyclotomic_polynomial, factorial, PowerSeriesRing, ZZ

p, n, a, v = (int(x) for x in sys.argv[1:5])
P = p ** a
N = P * n
m = (N - 1) // 2
e = (p - 1) * p ** (v - 1)
D = 2 * P - p ** (a - 1)
B = a * (p - 1) * p ** (a - 1) // 2
print("(p,n,a,v) = (%d,%d,%d,%d): P=%d, N=%d, m=%d, e=%d, D=%d, B=%d, e >= D+B: %s"
      % (p, n, a, v, P, N, m, e, D, B, e >= D + B), flush=True)
assert e >= D + B

F = GF(p)
Rx = PolynomialRing(F, 'x')
phi = Rx(cyclotomic_polynomial(n))
f0 = sorted(phi.factor(), key=lambda fa: [int(c) for c in fa[0].list()])[0][0]  # un factor (convencion lexicografica)
K = F.extension(f0, 'al')
al = K.gen()
S_ = PowerSeriesRing(K, 't', default_prec=e)
t = S_.gen()
h = sum(t ** j / F(factorial(j)) for j in range(p))
E = (h / h(-t)).sqrt()
assert E[0] == 1
z = al * E
zi = z ** -1


def val(s):
    s = s.add_bigoh(e)
    return s.valuation() if s != 0 else e


def br(k):
    return z ** k - zi ** k


# e_k(z^{2i}), -m <= i <= m: polinomio en y truncado en grado m
coef = [S_(1)] + [S_(0)] * m
z2 = z ** 2
w = zi ** (2 * m)
for i in range(-m, m + 1):
    for k in range(m, 0, -1):
        coef[k] = coef[k] + w * coef[k - 1]
    w = w * z2
print("c_k calculados (k <= %d)" % m, flush=True)

en_Fp = lambda s: all(K(c) in F and K(c).polynomial().degree() <= 0 for c in s.add_bigoh(D).list())  # noqa: E731
BN = br(N)
fa = 0
na = nb = 0
for k in range(1, m + 1):
    ck = coef[k]
    if k % n == 0:
        l = k // n
        esperado = P - p ** ZZ(l).valuation(p)
        ok = val(ck) == esperado and en_Fp(ck)
        na += 1
    else:
        ok = val(((-1) ** (k - 1) * ck * br(k) - BN).add_bigoh(D)) >= D
        nb += 1
    if not ok:
        fa += 1
        print("  FALLO k=%d" % k, flush=True)
print("(a) n|k: %d casos; (b) n no divide k: %d casos; fallos: %d" % (na, nb, fa), flush=True)

k2, h2 = 2 * n * p ** (a - 1), n * p ** (a - 1)
if k2 <= m:
    SS = BN / 2
    term = SS ** 2 * (2 / br(k2)) * ((z ** h2 + zi ** h2) / br(h2))
    print("(c) k=%d, h=%d: v(S^2 U_k Y_h) = %d  (2P-2p^{a-1} = %d, D = %d); c_k sin alpha mod t^D: %s"
          % (k2, h2, val(term), 2 * P - 2 * p ** (a - 1), D, en_Fp(coef[k2])), flush=True)
print("RESULTADO:", "OK" if fa == 0 else "FALLA", flush=True)
