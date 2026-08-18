# -*- coding: utf-8 -*-
# LAS DOS AFIRMACIONES DE SU VUELTA 24, VERIFICADAS ANTES DE TOCAR NADA.   16 de agosto de 2026.
#
# A) SU ERROR MATEMATICO 1: el filtro impar del Lema 3.1 (simplectico) NO es el de la Prop 5.3
#    (ortogonal impar).  Su testigo:  t=5, m'=2, eta=(1,0):
#         sp_{(1,0)}(xi, xi^2)              = xi + xi^-1 + xi^2 + xi^-2   = -1
#         o^{B_2}_{(1,0)}(1, xi^{+-1}, xi^{+-2}) = 1 + eso                = 0
#    Si es asi, el Lema 3.1 NO puede seguir diciendo "and for odd t the same with m' in place of m",
#    y hay que separarlo por paridad.  Se comprueba en TODA la caja, no solo en su eta.
#
# B) SU CONTRAEJEMPLO A "bajo anulacion, e <= 4":
#         t=6, r=3, beta = (12,11,10,9,8,7,5,4,3,2,1,0),  lambda = (1^6),
#    con e = 6 y |g_com| = 4, y Phi = 0.  El da ademas una prueba de una linea:
#         s_{(1^6)} = e_6,  prod_{alpha en mu_6}(1 + alpha u) = 1 - u^6,
#         y los tres pares libres aportan coeficiente lider 1 en u^6, luego [u^6] = 1 - 1 = 0.
#    Se comprueba por DOS rutas: el bialternante, y el coeficiente de u^6 del producto.
#
# CONTROLES
#   C0  en (A), las dos rutas se calculan con Freudenthal, que no sabe de nuestra discusion.
#   C1  en (B), el bialternante y la cuenta de e_6 no comparten una linea de codigo.
#   C2  se imprime el numero de eta donde los dos filtros DIFIEREN, no solo si difieren.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage check_24in.sage

import json
from collections import Counter

print("=" * 112)
print("A)  ¿SON DISTINTOS EL FILTRO SIMPLECTICO Y EL ORTOGONAL IMPAR?")
print("=" * 112)

_CH = {}
def car(typ, rk, mu):
    key = (typ, rk, tuple(int(v) for v in mu))
    if key not in _CH:
        W = WeylCharacterRing("%s%d" % (typ, rk))
        el = W(W.space().from_vector(vector([Integer(v) for v in mu])))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + int(mult)
        _CH[key] = d
    return _CH[key]


def evalua(typ, rk, eta, t):
    """el caracter de tipo typ evaluado en el toro x_i = xi^i."""
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car(typ, rk, eta).items():
        s += mult * z ** (sum((i + 1) * wt[i] for i in range(rk)) % t)
    return QQ(s) if s in QQ else None


A = []
for t in (3, 5, 7, 9):
    mp = (t - 1) // 2
    n = dif = 0
    primero = None
    for k in range(0, 3 * t + 1):
        for e in Partitions(k, max_length=mp):
            eta = tuple(list(e) + [0] * (mp - len(e)))
            vC = evalua("C", mp, eta, t)      # sp_eta(xi, ..., xi^{m'})  -- el Lema 3.1
            vB = evalua("B", mp, eta, t)      # o_eta(1, xi^{+-1..m'})    -- la Prop 5.3
            if vC is None or vB is None:
                continue
            n += 1
            if vC != vB:
                dif += 1
                if primero is None:
                    primero = (eta, vC, vB)
    print("  t=%2d  m'=%d :  eta probados %4d  |  los dos filtros DIFIEREN en %4d  (%.1f %%)"
          % (t, mp, n, dif, 100.0 * dif / n))
    if primero:
        print("        primer testigo: eta=%s   sp=%s   o=%s" % (str(primero[0]), primero[1], primero[2]))
    A.append({"t": int(t), "mp": int(mp), "n": int(n), "difieren": int(dif),
              "testigo": [list(map(int, primero[0])), int(primero[1]), int(primero[2])] if primero else None})

# su eta exacto
print("")
print("  SU TESTIGO, t=5, eta=(1,0):  sp = %s   o = %s"
      % (evalua("C", 2, (1, 0), 5), evalua("B", 2, (1, 0), 5)))
A_j = [2 * 1 + 2 * (2 - 0 - 1) + 1, 2 * 0 + 2 * (2 - 1 - 1) + 1]
print("        y su A_j = 2 eta_j + 2(m'-j) + 1 = %s,  con A_1 = %d = 0 (mod 5) -> el peso muere"
      % (str(A_j), A_j[0]))

print("")
print("=" * 112)
print("B)  SU CONTRAEJEMPLO A  'bajo anulacion, e <= 4'")
print("=" * 112)

t, r = 6, 3
N = t + 2 * r
beta = (12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 0)
delta = list(range(N - 1, -1, -1))
lam = tuple(beta[i] - delta[i] for i in range(N))
print("  beta   = %s" % str(beta))
print("  lambda = %s   (|lambda| = %d)" % (str(tuple(v for v in lam if v > 0)), sum(lam)))
cl = Counter(v % t for v in beta)
E = sorted(i for i in cl if cl[i] >= 2)
print("  clases mod %d: %s   ->  clases con >= 2 (el exceso): %s,  e = %d"
      % (t, dict(sorted(cl.items())), str(E), len(E)))

# RUTA 1: el bialternante
K = CyclotomicField(t)
zeta = K.gen()
L = LaurentPolynomialRing(K, r, 'z')
zs = L.gens()
x = [L(K(zeta) ** k) for k in range(t)] + [g ** e for g in zs for e in (1, -1)]
den = matrix(L, N, N, lambda i, j: x[i] ** delta[j]).determinant()
num = matrix(L, N, N, lambda i, j: x[i] ** beta[j]).determinant()
print("  RUTA 1  bialternante:  numerador == 0 ?  %s" % (num == 0))

# RUTA 2: su cuenta de una linea.  s_{(1^6)} = e_6, y e_6 del alfabeto entero
#         = coeficiente de u^6 en prod_{x en alfabeto} (1 + x u)
P = PolynomialRing(K, 'u')
u = P.gen()
Lu = LaurentPolynomialRing(K, r, 'z')
Pu = PolynomialRing(Lu, 'u')
uu = Pu.gen()
prod = Pu(1)
for k in range(t):
    prod *= (1 + K(zeta) ** k * uu)
for g in Lu.gens():
    prod *= (1 + g * uu) * (1 + g ** (-1) * uu)
c6 = prod[6]
print("  RUTA 2  el coeficiente de u^6 en prod (1 + x u):  %s" % ("0" if c6 == 0 else str(c6)))
print("          (y prod sobre la orbita = 1 - u^%d, luego su termino u^6 es -1;" % t)
print("           los tres pares libres aportan +1 en u^6: 1 - 1 = 0)")

json.dump({"A_filtros_distintos": A,
           "B_contraejemplo": {"beta": list(beta), "e": int(len(E)),
                               "numerador_cero": bool(num == 0),
                               "coef_u6_cero": bool(c6 == 0)}},
          open("check_24in_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("DONE")
