# -*- coding: utf-8 -*-
# .ES NUESTRO nu UN MULTIPLETE GKRS EN EL SENTIDO DE LANDWEBER-SJAMAAR?
#
# Verificado en el PDF de arXiv:1107.3578 (Selecta Math. 19 (2013) 49-95), textual:
#
#   J_M^op = sum_{w in W^H} det(w) w^{-1}          (p. 28)
#   J_G = J_H J_M^op                               (4.3)
#   j_H^*(e(D_M)^*) . d_G  =  d_H . J_M^op         (4.4)
#   Teorema 4.2.2:  sum_{w in W^H} det(w) f^*(a_w) = 0,   con  a_w = d_H(w^{-1}(a))
#   y la prueba:    f^* j_H^*(e(D_M)^*) = prod_{alpha in R_M^+} (1 - 1) = 0
#
# f^* es el mapa de OLVIDO de la equivariancia: e^alpha -> 1.  En nuestro lenguaje eso es evaluar el
# caracter en la identidad, o sea tomar dimensiones.  Asi que su teorema, traido a nuestras
# coordenadas, predice
#
#       sum_mu  nu(Lambda, mu) . dim chi^{D_r}_mu  =  0.
#
# Si sale, nuestro nu ES un multiplete GKRS en su sentido y su marco es el nuestro.  Y entonces queda
# EXACTAMENTE localizada la distancia entre su teorema y nuestra conjetura:
#
#   ellos    la suma alternante TOTAL se anula tras olvidar la equivariancia
#   nosotros cada FIBRA de la extraccion de coeficientes se anula por separado
#
# Lo segundo es estrictamente mas fuerte y no se sigue de lo primero.  Este gate mide lo primero,
# para poder decir en el paper donde acaba lo publicado y donde empieza lo nuestro.
#
# CONTROLES
#   M0 (FATAL)  sum nu . dim  ==  0, peso a peso.
#   M1          .y es no trivial?  Se imprime cuantos sumandos tiene y la mayor dimension, para que
#               "0 = 0" no sea 0 = suma vacia.
#   M2          la fibra: sum sobre UNA progresion de nu . dim, .tambien es 0?  NO tiene por que:
#               eso seria mucho mas fuerte.  Se mide para ver si por casualidad tambien pasa.
#   N1 senuelo  los mismos soportes con los signos de nu barajados: M0 debe hundirse.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_multiplete.py

import itertools
import json
import random
from collections import Counter

from divided_differences import CASOS, nu_de, nu_extendida, enderezar_D


def dim_D(X, r):
    """dim chi^{D_r}_mu con X = 2(mu + rho_{D_r}).  Formula de Weyl para D_r."""
    X0 = [2 * (r - 1 - i) for i in range(r)]
    num = den = 1
    for i in range(r):
        for j in range(i + 1, r):
            num *= (X[i] ** 2 - X[j] ** 2)
            den *= (X0[i] ** 2 - X0[j] ** 2)
    if den == 0:
        return None
    if num % den != 0:
        return None
    return num // den


print("=" * 100)
print(".ES nu UN MULTIPLETE GKRS?   suma alternante tras olvidar la equivariancia")
print("=" * 100)
print("")

M0 = M0n = 0
N1 = N1n = 0
M2 = M2n = 0
malos = []
tam = Counter()
rnd = random.Random(20260816)

for (t, r, cota) in CASOS:
    Rp = (t - 1) // 2 + r
    m0 = m0n = 0
    maxdim = 0
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu = nu_de(list(Lam), t, r)
        if not nu:
            continue
        dims = {}
        malo = False
        for Y, v in nu.items():
            d = dim_D(list(Y), r)
            if d is None:
                malo = True
                break
            dims[Y] = d
        if malo:
            continue
        s = sum(v * dims[Y] for Y, v in nu.items())
        m0n += 1
        M0n += 1
        tam[len(nu)] += 1
        maxdim = max([maxdim] + list(dims.values()))
        if s == 0:
            m0 += 1
            M0 += 1
        elif len(malos) < 3:
            malos.append({"t": t, "r": r, "Lambda": list(Lam),
                          "nu": {str(k): v for k, v in nu.items()},
                          "dims": {str(k): v for k, v in dims.items()}, "suma": s})
        # N1  senuelo: los mismos soportes, signos barajados
        falso = {k: rnd.choice([1, -1]) for k in nu}
        N1n += 1
        if sum(v * dims[Y] for Y, v in falso.items()) == 0:
            N1 += 1
        # M2  .se anula tambien FIBRA a FIBRA con las dimensiones dentro?
        M = max(max(abs(v) for v in k) for k in nu)
        for X in itertools.product(range(-M, M + 1), repeat=r):
            if any(X[j] <= X[j + 1] for j in range(r - 2)):
                continue
            if not (X[r - 2] > abs(X[r - 1])):
                continue
            ks = []
            for j in range(r):
                L, k = [], 1
                while X[j] + t * k <= M:
                    if X[j] + t * k >= -M:
                        L.append(k)
                    k += 2
                ks.append(L)
            if any(not L for L in ks):
                continue
            tot = 0
            nz = 0
            for k in itertools.product(*ks):
                Y = tuple(X[j] + t * k[j] for j in range(r))
                v = nu_extendida(nu, Y)
                if v:
                    nz += 1
                    e = enderezar_D(Y)
                    d = dim_D(list(e[0]), r)
                    if d is None:
                        nz = -1
                        break
                    tot += v * d
            if nz >= 2:
                M2n += 1
                if tot == 0:
                    M2 += 1
    print("  t=%d r=%d :  sum nu.dim == 0 en %4d de %4d   (mayor dimension vista: %d)"
          % (t, r, m0, m0n, maxdim))

print("")
print("-" * 100)
print("  M0  FATAL  sum_mu nu . dim = 0                 : %d de %d" % (M0, M0n))
print("  M1  tamano del soporte de nu (para que 0 no sea la suma vacia): %s"
      % dict(sorted(tam.items())))
print("  M2  .tambien fibra a fibra, con dimensiones?   : %d de %d" % (M2, M2n))
print("      (no tiene por que: seria mucho mas fuerte que su Teorema 4.2.2)")
print("")
print("  N1  SENUELO signos barajados                   : %d de %d" % (N1, N1n))
if malos:
    print("")
    print("  !! pesos donde NO se anula:")
    for m in malos:
        print("    " + json.dumps(m)[:400])
print("")
print("  LECTURA: si M0 sale entero y N1 se hunde, nuestro nu es un multiplete GKRS en el sentido de")
print("  LS 4.2.2, su marco es el nuestro, y la distancia exacta entre su teorema y nuestra conjetura")
print("  es la que separa 'la suma total se anula' de 'cada fibra se anula'.")

json.dump({"M0": [M0, M0n], "M2": [M2, M2n], "N1_senuelo": [N1, N1n],
           "tam_soporte": {str(k): v for k, v in sorted(tam.items())}, "malos": malos},
          open("_probe_multiplete_DUMP.json", "w"), indent=1)
print("")
print("=" * 100)
print("DONE")
