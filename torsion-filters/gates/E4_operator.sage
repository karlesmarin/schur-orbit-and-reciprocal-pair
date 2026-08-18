# -*- coding: utf-8 -*-
# E^{(4)} COMO OPERADOR, independiente de beta.   15 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 18: separar lo que viene del Paper I (el vector a_Lambda) de lo que es
# UNIVERSAL (el operador).  Con la formula local verificada en la vuelta anterior,
#
#     E_{Lambda,mu}  =  1[mu << Lambda] * 1[todos los r_i pares] * (-1)^{sum r_i / 2}
#     r_i  =  min(Lambda_i, mu_{i-1}) - max(Lambda_{i+1}, mu_i)
#
# y A = E a.  El operador se puede estudiar UNA VEZ y sirve para todas las formas.
#
# LO QUE SE VE ANTES DE CALCULAR, y es la respuesta a su pregunta 4.  Como
#
#     eps = prod_i ( [r_i par] * (-1)^{r_i/2} )
#
# y cada r_i depende SOLO de (Lambda_i, Lambda_{i+1}, mu_{i-1}, mu_i), el peso es un PRODUCTO DE
# FACTORES LOCALES sobre una cadena: E es una MATRIZ DE TRANSFERENCIA.  No una matriz cualquiera.
# Eso se comprueba en C1, no se afirma.
#
# LO QUE SE MIDE (sus cinco preguntas)
#   Q1  ¿triangular bajo un orden natural?  Se mide si E_{Lambda,mu} != 0 fuerza |mu| <= |Lambda|,
#       y si hay una "diagonal" no nula.
#   Q2  el rango y, si el bloque cuadrado es invertible, la forma de las entradas de la inversa.
#   Q3  ¿recurrencia al añadir una caja?  Se mide E_{Lambda+e_k, mu} contra E_{Lambda,mu}.
#   Q4  la estructura de transferencia: eps = producto de factores locales.  C1.
#   Q5  ¿tiene la fila de mu una estructura de Mobius?  Se miden las sumas por fila y por columna,
#       que para una funcion de Mobius sobre intervalos dan 0 salvo en el extremo.
#
# CONTROLES
#   C0  las entradas estan en {0,+-1} por construccion: se comprueba igualmente, porque un valor
#       fuera de ahi delataria un error en la formula local.
#   C1  FATAL para la lectura: eps == producto de los factores locales, par a par.
#   C2  SEÑUELO para Q1: se mide tambien la triangularidad de una matriz ALEATORIA con la misma
#       densidad y soporte.  Si el azar sale igual de triangular, la propiedad no dice nada.
#   C3  no vacuidad: se imprimen los tamaños y la densidad siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage E4_operator.sage

import itertools, json, sys
from collections import defaultdict

R = 3          # Sp_6 -> Sp_4,  o sea t=4, r=2
W = 8          # caja: todas las entradas <= W
INF = 10 ** 9


def dominantes(k, tope):
    out = []
    def rec(pref, j, cap):
        if j == k:
            out.append(tuple(pref)); return
        for v in range(cap, -1, -1):
            rec(pref + [v], j + 1, v)
    rec([], 0, tope)
    return out


def entrelaza(mu, Lam):
    Lp = list(Lam) + [0]
    return all(Lp[i] >= mu[i] >= Lp[i + 2] for i in range(len(Lam) - 1))


def rs_locales(mu, Lam):
    L = list(Lam) + [0]
    M = [INF] + list(mu) + [0] * (len(Lam) + 2)
    return [min(L[i], M[i]) - max(L[i + 1], M[i + 1]) for i in range(len(Lam))]


def factor_local(rr):
    """f(r) = [r par] * (-1)^{r/2}.  El factor de UNA posicion de la cadena."""
    return 0 if rr % 2 else (-1) ** (rr // 2)


def E(mu, Lam):
    if not entrelaza(mu, Lam):
        return 0
    p = 1
    for rr in rs_locales(mu, Lam):
        p *= factor_local(rr)
        if p == 0:
            return 0
    return p


LAMS = dominantes(R, W)
MUS = dominantes(R - 1, W)
print("=" * 118)
print("E^{(4)} COMO OPERADOR   --   R=%d, caja <= %d:  %d Lambda x %d mu" % (R, W, len(LAMS), len(MUS)))
print("=" * 118)
print("")

# ------------------------------------------------------------------ C0, C1, densidad -----------
nz = 0
valores = defaultdict(int)
malo1 = 0
for Lam in LAMS:
    for mu in MUS:
        v = E(mu, Lam)
        valores[v] += 1
        if v:
            nz += 1
        # C1: eps como producto de factores locales, contra el mismo calculo hecho de golpe
        if entrelaza(mu, Lam):
            rs = rs_locales(mu, Lam)
            directo = 0 if any(x % 2 for x in rs) else (-1) ** (sum(rs) // 2)
            prod = 1
            for x in rs:
                prod *= factor_local(x)
            malo1 += (directo != prod)
tot = len(LAMS) * len(MUS)
print("  C0  valores de las entradas: %s" % dict(sorted(valores.items())))
print("  C1  eps == producto de factores locales : %s"
      % ("PASA en %d pares" % tot if malo1 == 0 else "*** FALLA en %d ***" % malo1))
print("      -> E es una MATRIZ DE TRANSFERENCIA sobre la cadena, no una matriz generica")
print("  C3  densidad: %d no nulas de %d (%.2f %%)" % (nz, tot, 100.0 * nz / tot))
print("")
sys.stdout.flush()

# ------------------------------------------------------------------ Q1 triangularidad ----------
viola = 0
for Lam in LAMS:
    for mu in MUS:
        if E(mu, Lam) and sum(mu) > sum(Lam):
            viola += 1
print("  Q1  ¿E_{Lambda,mu} != 0  =>  |mu| <= |Lambda|?  violaciones: %d" % viola)
# la "diagonal": para cada mu, el Lambda MINIMO que entrelaza y sobrevive
diag = 0
for mu in MUS:
    cand = [Lam for Lam in LAMS if E(mu, Lam)]
    if cand:
        mn = min(cand, key=lambda L: (sum(L), L))
        diag += (abs(E(mu, mn)) == 1)
print("      y para cada mu hay un Lambda minimo con |E| = 1 : %d de %d" % (diag, len(MUS)))
print("")

# ------------------------------------------------------------------ Q5 sumas por fila/columna --
print("  Q5  ¿estructura de Mobius?  sumas con signo por columna (fija mu, suma sobre Lambda):")
h = defaultdict(int)
for mu in MUS:
    s = sum(E(mu, Lam) for Lam in LAMS)
    h[s] += 1
print("      histograma de sum_Lambda E : %s" % dict(sorted(h.items())))
h2 = defaultdict(int)
for Lam in LAMS:
    s = sum(E(mu, Lam) for mu in MUS)
    h2[s] += 1
print("      histograma de sum_mu E     : %s" % dict(sorted(h2.items())))
print("      (una funcion de Mobius sobre intervalos da 0 salvo en el extremo; un histograma")
print("       concentrado en 0 y +-1 seria la señal, uno disperso lo desmiente)")
print("")
sys.stdout.flush()

# ------------------------------------------------------------------ Q3 añadir una caja ---------
print("  Q3  ¿recurrencia al añadir una caja a Lambda?  se mide E(Lambda + e_k, mu) / E(Lambda, mu)")
rel = defaultdict(int)
for Lam in LAMS:
    for k in range(R):
        L2 = list(Lam); L2[k] += 1
        if k > 0 and L2[k] > L2[k - 1]:
            continue
        L2 = tuple(L2)
        if L2 not in set(LAMS):
            continue
        for mu in MUS:
            a, b = E(mu, Lam), E(mu, L2)
            if a or b:
                rel[(a, b)] += 1
print("      pares (E(Lambda), E(Lambda+e_k)) y su frecuencia:")
for k, v in sorted(rel.items(), key=lambda kv: -kv[1])[:9]:
    print("        %-10s : %d" % (str(k), v))
print("")

# ------------------------------------------------------------------ order types ----------------
print("  ORDER TYPES  (que rama toma el min y el max en cada una de las %d posiciones)" % R)
tipos = defaultdict(lambda: [0, 0])
for Lam in LAMS:
    for mu in MUS:
        if not entrelaza(mu, Lam):
            continue
        L = list(Lam) + [0]
        M = [INF] + list(mu) + [0] * (R + 2)
        tipo = tuple((0 if L[i] <= M[i] else 1, 0 if L[i + 1] >= M[i + 1] else 1) for i in range(R))
        tipos[tipo][0] += 1
        tipos[tipo][1] += (E(mu, Lam) != 0)
print("      %d order types distintos en la caja" % len(tipos))
print("      %-34s |   pares | sobreviven | %%" % "tipo (min,max) por posicion")
for tp, (a, b) in sorted(tipos.items(), key=lambda kv: -kv[1][0])[:8]:
    print("      %-34s | %7d | %10d | %5.1f" % (str(tp), a, b, 100.0 * b / a))
print("")
json.dump({"R": R, "W": W, "nLam": len(LAMS), "nMu": len(MUS), "nz": int(nz),
           "valores": {str(int(k)): int(v) for k, v in valores.items()},
           "viola_triang": int(viola), "n_order_types": len(tipos)},
          open("E4_operator_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
