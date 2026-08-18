# -*- coding: utf-8 -*-
# .ES EL PLEGADO AFIN?  DESCOMPOSICION DEL SIGNO EN LOS PARES QUE SE CANCELAN.
#
# El plegado afin se MIDE antes de preguntarle.  Y ademas el sitio donde hay que
# medirlo NO es el bloque libre.  La antiinvariancia afin del bloque libre ya murio:
# nu~(X + 2t e_j) = -nu~(X) falla fuera del soporte (t=3, Lambda=(2,1,0): nu~(9,5)=+1 y nu~(9,11)=0).
#
# El plegado afin vive en el bloque CONGELADO, y ahi ya esta escrito en el paper: delta(A) es el
# determinante del elemento de W(B_{m'}) |x t Z^{m'} que devuelve A al alcove, y vale 0 en la pared.
# Es decir, tau^B_t ES el plegado afin de nivel t.  Asi que la pregunta afilada es:
#
#     en cada par de terminos que se cancelan, .QUE factor del signo se voltea?
#
# El signo de un termino de la progresion se factoriza en tres piezas independientes:
#
#     contribucion  =  s  .  sg  .  E  .  delta
#
#       s      signo del enderezado por W(D_r) al llevar Y = X + t k a su forma canonica
#       sg     signo de la permutacion de barajado (transversal delante, complemento detras)
#       E      epsilon_t, constante para cada t -- se cancela sola, no puede ser el mecanismo
#       delta  delta(A), el signo del PLEGADO AFIN del bloque congelado   <== la hipotesis
#
# LO QUE SE FALSA
#   A1  en los pares que se cancelan, delta se voltea y los otros dos factores no.  (afin puro)
#   A2  reparto conjunto de (s'/s, sg'/sg, delta'/delta) -- si sale una sola combinacion, ESA es la
#       involucion, la nombre como la nombre.
#   A3  .cambia el transversal?  .cambia la quiralidad?  Reparto de (S == S', qui == qui').
#   A4  los bloques congelados A y A' de un par, .tienen las mismas clases plegadas?  (Trivial si
#       ambos son transversales validos -- se mide para PROBAR que es trivial y no confundirlo con
#       contenido.)
#   A5  y el reciproco, que es lo que haria de esto un teorema: dos transversales cuyos bloques
#       libres caen en la misma progresion, .SIEMPRE tienen delta opuesto?
#
# SENUELO
#   G1  el mismo reparto sobre pares de terminos de progresiones DISTINTAS (elegidos al azar entre
#       los que si contribuyen).  Si A1 saliera igual ahi, no estaria midiendo la progresion.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_afin.py

import itertools
import json
import random
from collections import Counter

from divided_differences import (CASOS, plegar, sgn_perm, eps_t, delta_dec, enderezar_D)


def nu_con_factores(Lam, t, r):
    """nu, y para cada punto canonico del soporte los TRES factores de su signo."""
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = {}
    for i, x in enumerate(V):
        d.setdefault(plegar(x, t)[0], []).append(i)
    if any(not d.get(j) for j in range(1, mp + 1)):
        return {}, {}, V
    E = eps_t(t, mp)
    val, fac = {}, {}
    for pick in itertools.product(*[d[j] for j in range(1, mp + 1)]):
        S = frozenset(pick)
        if len(S) != mp:
            continue
        Sc = [i for i in range(Rp) if i not in S]
        A = sorted([V[i] for i in S], reverse=True)
        dv = delta_dec(A, t, mp)
        if not dv:
            continue
        for qui in (1, -1):
            libre = sorted([V[i] for i in Sc], reverse=True)
            libre[-1] *= qui
            orden = sorted(S, key=lambda i: -V[i]) + sorted(Sc, key=lambda i: -V[i])
            sg = sgn_perm([orden.index(i) for i in range(Rp)])
            if qui == -1:
                sg = -sg
            Y = tuple(libre)
            val[Y] = val.get(Y, 0) + sg * E * dv
            # sg lleva ya la quiralidad dentro; se guarda tambien qui aparte para poder separarlos
            fac[Y] = {"S": tuple(sorted(S)), "qui": qui, "delta": dv, "sg": sg,
                      "clases": tuple(sorted(plegar(v, t)[0] for v in A)),
                      "A": tuple(A)}
    val = {k: v for k, v in val.items() if v != 0}
    fac = {k: v for k, v in fac.items() if k in val}
    return val, fac, V


def progresiones(nu, t, r):
    """Para cada X, los terminos (k, Y_crudo, s, clave canonica)."""
    if not nu:
        return {}
    M = max(max(abs(v) for v in k) for k in nu)
    out = {}
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
        lst = []
        for k in itertools.product(*ks):
            Y = tuple(X[j] + t * k[j] for j in range(r))
            e = enderezar_D(Y)
            if e is None or e[0] not in nu:
                continue
            lst.append((k, Y, e[1], e[0]))
        if lst:
            out[X] = lst
    return out


print("=" * 104)
print(".ES EL PLEGADO AFIN?  DESCOMPOSICION DEL SIGNO EN LOS PARES QUE SE CANCELAN")
print("=" * 104)
print("")

A1 = A1n = 0
A4 = A4n = 0
A5 = A5n = 0
conj = Counter()
A6 = Counter()
mismos = Counter()
G1 = Counter()
A7 = Counter()
A7n = Counter()
testigos = []
rnd = random.Random(20260816)

for (t, r, cota) in CASOS:
    Rp = (t - 1) // 2 + r
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu, fac, V = nu_con_factores(list(Lam), t, r)
        if not nu:
            continue
        pr = progresiones(nu, t, r)
        # bolsa para el senuelo: todos los terminos de todas las progresiones, con su X
        bolsa = [(X, tm) for X, lst in pr.items() for tm in lst]
        for X, lst in pr.items():
            if len(lst) < 2:
                continue
            # A7  .de donde salen los pares que NO se cancelan?  Mi lectura es que viven todos en
            # las progresiones de CUATRO terminos, donde el emparejamiento es otro.  Se mide.
            npar = ncan = 0
            for p in range(len(lst)):
                for q in range(p + 1, len(lst)):
                    f1, f2 = fac[lst[p][3]], fac[lst[q][3]]
                    if lst[p][3] == lst[q][3]:
                        continue
                    npar += 1
                    if (lst[p][2] * f1["sg"] * f1["delta"]
                            == -lst[q][2] * f2["sg"] * f2["delta"]):
                        ncan += 1
            A7n[len(lst)] += npar
            A7[len(lst)] += npar - ncan
            for p in range(len(lst)):
                for q in range(p + 1, len(lst)):
                    (k1, Y1, s1, c1) = lst[p]
                    (k2, Y2, s2, c2) = lst[q]
                    f1, f2 = fac[c1], fac[c2]
                    if c1 == c2:
                        mismos["misma clave canonica"] += 1
                        continue
                    A1n += 1
                    rs = s2 // s1
                    rsg = f2["sg"] // f1["sg"]
                    rd = f2["delta"] // f1["delta"]
                    conj[(rs, rsg, rd)] += 1
                    if rd == -1 and rs == 1 and rsg == 1:
                        A1 += 1
                    # A6  los invariantes de verdad no son mis tres factores: s y sg dependen de
                    # convenios (enderezado, orden del barajado) y su PRODUCTO es el determinante de
                    # Weyl, que no depende de nada.  El otro invariante es delta = tau(eta_w).
                    A6[(rs * rsg, rd)] += 1
                    mismos[("S igual" if f1["S"] == f2["S"] else "S distinto") + " / " +
                           ("qui igual" if f1["qui"] == f2["qui"] else "qui distinto")] += 1
                    A4n += 1
                    if f1["clases"] == f2["clases"]:
                        A4 += 1
                    A5n += 1
                    if s1 * f1["sg"] * f1["delta"] == -s2 * f2["sg"] * f2["delta"]:
                        A5 += 1
                    elif len(testigos) < 3:
                        testigos.append({"t": t, "r": r, "Lambda": list(Lam), "X": list(X),
                                         "V": V,
                                         "term1": {"k": list(k1), "clave": list(c1), "s": s1,
                                                   "S": list(f1["S"]), "qui": f1["qui"],
                                                   "sg": f1["sg"], "delta": f1["delta"],
                                                   "A": list(f1["A"])},
                                         "term2": {"k": list(k2), "clave": list(c2), "s": s2,
                                                   "S": list(f2["S"]), "qui": f2["qui"],
                                                   "sg": f2["sg"], "delta": f2["delta"],
                                                   "A": list(f2["A"])}})
        # G1  senuelo: pares de progresiones DISTINTAS
        if len(bolsa) >= 2:
            for _ in range(min(40, len(bolsa))):
                (Xa, ta), (Xb, tb) = rnd.sample(bolsa, 2)
                if Xa == Xb or ta[3] == tb[3]:
                    continue
                fa, fb = fac[ta[3]], fac[tb[3]]
                G1[(ta[2] // tb[2], fa["sg"] // fb["sg"], fa["delta"] // fb["delta"])] += 1

print("  A1  el par voltea SOLO delta (plegado afin puro)      : %d de %d" % (A1, A1n))
print("")
print("  A2  reparto conjunto de (s'/s, sg'/sg, delta'/delta) sobre los pares de una progresion:")
for key, n in sorted(conj.items(), key=lambda kv: -kv[1]):
    print("        s'/s=%+d  sg'/sg=%+d  delta'/delta=%+d  :  %5d" % (key[0], key[1], key[2], n))
print("")
print("  A6  cruzado en los invariantes de verdad:  det(w) = s.sg   contra   delta = tau(eta_w)")
print("      (s y sg dependen de convenios; su producto y delta, no)")
n_uno = 0
for key, n in sorted(A6.items(), key=lambda kv: -kv[1]):
    etiqueta = {(-1, 1): "voltea SOLO det", (1, -1): "voltea SOLO delta",
                (1, 1): "no voltea nada", (-1, -1): "voltea los dos"}[key]
    print("        det'/det=%+d  delta'/delta=%+d  -> %-18s : %5d" % (key[0], key[1], etiqueta, n))
    if key in ((-1, 1), (1, -1)):
        n_uno += n
print("      pares donde voltea EXACTAMENTE UNO de los dos : %d  (y los que cancelan son %d)"
      % (n_uno, A5))
print("")
print("  A7  pares que NO cancelan, por tamano de la progresion:")
for n in sorted(A7n):
    print("        progresiones de %d terminos : %4d pares, %4d sin cancelar" % (n, A7n[n], A7[n]))
print("")
print("  A3  .que cambia entre los dos terminos?")
for key, n in sorted(mismos.items(), key=lambda kv: -kv[1]):
    print("        %-34s : %5d" % (key, n))
print("")
print("  A4  los dos bloques congelados tienen las mismas clases plegadas : %d de %d" % (A4, A4n))
print("      (trivial: los dos son transversales validos.  Se mide para no confundirlo con")
print("       contenido; si NO saliera entero, el fallo estaria en delta_dec.)")
print("")
print("  A5  el par se CANCELA (s.sg.delta opuesto)             : %d de %d" % (A5, A5n))
print("")
print("  G1  SENUELO: el mismo reparto entre progresiones DISTINTAS:")
tg = sum(G1.values())
for key, n in sorted(G1.items(), key=lambda kv: -kv[1])[:6]:
    print("        s=%+d sg=%+d delta=%+d : %5d  (%5.1f%%)" % (key[0], key[1], key[2], n,
                                                              100.0 * n / tg if tg else 0))
if testigos:
    print("")
    print("  !! pares que NO se cancelan (si los hay, A5 no es 100%%):")
    for tt in testigos:
        print("    " + json.dumps(tt)[:420])
print("")
print("  LECTURA: si A2 se concentra en UNA fila, esa fila ES la involucion.  Si la fila es")
print("  (+1, +1, -1), el mecanismo es exactamente el plegado afin del bloque congelado y (L1)")
print("  se prueba con delta.  Si el que voltea es sg, el mecanismo es de barajado y el plegado")
print("  afin no tiene nada que ver.")

json.dump({"A1": [A1, A1n], "A4": [A4, A4n], "A5": [A5, A5n],
           "A7_por_tamano": {str(k): [A7n[k], A7[k]] for k in sorted(A7n)},
           "A6_det_vs_delta": {str(k): v for k, v in sorted(A6.items(), key=lambda kv: -kv[1])},
           "A2_conjunto": {str(k): v for k, v in sorted(conj.items(), key=lambda kv: -kv[1])},
           "A3_que_cambia": {str(k): v for k, v in sorted(mismos.items(), key=lambda kv: -kv[1])},
           "G1_senuelo": {str(k): v for k, v in sorted(G1.items(), key=lambda kv: -kv[1])},
           "testigos": testigos},
          open("_probe_afin_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
