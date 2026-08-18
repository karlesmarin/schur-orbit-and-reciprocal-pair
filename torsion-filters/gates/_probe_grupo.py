# -*- coding: utf-8 -*-
# OTRA DESCOMPOSICION DE W^1: LEER EL ELEMENTO DE GRUPO, NO ADIVINAR EL FACTOR.
#
# Otra descomposicion de W^1.  Pero la manera de probar otra no
# es proponer una tercera factorizacion del signo a ver si suena.  Cada objeto contribuyente ES una
# permutacion con signo del vector V:
#
#     w(V) = ( bloque congelado ordenado desc. ; bloque libre ordenado desc. con la quiralidad )
#
# y las dos entradas de un par que se cancela son dos reordenaciones con signo del MISMO multiconjunto
# V.  Luego existe un elemento honesto
#
#     w'' = w' . w^{-1}   en   W(B_{R'})
#
# que lleva una a la otra, sin que yo elija nada.  Nunca lo he mirado: en _probe_involucion.py medi
# solo |S \ S'|, que es la sombra de w'' sobre los subconjuntos y pierde el grupo entero.
#
# LO QUE SE FALSA
#   B1  w'' es una INVOLUCION (orden 2).  Si lo es, "la involucion" existe literalmente y solo hay
#       que nombrarla.
#   B2  w'' es una REFLEXION de B_{R'}:  s_{e_i-e_j} (transposicion),  s_{e_i+e_j} (transposicion con
#       los dos signos cambiados)  o  s_{e_i} (un cambio de signo).  Son las tres clases posibles.
#   B3  el tipo de ciclo de w'' -- si se concentra en uno solo, ESE es el mecanismo, se llame como se
#       llame, y ya no depende de como yo parta el signo.
#   B4  .los dos indices que mueve w'' estan en la misma clase plegada?  Es lo unico que puede
#       explicar que el bloque libre se desplace un multiplo de t.
#   B5  det(w'') contra el cociente de signos de nu: tienen que coincidir salvo el factor delta.
#
# SENUELO
#   H1  lo mismo sobre pares de progresiones DISTINTAS.  Si el tipo de ciclo se concentrara igual
#       ahi, no estaria midiendo la progresion sino la forma de W^1.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_grupo.py

import itertools
import json
import random
from collections import Counter

from divided_differences import (CASOS, plegar, sgn_perm, eps_t, delta_dec, enderezar_D)


def objetos(Lam, t, r):
    """Cada contribuyente, con su ARREGLO con signo de V y los datos que lo generan."""
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = {}
    for i, x in enumerate(V):
        d.setdefault(plegar(x, t)[0], []).append(i)
    if any(not d.get(j) for j in range(1, mp + 1)):
        return {}, {}, V
    E = eps_t(t, mp)
    val, obj = {}, {}
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
            # el arreglo: posicion -> (indice de V, signo).  orden ya es la lista de indices.
            signos = [1] * Rp
            if qui == -1:
                signos[Rp - 1] = -1
            obj[Y] = {"sigma": list(orden), "signos": signos, "S": tuple(sorted(S)),
                      "qui": qui, "delta": dv, "sg": sg}
    val = {k: v for k, v in val.items() if v != 0}
    obj = {k: v for k, v in obj.items() if k in val}
    return val, obj, V


def compone(o1, o2):
    """w'' = w' . w^{-1} sobre POSICIONES: p -> q con sigma'(q) = sigma(p).  Devuelve (perm, signos)."""
    Rp = len(o1["sigma"])
    donde = {idx: q for q, idx in enumerate(o2["sigma"])}
    perm = [donde[o1["sigma"][p]] for p in range(Rp)]
    sig = [o1["signos"][p] * o2["signos"][perm[p]] for p in range(Rp)]
    return perm, sig


def tipo_ciclo(perm):
    n, visto, tipos = len(perm), [False] * len(perm), []
    for i in range(n):
        if visto[i]:
            continue
        L, j = 0, i
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        tipos.append(L)
    return tuple(sorted(tipos, reverse=True))


def es_involucion(perm, sig):
    n = len(perm)
    for p in range(n):
        if perm[perm[p]] != p:
            return False
        if sig[p] * sig[perm[p]] != 1:
            return False
    return True


def clase_reflexion(perm, sig, Rp):
    """.Es una reflexion de B_n?  s_{e_i-e_j}, s_{e_i+e_j} o s_{e_i}."""
    movidos = [p for p in range(Rp) if perm[p] != p]
    negados = [p for p in range(Rp) if sig[p] == -1]
    if not movidos and len(negados) == 1:
        return "s_{e_i}  (cambio de signo)"
    if len(movidos) == 2 and perm[movidos[0]] == movidos[1] and perm[movidos[1]] == movidos[0]:
        if not negados:
            return "s_{e_i-e_j}  (transposicion)"
        if sorted(negados) == sorted(movidos):
            return "s_{e_i+e_j}  (transposicion con signo)"
    return "no es reflexion"


def progresiones(nu, t, r):
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
            lst.append((k, e[1], e[0]))
        if lst:
            out[X] = lst
    return out


print("=" * 104)
print("OTRA DESCOMPOSICION DE W^1: EL ELEMENTO w'' = w' . w^{-1} QUE UNE LOS DOS TERMINOS")
print("=" * 104)
print("")

B1 = B1n = B2 = B2n = B4 = B4n = B5 = B5n = 0
ciclos = Counter()
refl = Counter()
H1 = Counter()
mismaclase = Counter()
rnd = random.Random(20260816)

for (t, r, cota) in CASOS:
    Rp = (t - 1) // 2 + r
    mp = (t - 1) // 2
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu, obj, V = objetos(list(Lam), t, r)
        if not nu:
            continue
        pr = progresiones(nu, t, r)
        bolsa = [(X, tm) for X, lst in pr.items() for tm in lst]
        for X, lst in pr.items():
            if len(lst) < 2:
                continue
            for p in range(len(lst)):
                for q in range(p + 1, len(lst)):
                    (k1, s1, c1) = lst[p]
                    (k2, s2, c2) = lst[q]
                    if c1 == c2:
                        continue
                    o1, o2 = obj[c1], obj[c2]
                    # solo los pares que de verdad se cancelan: son los que un mecanismo debe explicar
                    if s1 * nu[c1] != -s2 * nu[c2]:
                        continue
                    perm, sig = compone(o1, o2)
                    B1n += 1
                    if es_involucion(perm, sig):
                        B1 += 1
                    cl = clase_reflexion(perm, sig, Rp)
                    refl[cl] += 1
                    B2n += 1
                    if cl != "no es reflexion":
                        B2 += 1
                    ciclos[(tipo_ciclo(perm), sum(1 for x in sig if x == -1))] += 1
                    movidos = [i for i in range(Rp) if perm[i] != i]
                    if len(movidos) == 2:
                        i1 = o1["sigma"][movidos[0]]
                        i2 = o1["sigma"][movidos[1]]
                        B4n += 1
                        cl1, cl2 = plegar(V[i1], t)[0], plegar(V[i2], t)[0]
                        if cl1 == cl2:
                            B4 += 1
                        mismaclase[(cl1 == cl2, cl1, cl2)] += 1
                    det = sgn_perm(perm)
                    for x in sig:
                        det *= x
                    B5n += 1
                    if det * (o2["delta"] // o1["delta"]) == -1:
                        B5 += 1
        if len(bolsa) >= 2:
            for _ in range(min(40, len(bolsa))):
                (Xa, ta), (Xb, tb) = rnd.sample(bolsa, 2)
                if Xa == Xb or ta[2] == tb[2]:
                    continue
                perm, sig = compone(obj[ta[2]], obj[tb[2]])
                H1[(tipo_ciclo(perm), sum(1 for x in sig if x == -1))] += 1

print("  B1  w'' es una INVOLUCION (orden 2)                : %d de %d" % (B1, B1n))
print("  B2  w'' es una REFLEXION de B_{R'}                 : %d de %d" % (B2, B2n))
print("      reparto:")
for k, n in sorted(refl.items(), key=lambda kv: -kv[1]):
    print("        %-34s : %5d" % (k, n))
print("")
print("  B3  tipo de ciclo de w'' y numero de signos cambiados:")
tc = sum(ciclos.values())
for k, n in sorted(ciclos.items(), key=lambda kv: -kv[1])[:10]:
    print("        ciclos %-14s signos- %d : %5d  (%5.1f%%)" % (str(k[0]), k[1], n,
                                                               100.0 * n / tc if tc else 0))
print("        tipos distintos : %d" % len(ciclos))
print("")
print("  B4  cuando w'' mueve dos posiciones, .misma clase plegada? : %d de %d" % (B4, B4n))
print("")
print("  B5  det(w'') . (delta'/delta) = -1                  : %d de %d" % (B5, B5n))
print("")
print("  H1  SENUELO: tipo de ciclo entre progresiones DISTINTAS:")
th = sum(H1.values())
for k, n in sorted(H1.items(), key=lambda kv: -kv[1])[:6]:
    print("        ciclos %-14s signos- %d : %5d  (%5.1f%%)" % (str(k[0]), k[1], n,
                                                               100.0 * n / th if th else 0))
print("        tipos distintos : %d" % len(H1))
print("")
print("  LECTURA: si B1 sale entero, la involucion existe literalmente y solo falta nombrarla.  Si")
print("  ademas B3 se concentra en un tipo, ese tipo ES el mecanismo, y ya no depende de como yo")
print("  parta el signo -- que es lo que fallaba en la descomposicion anterior.")

json.dump({"B1": [B1, B1n], "B2": [B2, B2n], "B4": [B4, B4n], "B5": [B5, B5n],
           "reflexiones": {str(k): v for k, v in sorted(refl.items(), key=lambda kv: -kv[1])},
           "ciclos": {str(k): v for k, v in sorted(ciclos.items(), key=lambda kv: -kv[1])},
           "H1_senuelo": {str(k): v for k, v in sorted(H1.items(), key=lambda kv: -kv[1])}},
          open("_probe_grupo_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
