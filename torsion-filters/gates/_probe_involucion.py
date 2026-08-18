# -*- coding: utf-8 -*-
# LA INVOLUCION QUE CANCELA.  .Que mueve un termino de la progresion al siguiente?
#
# _probe_cancelacion.py mato la hipotesis de la CAJA: 282 de 333, y el testigo {(1,5),(3,3)} se mueve
# en diagonal.  Asi que la estructura no es un producto de elecciones independientes.  Lo que queda
# en pie es lo que de verdad hace falta para (L1): una INVOLUCION LIBRE QUE INVIERTE EL SIGNO.
#
# Donde vive.  Por lem:muinj el mapa (S, quiralidad) -> peso es INYECTIVO, asi que cada punto del
# soporte de nu tiene UN solo origen.  Dos puntos de la misma progresion difieren en t*(vector par),
# lo que en el lenguaje de los V significa: una entrada del complemento pasa de V_a a V_b con
# V_a = V_b (mod 2t).  Es decir, S y S' se diferencian en un INTERCAMBIO a <-> b entre el transversal
# y su complemento, con a y b en la misma clase plegada y a distancia par de niveles.
#
# LO QUE SE FALSA
#   I0  cada punto del soporte tiene exactamente un origen (lem:muinj, aqui re-medido).
#   I1  dos contribuyentes de la misma progresion tienen |S \ S'| = 1 (un solo intercambio)
#       o se conectan por una cadena de intercambios asi.
#   I2  V_a = V_b (mod 2t) para el par intercambiado.   <-- por que la progresion es la MISMA
#   I3  un intercambio de esos INVIERTE el signo de nu.  <-- ESTA es la cancelacion
#   I4  el conjunto contribuyente es una sola orbita de esos intercambios (conexo).
#
# SENUELO
#   E1  el mismo test exigiendo congruencia modulo t en vez de 2t: tiene que admitir pares que NO
#       estan en la misma progresion, o sea I2 con t debe ser mas laxo que con 2t.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_involucion.py

import itertools
import json
from collections import Counter

from divided_differences import (CASOS, plegar, sgn_perm, eps_t, delta_dec, enderezar_D,
                                 nu_extendida)


def nu_con_origen(Lam, t, r):
    """nu, y para cada punto del soporte la lista de (S, quiralidad) que lo produjeron."""
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = {}
    for i, x in enumerate(V):
        d.setdefault(plegar(x, t)[0], []).append(i)
    if any(not d.get(j) for j in range(1, mp + 1)):
        return {}, {}, V
    E = eps_t(t, mp)
    val, org = {}, {}
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
            org.setdefault(Y, []).append((S, qui, sg * E * dv))
    val = {k: v for k, v in val.items() if v != 0}
    return val, org, V


def contribuyentes(nu, t, r):
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
        viv = {}
        for k in itertools.product(*ks):
            Y = tuple(X[j] + t * k[j] for j in range(r))
            if nu_extendida(nu, Y):
                e = enderezar_D(Y)
                viv[k] = e[0]                      # la clave canonica del soporte
        if viv:
            out[X] = viv
    return out


print("=" * 104)
print("LA INVOLUCION QUE CANCELA")
print("=" * 104)
print("")

I0 = I0n = I1 = I1n = I2 = I2n = I3 = I3n = I4 = I4n = 0
E1 = E1n = 0
distS = Counter()
fallos = []

for (t, r, cota) in CASOS:
    Rp = (t - 1) // 2 + r
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu, org, V = nu_con_origen(list(Lam), t, r)
        if not nu:
            continue
        for Y, lst in org.items():
            if Y in nu:
                I0n += 1
                if len(lst) == 1:
                    I0 += 1
        con = contribuyentes(nu, t, r)
        for X, viv in con.items():
            if len(viv) < 2:
                continue
            claves = sorted(set(viv.values()))
            objs = [org[c][0] for c in claves if c in org and len(org[c]) == 1]
            if len(objs) != len(claves):
                continue                       # sin origen unico no se puede hablar de intercambio
            # aristas: pares que difieren en un solo intercambio
            aristas = []
            for p in range(len(objs)):
                for q in range(p + 1, len(objs)):
                    Sp, Sq = objs[p][0], objs[q][0]
                    dif = Sp ^ Sq
                    distS[len(dif) // 2] += 1
                    I1n += 1
                    if len(dif) == 2:
                        I1 += 1
                        a = next(iter(Sp - Sq))
                        b = next(iter(Sq - Sp))
                        I2n += 1
                        if (V[a] - V[b]) % (2 * t) == 0:
                            I2 += 1
                        E1n += 1
                        if (V[a] - V[b]) % t == 0:
                            E1 += 1
                        I3n += 1
                        if objs[p][2] == -objs[q][2]:
                            I3 += 1
                            aristas.append((p, q))
                        elif len(fallos) < 3:
                            fallos.append({"que": "un intercambio NO invierte el signo",
                                           "t": t, "r": r, "Lambda": list(Lam), "X": list(X),
                                           "V": V, "a": a, "b": b,
                                           "signos": [objs[p][2], objs[q][2]]})
            # I4  conexo por esas aristas
            I4n += 1
            vistos = {0}
            frente = [0]
            while frente:
                x = frente.pop()
                for (p, q) in aristas:
                    for (u, v) in ((p, q), (q, p)):
                        if u == x and v not in vistos:
                            vistos.add(v)
                            frente.append(v)
            if len(vistos) == len(objs):
                I4 += 1
            elif len(fallos) < 4:
                fallos.append({"que": "conjunto contribuyente NO conexo por intercambios",
                               "t": t, "r": r, "Lambda": list(Lam), "X": list(X),
                               "n_objs": len(objs), "alcanzados": len(vistos)})

print("  I0  cada punto del soporte tiene UN origen (lem:muinj)      : %d de %d" % (I0, I0n))
print("  I1  dos contribuyentes difieren en un solo intercambio      : %d de %d" % (I1, I1n))
print("      reparto de |S \\ S'| : %s" % dict(sorted(distS.items())))
print("  I2  y el par intercambiado cumple V_a = V_b (mod 2t)        : %d de %d" % (I2, I2n))
print("  I3  un intercambio asi INVIERTE el signo  <== la cancelacion: %d de %d" % (I3, I3n))
print("  I4  el conjunto contribuyente es conexo por intercambios    : %d de %d" % (I4, I4n))
print("")
print("  E1  SENUELO: la misma congruencia pero modulo t             : %d de %d" % (E1, E1n))
print("      (si E1 == I2, la condicion 2t no esta aportando nada)")
if fallos:
    print("")
    print("  !! primeros fallos:")
    for f in fallos:
        print("    " + json.dumps(f)[:280])
print("")
print("  LECTURA: I3 + I4 son la demostracion.  Si cada intercambio invierte el signo y el conjunto")
print("  contribuyente es una sola orbita, la suma sobre la progresion es 0 en cuanto hay dos")
print("  terminos, y vale +-1 cuando hay uno.  Eso es (L1) entera.")

json.dump({"I0": [I0, I0n], "I1": [I1, I1n], "I2": [I2, I2n], "I3": [I3, I3n], "I4": [I4, I4n],
           "E1_senuelo": [E1, E1n], "dist_intercambios": {str(k): v for k, v in sorted(distS.items())},
           "fallos": fallos},
          open("_probe_involucion_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
