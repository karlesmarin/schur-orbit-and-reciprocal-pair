# -*- coding: utf-8 -*-
# .CUANTA CANCELACION ES SOLO ANTISIMETRIA DE WEYL?
#
# _probe_involucion.py mato la involucion por intercambio de transversales (I3 166/266, I4 197/333),
# y de paso enseno el camino bueno: 47 pares de contribuyentes tenian el MISMO transversal.  Mirando
# un testigo a mano, t=3, X=(4,-2):
#
#       k=(1,5) -> Y=(7,13)        k=(3,3) -> Y=(13,7)
#
# es el MISMO punto del soporte alcanzado dos veces por la progresion, en dos ordenes distintos.  Al
# enderezar por W(D_r) uno lleva signo + y el otro signo -, y se cancelan sin que intervenga nada
# del transversal.  Pura antisimetria.
#
# Entonces la pregunta correcta no es "que involucion hay sobre los transversales" sino:
#
#       .se cancela TODO por esa via, o queda un residuo que necesita otro mecanismo?
#
# LO QUE SE MIDE, agrupando los terminos de cada progresion por su clave CANONICA:
#   W0  cuantas claves distintas toca una progresion.
#   W1  progresiones que se cancelan enteras DENTRO de cada clave (antisimetria pura).
#   W2  el residuo: progresiones donde alguna clave sobrevive.  Ahi hay contenido de verdad.
#   W3  tras quitar la cancelacion trivial, .queda a lo sumo UN termino?  Si si, (L1) esta probada:
#       c(X) = nu de ese termino, y |c| <= 1 se hereda de nu en {0,+-1}.
#
# SENUELO
#   F1  la misma agrupacion pero SIN el signo del enderezado: la cancelacion trivial tiene que
#       desaparecer casi entera, o no estoy midiendo la antisimetria sino el emparejamiento.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_weyl_doble.py

import itertools
import json
from collections import Counter, defaultdict

from divided_differences import CASOS, nu_de, enderezar_D


def terminos(nu, t, r):
    """Para cada X, la lista de terminos (k, clave canonica, signo del enderezado, valor de nu)."""
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
            if e is None:
                continue
            v = nu.get(e[0], 0)
            if v:
                lst.append((k, e[0], e[1], v))
        if lst:
            out[X] = lst
    return out


print("=" * 104)
print("CUANTA CANCELACION ES SOLO ANTISIMETRIA DE WEYL")
print("=" * 104)
print("")

W0 = Counter()
tot = trivial = residuo = 0
W3 = W3n = 0
resto_hist = Counter()
F1 = F1n = 0
testigos = []

for (t, r, cota) in CASOS:
    Rp = (t - 1) // 2 + r
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu = nu_de(list(Lam), t, r)
        if not nu:
            continue
        for X, lst in terminos(nu, t, r).items():
            if len(lst) < 2:
                continue
            tot += 1
            por_clave = defaultdict(lambda: 0)
            for (k, clave, s, v) in lst:
                por_clave[clave] += s * v
            W0[len(por_clave)] += 1
            vivos = {c: x for c, x in por_clave.items() if x != 0}
            if not vivos:
                trivial += 1
            else:
                residuo += 1
                if len(testigos) < 4:
                    testigos.append({"t": t, "r": r, "Lambda": list(Lam), "X": list(X),
                                     "terminos": [[list(k), list(c), s, v] for (k, c, s, v) in lst],
                                     "residuo_por_clave": {str(c): x for c, x in vivos.items()}})
            resto_hist[len(vivos)] += 1

            # F1  senuelo: sin el signo del enderezado
            pc = defaultdict(lambda: 0)
            for (k, clave, s, v) in lst:
                pc[clave] += v
            F1n += 1
            if all(x == 0 for x in pc.values()):
                F1 += 1

            # W3  tras quitar lo trivial, .queda a lo sumo un termino?
            W3n += 1
            if len(vivos) <= 1:
                W3 += 1

print("  progresiones con 2 o mas terminos                    : %d" % tot)
print("  claves canonicas distintas que toca cada una         : %s" % dict(sorted(W0.items())))
print("")
print("  W1  se cancelan ENTERAS dentro de cada clave         : %d de %d" % (trivial, tot))
print("  W2  dejan residuo en alguna clave                    : %d de %d" % (residuo, tot))
print("      claves con residuo no nulo : %s" % dict(sorted(resto_hist.items())))
print("  W3  tras la cancelacion trivial queda <= 1 termino    : %d de %d" % (W3, W3n))
print("")
print("  F1  SENUELO: lo mismo SIN el signo del enderezado     : %d de %d" % (F1, F1n))
print("      (si F1 fuera tan alto como W1, no estaria midiendo la antisimetria)")
if testigos:
    print("")
    print("  primeros casos CON residuo (el contenido que queda):")
    for tt in testigos:
        print("    t=%d r=%d Lambda=%s X=%s" % (tt["t"], tt["r"], tt["Lambda"], tt["X"]))
        for (k, c, s, v) in tt["terminos"]:
            print("        k=%-10s clave=%-12s signo=%+d  nu=%+d" % (str(k), str(c), s, v))
        print("        residuo: %s" % tt["residuo_por_clave"])
print("")
print("  LECTURA: si W1 == tot, toda la cancelacion es antisimetria y (L1) sale sola.  Si queda")
print("  residuo, ESE es el objeto que falta entender, y ya no es la division: es el soporte de nu.")

json.dump({"tot": tot, "trivial": trivial, "residuo": residuo,
           "claves": {str(k): v for k, v in sorted(W0.items())},
           "resto_hist": {str(k): v for k, v in sorted(resto_hist.items())},
           "W3": [W3, W3n], "F1_senuelo": [F1, F1n], "testigos": testigos},
          open("_probe_weyl_doble_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
