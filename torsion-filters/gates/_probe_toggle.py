# -*- coding: utf-8 -*-
# EL TOGGLE ESTRECHO: NO INTERCAMBIOS ARBITRARIOS, SINO EL PASO DENTRO DE UNA CLASE.
#
# QUE SE PRUEBA, Y POR QUE ESE Y NO OTRO.  Los intercambios arbitrarios ya estan descartados.  El
# candidato estrecho es: sustituir, en una UNICA clase residual, el elemento elegido por el
# siguiente de esa misma clase plegada, cuando sus doubled shifted exponents difieren en 2t.
#
# Y la lectura ingenua se cae antes de empezar, por una observacion de dos lineas: en el grupo
# afin de A_1 con paredes en 0 y t,  s_0(x) = -x,  s_1(x) = 2t - x,  luego  s_1 s_0 (x) = x + 2t.
# La traslacion por 2t es un producto de DOS reflexiones, asi que su signo afin es +1 y NO puede
# aportar el -1 por si sola.  El -1 tiene que venir de otro sitio: del representante w, del
# enderezado por W(D_r) o de delta_S.
#
# LO QUE SE FALSA, sobre los pares que de verdad se cancelan:
#   T1  el par ES el toggle exacto: S' = S con un indice a cambiado por b, misma clase plegada,
#       |V_a - V_b| = 2t EXACTAMENTE, y la misma quiralidad.
#   T2  lo mismo relajando a |V_a - V_b| multiplo de 2t.
#   T3  lo mismo relajando a multiplo de t (para ver si el 2 aporta algo o no).
#   T4  los pares que NO son toggle, .que son?  Reparto: quiralidad pura, |S \ S'| > 2, etc.
#   T5  y el reciproco, que es lo que haria de esto una involucion: aplicar el toggle a un objeto
#       contribuyente, .cae en la MISMA progresion y con signo opuesto?
#
# SENUELO
#   J1  el toggle con paso 3t (impar), que por la observacion afin de arriba deberia comportarse
#       distinto.  Si acertara igual, el 2t no seria el mecanismo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_toggle.py

import itertools
import json
from collections import Counter

from divided_differences import (CASOS, plegar, sgn_perm, eps_t, delta_dec, enderezar_D)


def objetos(Lam, t, r):
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = {}
    for i, x in enumerate(V):
        d.setdefault(plegar(x, t)[0], []).append(i)
    if any(not d.get(j) for j in range(1, mp + 1)):
        return {}, {}, V, d
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
            obj[Y] = {"S": frozenset(S), "qui": qui, "delta": dv, "sg": sg}
    val = {k: v for k, v in val.items() if v != 0}
    obj = {k: v for k, v in obj.items() if k in val}
    return val, obj, V, d


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


def es_toggle(S1, S2, V, t, paso, exacto):
    """.S2 sale de S1 cambiando UN indice por otro de la misma clase plegada a distancia paso?"""
    dif1, dif2 = S1 - S2, S2 - S1
    if len(dif1) != 1 or len(dif2) != 1:
        return False, None
    a, b = next(iter(dif1)), next(iter(dif2))
    if plegar(V[a], t)[0] != plegar(V[b], t)[0]:
        return False, None
    dv = abs(V[a] - V[b])
    if exacto:
        return (dv == paso), (a, b)
    return (dv % paso == 0), (a, b)


print("=" * 104)
print("EL TOGGLE ESTRECHO: cambiar UN elemento por el de su clase a distancia 2t")
print("=" * 104)
print("")

T1 = T2 = T3 = J1 = 0
Tn = 0
resto = Counter()
T5 = T5n = 0
T6 = [0]
T6n = [0]
T7 = [0]
T7n = [0]
forma = Counter()
signos_toggle = Counter()
testigos = []

for (t, r, cota) in CASOS:
    Rp = (t - 1) // 2 + r
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu, obj, V, clases = objetos(list(Lam), t, r)
        if not nu:
            continue
        for X, lst in progresiones(nu, t, r).items():
            if len(lst) < 2:
                continue
            for p in range(len(lst)):
                for q in range(p + 1, len(lst)):
                    (k1, s1, c1) = lst[p]
                    (k2, s2, c2) = lst[q]
                    if c1 == c2:
                        continue
                    if s1 * nu[c1] != -s2 * nu[c2]:
                        continue                    # solo los pares que cancelan
                    o1, o2 = obj[c1], obj[c2]
                    Tn += 1
                    mismo_qui = o1["qui"] == o2["qui"]
                    ok1, par = es_toggle(o1["S"], o2["S"], V, t, 2 * t, True)
                    ok2, _ = es_toggle(o1["S"], o2["S"], V, t, 2 * t, False)
                    ok3, _ = es_toggle(o1["S"], o2["S"], V, t, t, False)
                    okJ, _ = es_toggle(o1["S"], o2["S"], V, t, 3 * t, True)
                    T1 += 1 if (ok1 and mismo_qui) else 0
                    T2 += 1 if (ok2 and mismo_qui) else 0
                    T3 += 1 if (ok3 and mismo_qui) else 0
                    J1 += 1 if (okJ and mismo_qui) else 0
                    # T6  el toggle AMPLIADO: un intercambio dentro de una sola clase plegada
                    # (a distancia multiplo de t) y/o un cambio de quiralidad.  Se mide entero, en
                    # vez de sumar a mano las categorias de T4.
                    T6n[0] += 1
                    dif = o1["S"] ^ o2["S"]
                    intercambio_ok = (len(dif) == 0) or ok3
                    if intercambio_ok:
                        T6[0] += 1
                    # T7  el toggle COMPLETO.  Hipotesis: los 58 que fallan T6 son intercambios
                    # dentro de la misma clase plegada pero con el SIGNO DE PLEGADO OPUESTO (y
                    # entonces V_a - V_b no es multiplo de t, sino congruente a 2j), compensados por
                    # el cambio de quiralidad.  Si es asi, las tres piezas cierran 313 de 313.
                    T7n[0] += 1
                    if len(dif) == 0:
                        if not mismo_qui:
                            T7[0] += 1
                            forma["quiralidad pura"] += 1
                    elif len(dif) == 2:
                        a_ = next(iter(o1["S"] - o2["S"]))
                        b_ = next(iter(o2["S"] - o1["S"]))
                        ca, ea = plegar(V[a_], t)
                        cb, eb = plegar(V[b_], t)
                        if ca == cb:
                            if ea == eb and mismo_qui:
                                T7[0] += 1
                                forma["misma clase, mismo signo de plegado, misma quiralidad"] += 1
                            elif ea != eb and not mismo_qui:
                                T7[0] += 1
                                forma["misma clase, signo de plegado OPUESTO, quiralidad volteada"] += 1
                            else:
                                forma["misma clase pero signos descoordinados"] += 1
                        else:
                            forma["clase distinta"] += 1
                    else:
                        forma["mas de un intercambio"] += 1
                    if ok1 and mismo_qui:
                        signos_toggle[(o2["delta"] // o1["delta"],
                                       (o2["sg"] // o1["sg"]) * (s2 // s1))] += 1
                    else:
                        if o1["S"] == o2["S"]:
                            resto["quiralidad pura (S igual)"] += 1
                        elif len(o1["S"] ^ o2["S"]) == 2 and not mismo_qui:
                            resto["un cambio PERO quiralidad distinta"] += 1
                        elif len(o1["S"] ^ o2["S"]) == 2:
                            resto["un cambio, distancia no 2t"] += 1
                            if len(testigos) < 3:
                                a = next(iter(o1["S"] - o2["S"]))
                                b = next(iter(o2["S"] - o1["S"]))
                                testigos.append({"t": t, "Lambda": list(Lam), "X": list(X),
                                                 "V": V, "V_a": V[a], "V_b": V[b],
                                                 "dif": V[a] - V[b], "2t": 2 * t})
                        else:
                            resto["|S \\ S'| = %d" % (len(o1["S"] ^ o2["S"]) // 2)] += 1

print("  sobre los %d pares que se cancelan (puntos distintos):" % Tn)
print("")
print("  T1  ES el toggle exacto (|V_a - V_b| = 2t, misma quiralidad) : %d de %d" % (T1, Tn))
print("  T2  relajado a multiplo de 2t                                : %d de %d" % (T2, Tn))
print("  T3  relajado a multiplo de t                                 : %d de %d" % (T3, Tn))
print("")
print("  J1  SENUELO: el mismo toggle con paso 3t                     : %d de %d" % (J1, Tn))
print("")
print("  T6  el toggle AMPLIADO (un intercambio en UNA clase a distancia multiplo de t)")
print("      : %d de %d" % (T6[0], T6n[0]))
print("")
print("  T7  el toggle COMPLETO, con la quiralidad ACOPLADA al signo de plegado : %d de %d"
      % (T7[0], T7n[0]))
for k, n in sorted(forma.items(), key=lambda kv: -kv[1]):
    print("        %-52s : %5d" % (k, n))
print("")
print("  NOTA sobre J1: los V_i son todos IMPARES, luego V_a - V_b es PAR; si ademas es multiplo")
print("      de t (impar), es multiplo de 2t.  Por eso T2 == T3 y por eso el paso 3t da 0 por")
print("      construccion.  El senuelo confirma la aritmetica, no discrimina hipotesis.")
print("")
print("  T4  los que NO son el toggle exacto, .que son?")
for k, n in sorted(resto.items(), key=lambda kv: -kv[1]):
    print("        %-34s : %5d" % (k, n))
print("")
print("  y en los que SI son el toggle, de donde sale el -1:")
for k, n in sorted(signos_toggle.items(), key=lambda kv: -kv[1]):
    print("        delta'/delta=%+d   (s.sg)'/(s.sg)=%+d  : %5d" % (k[0], k[1], n))
if testigos:
    print("")
    print("  testigos de 'un cambio pero la distancia no es 2t':")
    for tt in testigos:
        print("    t=%d Lambda=%s V=%s : V_a=%d V_b=%d dif=%d (2t=%d)"
              % (tt["t"], tt["Lambda"], tt["V"], tt["V_a"], tt["V_b"], tt["dif"], tt["2t"]))
print("")
print("  LECTURA: si T1 sale entero, la involucion es el toggle y la cancelacion se demuestra sin")
print("  Kac-Walton.  Si no, T4 dice exactamente que otras piezas hacen falta -- y la observacion")
print("  afin (2t = producto de dos reflexiones, signo +1) dice que el -1 no puede venir del salto.")

json.dump({"n_pares": Tn, "T1": T1, "T2": T2, "T3": T3, "J1_senuelo": J1,
           "T6_ampliado": [T6[0], T6n[0]], "T7_completo": [T7[0], T7n[0]],
           "forma": {str(k): v for k, v in sorted(forma.items(), key=lambda kv: -kv[1])}, "resto": {str(k): v for k, v in sorted(resto.items(), key=lambda kv: -kv[1])},
           "signos_en_el_toggle": {str(k): v for k, v in sorted(signos_toggle.items())},
           "testigos": testigos},
          open("_probe_toggle_DUMP.json", "w"), indent=1)
print("")
print("=" * 104)
print("DONE")
