# -*- coding: utf-8 -*-
# LOS TESTIGOS DE LA VUELTA 28, verbatim.   16 de agosto de 2026.
#
# Hay dos contraejemplos concretos y el paper los cita ahora como \stverif.  Un
# numero citado que no se ha corrido es exactamente lo que no se hace, asi que se corren aparte,
# rapido, sin depender de la corrida grande de galois_sign.sage.
#
#   (i)   t=5, n=2, eta=(0,0):   tau^B_5(0,0) = +1   pero   tau~^C_5(2,2) = -1.
#         Los LOCUS coinciden; los valores no.  Esto es lo que mata "the two filter rules differ
#         by nothing else".
#   (ii)  t=4, m=1:   tau(a=1) = +1   pero   tau(a=3) = -1,  y 3 es unidad mod 4.
#         Esto es lo que mata "una evaluacion por orbita y el filtro queda conocido".
#
# Y de propina, dos controles nuestros que el paper afirma y conviene tener sueltos:
#   (iii) la constante de normalizacion de tipo B,  eps_t,  en t = 3, 5, 7.
#   (iv)  los cardinales de R^B y R^C en (6,2), que sostienen 3.6(iii): 0 y 8.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage v28_witnesses.sage

import json
import sys
import itertools

def car(typ, rk, mu):
    W = WeylCharacterRing("%s%d" % (typ, rk))
    el = W(W.space().from_vector(vector([Integer(v) for v in mu])))
    d = {}
    for wt, mult in el.weight_multiplicities().items():
        k = tuple(int(v) for v in wt.to_vector())
        d[k] = d.get(k, 0) + int(mult)
    return d


def tau_eval(typ, rk, eta, t):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car(typ, rk, eta).items():
        s += mult * z ** (sum((i + 1) * int(wt[i]) for i in range(rk)) % t)
    return int(QQ(s)) if s in QQ else None


def plegar(v, t):
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def delta_de(a, rk, t):
    cl, ep = [], []
    for v in a:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, rk + 1)):
        return 0
    perm = [cl[i] - 1 for i in range(rk)]
    sg = 1
    visto = [False] * rk
    for i in range(rk):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            sg = -sg
    for e in ep:
        sg *= e
    return int(sg)


OUT = {}
print("=" * 110)
print("LOS TESTIGOS DE LA VUELTA 28")
print("=" * 110)

# (i)
b = tau_eval("B", 2, (0, 0), 5)
c = tau_eval("C", 2, (2, 2), 5)
print("")
print("  (i)  t=5, n=2 :  tau^B_5(0,0) = %+d   |   tau~^C_5(2,2) = %+d" % (b, c))
print("       -> el soporte se traslada por (t-1)/2 = 2 ; el VALOR cambia de signo.")
print("       (lo predicho es  +1  y  -1)")
OUT["i"] = {"tauB_00": int(b), "tauC_22": int(c), "predicho": [int(1), int(-1)],
            "coincide": bool(int(b) == 1 and int(c) == -1)}

# (ii)   a = eta_1 + 1  en C_1  ->  a=1 es eta=(0),  a=3 es eta=(2)
w1 = tau_eval("C", 1, (0,), 4)
w2 = tau_eval("C", 1, (2,), 4)
print("")
print("  (ii) t=4, m=1 :  tau(a=1) = %+d   |   tau(a=3) = %+d     (3 es unidad mod 4)" % (w1, w2))
print("       -> el soporte es invariante bajo unidades; el VALOR no.")
print("       (lo predicho es  +1  y  -1)")
OUT["ii"] = {"tau_a1": int(w1), "tau_a3": int(w2), "predicho": [int(1), int(-1)],
             "coincide": bool(int(w1) == 1 and int(w2) == -1)}

# (iii)  la constante eps_t de tipo B
print("")
print("  (iii) la constante de normalizacion de tipo B :  tau = eps_t . delta")
OUT["iii"] = {}
for t in [3, 5, 7, 9]:
    rk = (t - 1) // 2
    eps = None
    n_ok = n = 0
    for eta in itertools.product(range(3), repeat=rk):
        if any(eta[j] < eta[j + 1] for j in range(rk - 1)):
            continue
        tv = tau_eval("B", rk, eta, t)
        a0 = [2 * int(eta[j]) + 2 * (rk - j) - 1 for j in range(rk)]
        d0 = delta_de(a0, rk, t)
        if eps is None and tv and d0:
            eps = int(tv) * int(d0)
    for eta in itertools.product(range(3), repeat=rk):
        if any(eta[j] < eta[j + 1] for j in range(rk - 1)):
            continue
        tv = tau_eval("B", rk, eta, t)
        a0 = [2 * int(eta[j]) + 2 * (rk - j) - 1 for j in range(rk)]
        d0 = delta_de(a0, rk, t)
        n += 1
        if eps is not None and d0 * eps == tv:
            n_ok += 1
    print("        t=%2d  B_%d :  eps_t = %+d    control  %d de %d" % (t, rk, eps, n_ok, n))
    OUT["iii"][str(t)] = {"eps_t": (int(eps) if eps is not None else None),
                          "control": [int(n_ok), int(n)]}

# (iv)  los cardinales que sostienen 3.6(iii)
print("")
print("  (iv) |R^B| y |R^C| en los casos pares que el paper cita")
OUT["iv"] = []
for (t, n) in [(6, 2), (8, 2), (10, 2), (12, 2), (14, 2), (8, 3), (10, 3), (12, 3), (14, 3)]:
    pts = itertools.product(range(t), repeat=n)
    RB = RC = 0
    for x in pts:
        pares = all((x[i] - x[j]) % t != 0 and (x[i] + x[j]) % t != 0
                    for i in range(n) for j in range(i + 1, n))
        if not pares:
            continue
        if all(v % t != 0 for v in x):
            RB += 1
        if all((2 * v) % t != 0 for v in x):
            RC += 1
    print("        (t,n) = (%2d,%d) :  |R^B| = %4d   |R^C| = %4d" % (t, n, RB, RC))
    OUT["iv"].append({"t": int(t), "n": int(n), "RB": int(RB), "RC": int(RC),
                      "distintos": bool(RB != RC)})

json.dump(OUT, open("v28_witnesses_DUMP.json", "w"), indent=1)
print("")
print("=" * 110)
print("DONE")
