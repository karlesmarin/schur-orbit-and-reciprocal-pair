# -*- coding: utf-8 -*-
# LOS SUPERVIVIENTES SON EL GRUPO DE WEYL, Y tau ES SU CARACTER SIGNO.   16 de agosto de 2026.
#
# DE DONDE SALE.  Cruzar tres cosas que ya estan probadas y nunca se habian compuesto:
#
#   (a)  Cor. "el filtro es una funcion sobre (Z/t)^rango"  --  tau solo ve  a mod t.
#   (b)  Lema 3.1: si sobrevive, las clases  min(c_j, t-c_j)  son una PERMUTACION de {1,...,m},
#        y  tau = sgn(sigma) prod_j eps_j  con eps_j = +-1 segun de que lado del corte cae c_j.
#   (c)  el teorema de fusion: regular pliega al alcove con el signo del elemento que lo pliega.
#
# Componiendo: un patron de residuos superviviente ES el dato de (una permutacion de {1..m}, un
# vector de signos), o sea un elemento de  W(C_m) = {+-1}^m semidirecto S_m,  y tau es su SIGNO.
# Prediccion exacta, escrita antes de correr:
#
#     numero de clases de residuo supervivientes  =  |W| = 2^m m!   de las  t^m  posibles,
#     y  tau  toma  +1  en la mitad y  -1  en la otra mitad,  |W|/2 cada una.
#
# En el impar lo mismo con  W(B_m') = 2^{m'} m'!  sobre los residuos de  A_j = 2(eta_j + rho_j).
#
# Si sale, el conteo de supervivientes deja de ser un dato y pasa a ser el orden de un grupo, y el
# signo deja de ser una formula para ser un caracter.
#
# LO QUE SE MIDE
#   W1  cuantas clases de residuo (c_1..c_m) mod t sobreviven, contra 2^m m!.
#   W2  el reparto de tau entre +1 y -1 sobre esas clases, contra |W|/2 y |W|/2.
#   W3  la biyeccion explicita: para cada clase superviviente se extrae (sigma, eps) y se comprueba
#       que  tau = sgn(sigma) prod eps.
#
# CONTROLES
#   C0  tau se calcula por Freudenthal sobre un eta DOMINANTE de cada clase, no por la formula del
#       signo: si la formula y el caracter discrepan, es la formula la que esta mal.
#   C1  se comprueba que cada clase superviviente TIENE representante dominante; si alguna no lo
#       tuviera, el conteo de clases no seria el conteo de eta y habria que decirlo.
#   C2  SEÑUELO: el mismo conteo con las clases modulo t-1.  No tiene que dar 2^m m!.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage survivors_are_weyl.sage

import json
import itertools
from collections import Counter

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


def tau(typ, rk, eta, t):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car(typ, rk, eta).items():
        s += mult * z ** (sum((i + 1) * wt[i] for i in range(rk)) % t)
    return QQ(s) if s in QQ else None


print("=" * 116)
print("¿SON LOS SUPERVIVIENTES EL GRUPO DE WEYL, Y tau SU CARACTER SIGNO?")
print("=" * 116)
print("   t | G    |  t^rk | clases que sobreviven | 2^rk rk! | reparto (+1,-1) | signo == caracter |"
      " SEÑUELO t-1")
print("   " + "-" * 110)

RES = []
for t in range(3, 11):
    if t % 2 == 0:
        typ, rk = "C", (t - 2) // 2
        # a_j = eta_j + m - j + 1  (enteros)
        def shift(eta, j):
            return int(eta[j]) + rk - j
        doble = 1
    else:
        typ, rk = "B", (t - 1) // 2
        # A_j = 2 eta_j + 2(m'-j) + 1
        def shift(eta, j):
            return 2 * int(eta[j]) + 2 * (rk - j - 1) + 1
        doble = 2
    if rk < 1:
        continue

    # todos los eta dominantes de una caja generosa, agrupados por su clase de residuos
    clases = {}
    for k in range(0, 4 * t + 1):
        for e in Partitions(k, max_length=rk):
            eta = tuple(list(e) + [0] * (rk - len(e)))
            c = tuple(shift(eta, j) % t for j in range(rk))
            if c not in clases:
                clases[c] = eta

    vivos = {}
    for c, eta in clases.items():
        v = tau(typ, rk, eta, t)
        if v is not None and v != 0:
            vivos[c] = v

    W_orden = 2 ** rk * factorial(rk)
    reparto = Counter(vivos.values())

    # W3: la biyeccion (sigma, eps) y el signo predicho
    ok_signo = 0
    for c, v in vivos.items():
        cl = [min(x, t - x) for x in c]
        # sigma: el orden en que aparecen las clases 1..rk;  eps: de que lado del corte
        eps = [1 if x <= (t // 2) else -1 for x in c]
        if doble == 2:
            # en tipo B el corte natural es sobre A_j, y las clases son 1..rk en el mismo sentido
            eps = [1 if x <= (t - 1) // 2 else -1 for x in c]
        perm = Permutation([cl.index(i + 1) + 1 for i in range(rk)]) if sorted(cl) == list(range(1, rk + 1)) else None
        if perm is None:
            continue
        pred = perm.signature() * prod(eps)
        if pred == v:
            ok_signo += 1

    # SEÑUELO: clases modulo t-1
    cl2 = {}
    for k in range(0, 4 * t + 1):
        for e in Partitions(k, max_length=rk):
            eta = tuple(list(e) + [0] * (rk - len(e)))
            c = tuple(shift(eta, j) % (t - 1) for j in range(rk))
            if c not in cl2:
                cl2[c] = eta
    viv2 = sum(1 for c, eta in cl2.items() if (tau(typ, rk, eta, t) or 0) != 0)

    print("   %2d | %s%-2d  | %5d | %21d | %8d | %-15s | %-17s | %d"
          % (t, typ, rk, t ** rk, len(vivos), W_orden,
             "%d, %d" % (reparto.get(1, 0), reparto.get(-1, 0)),
             "%d de %d" % (ok_signo, len(vivos)), viv2))
    RES.append({"t": int(t), "tipo": typ, "rango": int(rk), "clases_totales": int(t ** rk),
                "clases_vivas": int(len(vivos)), "orden_W": int(W_orden),
                "mas1": int(reparto.get(1, 0)), "menos1": int(reparto.get(-1, 0)),
                "signo_ok": int(ok_signo), "senuelo": int(viv2)})

print("")
print("=" * 116)
print("  LECTURA, escrita ANTES de correr:")
print("   * si las clases vivas son exactamente 2^rk rk! = |W| y el reparto es mitad y mitad, el")
print("     conteo de supervivientes ES el orden del grupo de Weyl y tau ES su caracter signo.")
print("   * si el señuelo mod t-1 da otro numero, la periodicidad es genuinamente t.")
json.dump(RES, open("survivors_are_weyl_DUMP.json", "w"), indent=1)
print("=" * 116)
print("DONE")
