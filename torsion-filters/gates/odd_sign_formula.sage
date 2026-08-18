# -*- coding: utf-8 -*-
# EL SIGNO DEL FILTRO IMPAR, EN FORMA CERRADA.   16 de agosto de 2026.
#
# POR QUE HACE FALTA AHORA.  La Figura 2 hay que redibujarla con la regla ORTOGONAL (su vuelta 25),
# y el panel colorea por el SIGNO de tau.  La Prop. 5.3 dice cuando tau^B != 0 y que |tau^B| = 1,
# pero NO da el signo.  Antes de dibujarlo hay que saberlo, y con la formula correcta: si dibujamos
# el signo por la formula del caso par, repetimos exactamente el defecto que el acaba de cazar.
#
# LA CANDIDATA, por analogia con el Lema 3.1 y con el plegado afin: con
#
#     A_j = 2 eta_j + 2(m'-j) + 1,   c_j = A_j mod t,   cl_j = min(c_j, t-c_j),
#     eps_j = +1 si c_j <= m'  y  -1 si no,
#     sigma = la permutacion que ordena los cl_j,
#
#     tau^B_t(eta)  =?  sgn(sigma) prod_j eps_j.
#
# Es el caracter SIGNO del grupo de Weyl leido en el patron de residuos, que es lo que el plegado
# afin predice.  Se comprueba contra el caracter calculado por Freudenthal, que no sabe nada de esto.
#
# CONTROLES
#   C0  se comparan VALORES, no solo anulaciones.
#   C1  SEÑUELO: la misma formula sin el prod eps (solo sgn sigma).  Tiene que fallar.
#   C2  n impreso siempre; y se cuenta aparte cuantos sobreviven.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage odd_sign_formula.sage

import json

_CH = {}
def car(rk, mu):
    key = (rk, tuple(int(v) for v in mu))
    if key not in _CH:
        W = WeylCharacterRing("B%d" % rk)
        el = W(W.space().from_vector(vector([Integer(v) for v in mu])))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + int(mult)
        _CH[key] = d
    return _CH[key]


def tau_B(rk, eta, t):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car(rk, eta).items():
        s += mult * z ** (sum((i + 1) * wt[i] for i in range(rk)) % t)
    return QQ(s) if s in QQ else None


def formula(eta, t, mp, con_eps=True):
    A = [2 * int(eta[j]) + 2 * (mp - j - 1) + 1 for j in range(mp)]
    c = [x % t for x in A]
    if 0 in c:
        return 0
    cl, sg = [], 1
    for x in c:
        if x <= mp:
            cl.append(x)
        else:
            cl.append(t - x)
            if con_eps:
                sg *= -1
    if len(set(cl)) != mp:
        return 0
    perm = [mp - cl[j] for j in range(mp)]
    inv = sum(1 for i in range(mp) for j in range(i + 1, mp) if perm[i] > perm[j])
    return sg * (-1) ** inv


print("=" * 108)
print("EL SIGNO DEL FILTRO IMPAR:  ¿tau^B = sgn(sigma) prod eps ?")
print("=" * 108)
print("   t | m' | eta probados | vivos | la formula acierta | SEÑUELO sin prod eps | signo al azar")
print("   " + "-" * 100)

RES = []
for t in (3, 5, 7):
    mp = (t - 1) // 2
    n = ok = viv = sen = 0
    import random as _r
    _r.seed(int(11))
    azar = 0
    for k in range(0, 2 * t + 1):
        for e in Partitions(k, max_length=mp):
            eta = tuple(list(e) + [0] * (mp - len(e)))
            v = tau_B(mp, eta, t)
            if v is None:
                continue
            n += 1
            f = formula(eta, t, mp, True)
            fs = formula(eta, t, mp, False)
            if f == v:
                ok += 1
            if fs == v:
                sen += 1
            if v != 0:
                viv += 1
                if v == (1 if _r.randint(0, 1) else -1):
                    azar += 1
    print("   %2d | %2d | %12d | %5d | %18s | %20s | %d de %d"
          % (t, mp, n, viv, "%d de %d" % (ok, n), "%d de %d" % (sen, n), azar, viv))
    RES.append({"t": int(t), "mp": int(mp), "n": int(n), "vivos": int(viv),
                "formula_ok": int(ok), "senuelo_ok": int(sen), "azar": int(azar)})

print("")
print("=" * 108)
print("  LECTURA, escrita ANTES de correr: si la formula acierta el 100 %, el signo impar queda")
print("  cerrado y la Figura 2 se puede colorear por el; si no, el panel impar se colorea solo por")
print("  supervivencia y se dice en el pie.")
json.dump(RES, open("odd_sign_formula_DUMP.json", "w"), indent=1)
print("=" * 108)
print("DONE")
