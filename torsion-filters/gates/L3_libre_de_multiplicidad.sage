# -*- coding: utf-8 -*-
# .ES LIBRE DE MULTIPLICIDAD LA RESTRICCION IMPAR?   16 de agosto de 2026.
#
# L3_multiplicidad.sage encontro que (L3) es VACIA en la poblacion donde se enuncio: los 476
# coeficientes a^B_Lambda valen 1, todos, y por tanto que el Lambda que aporta valga 1 no distingue
# nada.  Los dos senuelos aciertan 132 de 132.
#
# Asi que la pregunta buena no es (L3), es esta:
#
#      .es  Phi_{1,R'} = sum_Lambda a^B_Lambda o_Lambda  libre de multiplicidad SIEMPRE?
#
# Aqui se empuja el rango todo lo que aguante, porque un enunciado vacio en una caja pequena puede
# dejar de serlo una talla arriba -- y si se rompe, (L3) recupera contenido exactamente ahi.
#
# No hace falta branching ni filtro: basta pelar Phi en B_{R'}.  Es mucho mas barato que el gate
# anterior, asi que el rango puede ser bastante mayor.
#
# CONTROLES
#   C0  el pelado tiene que cerrar (resto vacio); si no, el dato no vale.
#   C1  se cuenta tambien cuantos Lambda hay por forma, para que "todos valen 1" no sea "hay uno".
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage L3_libre_de_multiplicidad.sage

import json
import sys
from collections import Counter

CASOS = [(3, 2, 13)]   # el caso portante, aislado para que la corrida TERMINE y vuelque


def phi(beta, tt, nvar):
    Nn = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    Lr = LaurentPolynomialRing(K, nvar, 'z')
    xx = [Lr(K(zeta) ** k) for k in range(tt)] + [g ** e for g in Lr.gens() for e in (1, -1)]
    delta = list(range(Nn - 1, -1, -1))

    def alt(expo):
        return matrix(Lr, Nn, Nn, lambda i, j: xx[i] ** expo[j]).determinant()

    den = alt(delta)
    if den == 0:
        return None
    q = alt(list(beta)) / den
    try:
        q = Lr(q)
    except Exception:
        return None
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        k = tuple(int(v) for v in e) if hasattr(e, '__iter__') else (int(e),)
        if c not in QQ:
            return None
        if QQ(c) != 0:
            out[k] = QQ(c)
    return out


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


def pelar_B(P, rk, tope=20000):
    P = {e: QQ(c) for e, c in P.items() if c != 0}
    out = {}
    for _ in range(tope):
        if not P:
            return out, {}
        dom = [e for e in P if list(e) == sorted(e, reverse=True) and min(e) >= 0]
        if not dom:
            return out, P
        mu = max(dom, key=lambda e: (sum(e), e))
        c = P[mu]
        out[mu] = out.get(mu, 0) + c
        for k, v in car("B", rk, mu).items():
            nv = P.get(k, 0) - c * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P


print("=" * 108)
print(".ES LIBRE DE MULTIPLICIDAD LA RESTRICCION IMPAR  Phi_{1,R'} -> B_{R'}?")
print("=" * 108)
sys.stdout.flush()

RES = []
for (t, r, tope) in CASOS:
    Rp = (t - 1) // 2 + r
    N = t + 2 * r
    n = cierra = 0
    coef = Counter()
    nlam = Counter()
    peores = []
    for b in [tuple(sorted(c, reverse=True)) + (0,) for c in Combinations(range(1, tope + 1), N - 1)]:
        P = phi(b, 1, Rp)
        if not P:
            continue
        aB, rest = pelar_B(P, Rp)
        n += 1
        if rest:
            continue
        cierra += 1
        aB = {k: v for k, v in aB.items() if v != 0}
        nlam[len(aB)] += 1
        for L, a in aB.items():
            coef[int(a)] += 1
            if abs(int(a)) != 1 and len(peores) < 4:
                peores.append({"beta": [int(x) for x in b], "Lambda": [int(x) for x in L],
                               "a": int(a)})
    tot = sum(coef.values())
    unos = coef.get(1, 0)
    print("")
    print("  t=%d r=%d  beta_i <= %d   (%d formas, %d con el pelado cerrado)" % (t, r, tope, n, cierra))
    print("     reparto de a^B_Lambda : %s" % dict(sorted(coef.items())))
    print("     |a| = 1 en           : %d de %d = %.2f%%" % (unos, tot, 100.0 * unos / tot if tot else 0))
    print("     numero de Lambda por forma : %s" % dict(sorted(nlam.items())))
    if peores:
        print("     !! coeficientes distintos de 1: %s" % json.dumps(peores[:3]))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "tope": int(tope), "n": int(n), "cierra": int(cierra),
                "coef": {str(k): int(v) for k, v in coef.items()},
                "n_lambda": {str(k): int(v) for k, v in nlam.items()},
                "peores": peores})

json.dump(RES, open("L3_libre_de_multiplicidad_DUMP.json", "w"), indent=1)
print("")
print("=" * 108)
print("DONE")
