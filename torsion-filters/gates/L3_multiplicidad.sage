# -*- coding: utf-8 -*-
# (L3): EL Lambda QUE APORTA TIENE MULTIPLICIDAD UNO.   16 de agosto de 2026.
#
# (L3) dice que el unico Lambda con c(Lambda, mu_max) != 0 aparece en la restriccion de Littlewood
# con  a^B_Lambda = 1.  Es la unica de las tres patas que no se ha tocado.
#
# LA PREGUNTA QUE VA PRIMERO, y decide si (L3) dice algo:
#     .cuantos de TODOS los a^B_Lambda valen 1?
# Si casi todos valen 1, que el que aporta valga 1 no es informacion: es el caso general.  Eso es un
# control de vacuidad y va antes que la medida, no despues.
#
# LO QUE SE MIDE, por forma beta
#   C0  VACUIDAD: el reparto de a^B_Lambda sobre TODOS los Lambda, y la fraccion que vale 1.
#   L3  a^B_Lambda = 1 para el Lambda que aporta en mu_max.
#   L3b .es ese Lambda el maximo en dominancia?  (H1 esta muerta; se mide cuanto)
#   X1  CRUCE con prop:transversal: .hay relacion entre a_Lambda y el numero de transversales
#       prod_j n_j de ese Lambda?  Se tabula el par (a_Lambda, prod n_j).
#   X2  y con el top: .el Lambda que aporta es el de menos transversales?
#
# CONTROLES
#   D1  SENUELO: coger el Lambda de mayor a_Lambda en vez del que aporta.  Si (L3) fuera vacia,
#       este señuelo acertaria tanto como ella.
#   D2  SENUELO: un Lambda al azar del soporte (se toma el primero en orden lexicografico).
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage L3_multiplicidad.sage

import json
import sys
import itertools
from collections import Counter, defaultdict

CASOS = [(3, 2, 9), (5, 2, 10), (3, 3, 9)]


def plegar(v, t):
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


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


def dominante(e, typ):
    f = list(e)
    if typ in ("B", "C"):
        return f == sorted(f, reverse=True) and min(f) >= 0
    if len(f) < 2:
        return f[0] >= 0
    return all(f[i] >= f[i + 1] for i in range(len(f) - 2)) and f[-2] >= abs(f[-1])


def pelar(P, typ, rk, tope=9000):
    P = {e: QQ(c) for e, c in P.items() if c != 0}
    out = {}
    for _ in range(tope):
        if not P:
            return out, {}
        dom = [e for e in P if dominante(e, typ)]
        if not dom:
            return out, P
        mu = max(dom, key=lambda e: (sum(abs(v) for v in e), e))
        c = P[mu]
        out[mu] = out.get(mu, 0) + c
        for k, v in car(typ, rk, mu).items():
            nv = P.get(k, 0) - c * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P


_BR = {}
def branch(Rp, mp, rr, Lam):
    key = (Rp, mp, rr, tuple(int(v) for v in Lam))
    if key not in _BR:
        W = WeylCharacterRing("B%d" % Rp)
        X = WeylCharacterRing("B%dxD%d" % (mp, rr))
        br = branching_rule("B%d" % Rp, "B%dxD%d" % (mp, rr), "orthogonal_sum")
        el = W(W.space().from_vector(vector([Integer(v) for v in Lam])))
        d = {}
        for wt, c in el.branch(X, rule=br).monomial_coefficients().items():
            v = [int(u) for u in wt.to_vector()]
            d[(tuple(v[:mp]), tuple(v[mp:]))] = int(c)
        _BR[key] = d
    return _BR[key]


_TAU = {}
def tauB(mp, eta, t):
    key = (mp, tuple(int(v) for v in eta), t)
    if key not in _TAU:
        K = CyclotomicField(t)
        z = K.gen()
        s = K(0)
        for wt, mult in car("B", mp, eta).items():
            s += mult * z ** (sum((i + 1) * int(wt[i]) for i in range(mp)) % t)
        _TAU[key] = int(QQ(s)) if s in QQ else None
    return _TAU[key]


def domina(a, b):
    return all(sum(a[:k + 1]) >= sum(b[:k + 1]) for k in range(len(b)))


def n_transversales(Lam, t, r):
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    v = [2 * int(Lam[i]) + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = {}
    for i, x in enumerate(v):
        d.setdefault(plegar(x, t)[0], []).append(i)
    p = 1
    for j in range(1, mp + 1):
        p *= len(d.get(j, []))
    return p


print("=" * 112)
print("(L3): EL Lambda QUE APORTA TIENE MULTIPLICIDAD UNO")
print("=" * 112)
sys.stdout.flush()

RES = []
for (t, r, tope) in CASOS:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    N = t + 2 * r
    n = l3 = l3b = d1 = d2 = 0
    todos = Counter()
    del_que_aporta = Counter()
    cruce = Counter()
    x2 = 0
    for b in [tuple(sorted(c, reverse=True)) + (0,) for c in Combinations(range(1, tope + 1), N - 1)]:
        P1 = phi(b, t, r)
        if not P1:
            continue
        A1, rest = pelar(P1, "D", r)
        A1 = {k: QQ(v) for k, v in A1.items() if v != 0}
        if rest or not A1:
            continue
        Ap = {}
        for mu, c in A1.items():
            Ap[tuple(list(mu[:-1]) + [abs(mu[-1])])] = c
        S = list(Ap)
        maxi = [m for m in S if not any(nn != m and domina(nn, m) for nn in S)]
        if len(maxi) != 1:
            continue
        mm = maxi[0]
        P2 = phi(b, 1, Rp)
        aB, rest2 = pelar(P2, "B", Rp)
        if rest2 or not aB:
            continue
        n += 1
        for L, a in aB.items():
            todos[int(a)] += 1
        # el Lambda que aporta
        cL = {}
        for Lam in aB:
            s = 0
            for (eta, mu), cc in branch(Rp, mp, r, Lam).items():
                mup = tuple(list(mu[:-1]) + [abs(mu[-1])])
                if mup != mm:
                    continue
                tv = tauB(mp, eta, t)
                if tv:
                    s += int(cc) * int(tv)
            if s:
                cL[Lam] = s
        if len(cL) != 1:
            continue
        Lstar = list(cL)[0]
        a_star = int(aB[Lstar])
        del_que_aporta[a_star] += 1
        l3 += 1 if a_star == 1 else 0
        Ls = list(aB)
        mx = [L for L in Ls if not any(M != L and domina(M, L) for M in Ls)]
        l3b += 1 if (len(mx) == 1 and mx[0] == Lstar) else 0
        cruce[(a_star, n_transversales(Lstar, t, r))] += 1
        # X2: .es el de menos transversales?
        nts = {L: n_transversales(L, t, r) for L in aB}
        mn = min(v for v in nts.values() if v > 0) if any(v > 0 for v in nts.values()) else None
        if mn is not None and nts.get(Lstar) == mn:
            x2 += 1
        # senuelos
        d1 += 1 if int(aB[max(aB, key=lambda L: int(aB[L]))]) == 1 else 0
        d2 += 1 if int(aB[sorted(aB)[0]]) == 1 else 0
    tot = sum(todos.values())
    unos = todos.get(1, 0)
    print("")
    print("  t=%d r=%d  (%d formas)" % (t, r, n))
    print("     C0  VACUIDAD: reparto de a^B_Lambda sobre TODOS : %s" % dict(sorted(todos.items())))
    print("         fraccion que vale 1                         : %d de %d = %.1f%%"
          % (unos, tot, 100.0 * unos / tot if tot else 0))
    print("     L3  a^B = 1 para el Lambda que aporta           : %3d de %3d" % (l3, sum(del_que_aporta.values())))
    print("         reparto de a^B del que aporta               : %s" % dict(sorted(del_que_aporta.items())))
    print("     L3b .es ese Lambda el maximo en dominancia?     : %3d de %3d" % (l3b, sum(del_que_aporta.values())))
    print("     X2  .es el de MENOS transversales?              : %3d de %3d" % (x2, sum(del_que_aporta.values())))
    print("     X1  cruce (a^B, num. transversales)             : %s" % dict(sorted(cruce.items())))
    print("     D1  SENUELO el Lambda de mayor a^B tiene a^B=1  : %3d de %3d" % (d1, n))
    print("     D2  SENUELO el primero en orden lex tiene a^B=1 : %3d de %3d" % (d2, n))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "n": int(n),
                "reparto_todos": {str(k): int(v) for k, v in todos.items()},
                "reparto_aporta": {str(k): int(v) for k, v in del_que_aporta.items()},
                "L3": int(l3), "L3b": int(l3b), "X2": int(x2),
                "cruce": {str(k): int(v) for k, v in cruce.items()},
                "D1": int(d1), "D2": int(d2)})

json.dump(RES, open("L3_multiplicidad_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr:")
print("   * si casi todos los a^B valen 1, (L3) es casi vacia y hay que decirlo asi, no presentarla")
print("     como una de las tres patas.")
print("   * si hay una fraccion apreciable de a^B > 1 y el que aporta vale 1 siempre, (L3) tiene")
print("     contenido y el cruce X1/X2 dice de que tipo.")
print("=" * 112)
print("DONE")
