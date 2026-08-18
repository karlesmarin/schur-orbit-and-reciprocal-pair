# -*- coding: utf-8 -*-
# EL TOP TRANSVERSAL CONTRA LA LEY DE mu_max.   16 de agosto de 2026.
#
# DE DONDE SALE.  Quedan dos enunciados de "el top se lee de la entrada", en niveles distintos:
#
#   prop:transversal(v)  el top de nu(Lambda,.) lo da S_min, el transversal de v mas pequeno.
#                        Es POR Lambda, y esta probado.
#   eq:mumax             el top del compuesto Phi_{t,r}(beta) se lee de beta.  Es por BETA.
#
# Entre los dos esta (L2): "en mu_max exactamente un Lambda tiene c != 0".  Si el top de cada
# c(Lambda,.) se lee sin buscar, entonces
#
#        mu_max(beta)  =  max_{Lambda : a_Lambda != 0}  top c(Lambda, .)
#
# y (L2) deja de ser una busqueda para ser una COMPARACION de vectores explicitos.  Eso es lo que
# se mide.
#
# LO QUE SE MIDE, por forma beta
#   M1  FATAL: el mu_max medido (pelando Phi en D_r) coincide con el maximo de los top c(Lambda,.)
#       predichos por el transversal.
#   M2  ese maximo es UNICO entre los Lambda con a_Lambda != 0.   <-- la forma de (L2)
#   M3  el Lambda que alcanza el maximo es el mismo que el unico con c(Lambda, mu_max) != 0.
#       <-- (L2) propiamente dicha, cruzada con el transversal
#   M4  top c(Lambda,.) == enderezar( top nu(Lambda,.) - t.1 ), o sea el primer paso de la division.
#
# SENUELOS
#   D1  "el Lambda maximo en dominancia es el que aporta" -- ya muerta en odd_extremal.sage.
#       Si acertara, el gate no distingue nada.
#   D2  predecir con el transversal S_max (indice menor de cada clase) en vez de S_min.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage mumax_vs_transversal.sage

import json
import sys
import itertools
from collections import defaultdict

CASOS = [(3, 2, 9), (5, 2, 10), (3, 3, 9)]


def plegar(v, t):
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def sgn_perm(perm):
    n, s, visto = len(perm), 1, [False] * len(perm)
    for i in range(n):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            s = -s
    return int(s)


def enderezar_D(x):
    r = len(x)
    a = [abs(int(v)) for v in x]
    if len(set(a)) != r:
        return None
    idx = sorted(range(r), key=lambda i: -a[i])
    s = sgn_perm(list(idx))
    y = [int(x[i]) for i in idx]
    neg = sum(1 for v in y if v < 0)
    cero = any(v == 0 for v in y)
    y = [abs(v) for v in y]
    if neg % 2 == 1 and not cero:
        y[-1] = -y[-1]
    return (tuple(y), s)


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


def clases(v, t):
    d = {}
    for i, x in enumerate(v):
        d.setdefault(plegar(x, t)[0], []).append(i)
    return d


def top_transversal(Lam, t, r, cual="min"):
    """top de nu(Lambda,.) por prop:transversal(v), sin recorrer W^1.  Devuelve x = 2(mu+rho_D)."""
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    v = [2 * int(Lam[i]) + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = clases(v, t)
    if any(not d.get(j) for j in range(1, mp + 1)):
        return None
    S = frozenset((max(d[j]) if cual == "min" else min(d[j])) for j in range(1, mp + 1))
    if len(S) != mp:
        return None
    Sc = [i for i in range(Rp) if i not in S]
    return tuple(sorted([v[i] for i in Sc], reverse=True))


def domina(a, b):
    return all(sum(a[:k + 1]) >= sum(b[:k + 1]) for k in range(len(b)))


print("=" * 112)
print("EL TOP TRANSVERSAL CONTRA LA LEY DE mu_max")
print("=" * 112)
sys.stdout.flush()

RES = []
for (t, r, tope) in CASOS:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    N = t + 2 * r
    rhoD2 = [2 * (r - j) - 2 for j in range(r)]
    n = m1 = m2 = m3 = m4 = d1 = d2 = 0
    coinc = 0   # formas donde M1 falla Y M2 falla, o M1 acierta Y M2 acierta
    fallo = None
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
        if rest2:
            continue
        n += 1

        # top c(Lambda,.) predicho: enderezar( top nu - t.1 )
        pred = {}
        for Lam in aB:
            x = top_transversal(Lam, t, r)
            if x is None:
                continue
            e = enderezar_D(tuple(x[j] - t for j in range(r)))
            if e is None:
                continue
            mu = tuple((e[0][j] - rhoD2[j]) // 2 for j in range(r))
            pred[Lam] = tuple(list(mu[:-1]) + [abs(mu[-1])])

        # M1 / M2
        cL = {}
        if pred:
            alturas = sorted(((sum(v), v, L) for L, v in pred.items()), reverse=True)
            top_pred = alturas[0][1]
            m1 += 1 if top_pred == mm else 0
            unico = (len(alturas) == 1 or alturas[0][:2] != alturas[1][:2])
            m2 += 1 if unico else 0
            coinc += 1 if ((top_pred == mm) == unico) else 0
            if top_pred != mm and fallo is None:
                fallo = {"beta": [int(x) for x in b], "mu_max": [int(x) for x in mm],
                         "pred": [int(x) for x in top_pred],
                         "todos": {str(list(L)): list(v) for L, v in list(pred.items())[:5]}}
            # M3: el argmax contra el unico Lambda con c != 0 en mu_max
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
            m3 += 1 if (len(cL) == 1 and list(cL)[0] == alturas[0][2]) else 0
            # D2: el transversal S_max
            pred2 = {}
            for Lam in aB:
                x = top_transversal(Lam, t, r, cual="max")
                if x is None:
                    continue
                e = enderezar_D(tuple(x[j] - t for j in range(r)))
                if e is None:
                    continue
                mu = tuple((e[0][j] - rhoD2[j]) // 2 for j in range(r))
                pred2[Lam] = tuple(list(mu[:-1]) + [abs(mu[-1])])
            if pred2:
                t2 = max(((sum(v), v) for v in pred2.values()))[1]
                d2 += 1 if t2 == mm else 0
        # D1 SENUELO, bien planteado: .es el Lambda maximo en dominancia el que APORTA en mu_max?
        # Esa es la hipotesis H1, muerta en odd_extremal.sage.  La version anterior de este senuelo
        # comparaba dos PREDICCIONES entre si, que no testa nada.
        Ls = list(aB)
        mx = [L for L in Ls if not any(M != L and domina(M, L) for M in Ls)]
        if len(mx) == 1 and len(cL) == 1 and mx[0] == list(cL)[0]:
            d1 += 1
    print("")
    print("  t=%d r=%d  (%d formas)" % (t, r, n))
    print("     M1  mu_max == max_Lambda top c(Lambda,.)   : %3d de %3d   <== FATAL" % (m1, n))
    print("     M2  ese maximo es unico                    : %3d de %3d   <== la forma de (L2)" % (m2, n))
    print("     M3  el argmax ES el unico con c != 0       : %3d de %3d" % (m3, n))
    print("     D1  SENUELO: el Lambda maximo en dominancia: %3d de %3d  (ya muerto)" % (d1, n))
    print("     D2  SENUELO: el transversal S_max          : %3d de %3d" % (d2, n))
    print("     ++  M1 acierta EXACTAMENTE cuando el maximo es unico : %3d de %3d" % (coinc, n))
    if fallo:
        print("     !! %s" % str(fallo)[:400])
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "n": int(n), "M1": int(m1), "M2": int(m2),
                "M3": int(m3), "D1": int(d1), "D2": int(d2), "coinc": int(coinc), "fallo": None})

json.dump(RES, open("mumax_vs_transversal_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr:")
print("   * si M1 y M3 salen limpias, (L2) deja de ser una busqueda: el Lambda que aporta en mu_max")
print("     es el que maximiza un vector que se LEE de Lambda, sin evaluar nada.")
print("   * si M1 acierta y M3 no, el top se predice pero el Lambda que aporta es otro, y entonces")
print("     (L2) no se reduce a esta comparacion.")
print("=" * 112)
print("DONE")
