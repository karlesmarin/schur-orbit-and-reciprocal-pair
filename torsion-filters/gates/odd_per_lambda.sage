# -*- coding: utf-8 -*-
# LA PRIMITIVIDAD IMPAR, LOCALIZADA POR Lambda.   16 de agosto de 2026.
#
# DE DONDE SALE.  odd_extremal.sage abrio la descomposicion en mu_max y mato mi hipotesis H1 (el
# Lambda maximo NO siempre es el que aporta).  Pero enseño algo mejor, y local:
#
#     para cada Lambda por separado,   c(Lambda, mu) := sum_eta B^odd_{Lambda;eta,mu} tau^B(eta)
#     parecia valer  0 o +-1,  y en mu_max exactamente UN Lambda daba +-1.
#
# Si eso es cierto, la conjetura impar se parte en dos piezas mucho mas manejables:
#
#     (L1)  c(Lambda, mu) en {0, +-1}  para TODO (Lambda, mu)   -- un enunciado sobre
#           branching + proyeccion de fusion, SIN el paso de Littlewood;
#     (L2)  en mu = mu_max, exactamente un Lambda tiene c != 0.
#
# (L1) es lo que se llamaria primitividad del espacio de multiplicidad de UN Lambda, y es un
# enunciado de teoria de representaciones puro.  (L2) es el que sigue siendo combinatorio.
#
# LO QUE SE MIDE
#   L1  sobre TODOS los pares (Lambda, mu) de una caja: el reparto de los valores de c.  Si aparece
#       un 2, (L1) es falsa y hay que decirlo.
#   L2  en mu_max de cada forma: cuantos Lambda tienen c != 0.
#   L3  y el control de siempre: sum_Lambda a_Lambda c(Lambda, mu_max) = A.
#
# CONTROLES
#   C0  L3 es fatal.
#   C1  el reparto completo de c, no un maximo.
#   C2  se cuenta tambien cuantos (Lambda,mu) tienen c = 0 con eta's que SI sobreviven: esos son los
#       que cancelan, y son la evidencia de que hay involucion local.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage odd_per_lambda.sage

import json
import sys
from collections import Counter, defaultdict


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
            s += mult * z ** (sum((i + 1) * wt[i] for i in range(mp)) % t)
        _TAU[key] = QQ(s) if s in QQ else None
    return _TAU[key]


def domina(a, b):
    return all(sum(a[:k + 1]) >= sum(b[:k + 1]) for k in range(len(b)))


def betas(N, tope):
    return [tuple(sorted(c, reverse=True)) + (0,) for c in Combinations(range(1, tope + 1), N - 1)]


print("=" * 118)
print("LA PRIMITIVIDAD IMPAR, LOCALIZADA POR Lambda")
print("=" * 118)

RES = []
for (t, r, tope) in [(3, 2, 9), (5, 2, 10), (3, 3, 9)]:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    N = t + 2 * r
    if r > 2 and t > 3:
        continue
    repC = Counter()
    nL_en_mumax = Counter()
    ctrl = nform = 0
    cancelan = 0
    peor = None
    for b in betas(N, tope):
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
        maxi = [m for m in S if not any(n != m and domina(n, m) for n in S)]
        if len(maxi) != 1:
            continue
        mm = maxi[0]
        P2 = phi(b, 1, Rp)
        aB, rest2 = pelar(P2, "B", Rp)
        if rest2:
            continue
        nform += 1
        # L1: c(Lambda, mu) para todos los pares
        cLm = defaultdict(lambda: 0)
        vivos_por_par = defaultdict(lambda: 0)
        for Lam in aB:
            for (eta, mu), cc in branch(Rp, mp, r, Lam).items():
                tv = tauB(mp, eta, t)
                if tv is None or tv == 0:
                    continue
                mp_plus = tuple(list(mu[:-1]) + [abs(mu[-1])])
                cLm[(Lam, mp_plus)] += int(cc) * int(tv)
                vivos_por_par[(Lam, mp_plus)] += 1
        for k, v in cLm.items():
            repC[int(v)] += 1
            if v == 0 and vivos_por_par[k] > 0:
                cancelan += 1
            if abs(v) > 2 and peor is None:
                peor = (tuple(int(x) for x in b), [int(x) for x in k[0]], [int(x) for x in k[1]], int(v))
        # L2: cuantos Lambda con c != 0 en mu_max
        nz = [L for L in aB if cLm.get((L, mm), 0) != 0]
        nL_en_mumax[len(nz)] += 1
        suma = sum(int(aB[L]) * cLm[(L, mm)] for L in aB)
        fac = 2 if mm[-1] != 0 else 1
        if QQ(suma) == fac * Ap[mm]:
            ctrl += 1
    print("")
    print("  t=%d r=%d  (%d formas con mu_max unico)" % (t, r, nform))
    print("     C0  control sum_Lambda a_L c(L,mu_max) = A : %d de %d" % (ctrl, nform))
    print("     L1  reparto de c(Lambda,mu) sobre TODOS los pares : %s" % dict(sorted(repC.items())))
    print("         (los +-2 son el factor quiral: mu con ultima coordenada no nula cuenta dos veces)")
    print("     L2  numero de Lambda con c != 0 en mu_max : %s" % dict(sorted(nL_en_mumax.items())))
    print("     C2  pares (Lambda,mu) con eta supervivientes que CANCELAN a 0 : %d" % cancelan)
    if peor:
        print("     !! un c con |c| > 2 : %s" % str(peor))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "n_formas": int(nform), "control": int(ctrl),
                "reparto_c": {str(k): int(v) for k, v in repC.items()},
                "n_Lambda_en_mumax": {str(k): int(v) for k, v in nL_en_mumax.items()},
                "pares_que_cancelan": int(cancelan)})

print("")
print("=" * 118)
print("  LECTURA, escrita ANTES de correr:")
print("   * si el reparto de c vive en {0,+-1,+-2} (con el +-2 quiral) y en mu_max hay UN SOLO")
print("     Lambda con c != 0, la conjetura impar se parte en (L1) representacion y (L2) combinatoria")
print("     -- y (L1) es demostrable con branching + fusion, sin Littlewood.")
print("   * si aparece |c| > 2, (L1) es falsa y el impar no se localiza.")
json.dump(RES, open("odd_per_lambda_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
