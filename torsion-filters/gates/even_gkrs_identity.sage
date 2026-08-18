# -*- coding: utf-8 -*-
# E6: LA IDENTIDAD GKRS PAR, COMPLETA.   16 de agosto de 2026.
#
# even_transversal.sage declaro este test en su cabecera y NO lo implemento.  Un guion que promete
# un control y no lo corre es peor que uno que no lo promete, asi que va aparte y entero.
#
# LA CUENTA, antes de programar nada.  Para  Sp_{2R} > Sp_{2m} x Sp_{2r}  las raices positivas
# complementarias son  e_i +- f_j  (i <= m, j <= r) y NO hay  f_j  suelta.  El denominador relativo
# es  prod (e^{a/2} - e^{-a/2})  sobre ellas; congelando  y_i = xi^i  con  t = 2m+2  el alfabeto
# {xi^{+-i}} son todas las raices t-esimas MENOS +-1, luego por coordenada libre
#
#     prod_{i=1..m} (1 - xi^i z)(1 - xi^{-i} z)  =  (1 - z^t)/(1 - z^2)
#          =  monomio . (z^{t/2} - z^{-t/2})/(z - z^{-1})
#          =  monomio . ( z^{t/2-1} + z^{t/2-3} + ... + z^{-(t/2-1)} )        [t/2 terminos]
#
# Asi que en la base de ALTERNANTES de C_r la multiplicacion por el denominador NO son 2^r
# desplazamientos como en el impar: son (t/2)^r, con exponentes  t/2-1, t/2-3, ..., -(t/2-1),
# TODOS con coeficiente +1 -- no hay signos alternos.  Y el enderezado es por W(C_r), donde todos
# los cambios de signo estan permitidos (a diferencia de D_r) y un x_i = 0 es singular.
#
# LO QUE SE MIDE
#   E6  ( sum_mu c(L,mu) A_mu ) . Delta^par  ==  +- nu     en la base de alternantes de C_r.  FATAL.
#       c se calcula por BRANCHING C_R -> C_m x C_r con tau^C evaluado del caracter, y nu por
#       enumeracion de W^1 -- dos caminos que no comparten codigo salvo tau^C.
#   E6b el signo global, .es constante dentro de cada (t,r)?
#
# SENUELOS
#   D1  el denominador del IMPAR: desplazamientos {+t/2, -t/2} con signo alterno.  Debe fallar.
#   D2  el mismo numero de terminos pero con los exponentes de t+2.  Debe fallar.
#   D3  los exponentes correctos pero enderezando por W(D_r) en vez de W(C_r).  Debe fallar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage even_gkrs_identity.sage

import json
import sys
import itertools
from collections import Counter, defaultdict

CASOS = [(4, 2, 5), (6, 2, 4), (8, 2, 3), (6, 3, 3), (10, 2, 2)]


# ------------------------------------------------------------------ W(C_r): enderezar
def enderezar_C(x):
    """x en Z^r.  Devuelve (x dominante regular, det w) o None si es singular.
       En C_r todo cambio de signo esta permitido; x_i = 0 o |x_i| = |x_j| es singular."""
    r = len(x)
    a = [abs(int(v)) for v in x]
    if 0 in a or len(set(a)) != r:
        return None
    idx = sorted(range(r), key=lambda i: -a[i])
    s, visto = 1, [False] * r
    for i in range(r):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = idx[j]
            L += 1
        if L % 2 == 0:
            s = -s
    y = [int(x[i]) for i in idx]
    for v in y:
        if v < 0:
            s = -s
    return (tuple(abs(v) for v in y), int(s))


def enderezar_D(x):
    """el senuelo D3: enderezar por W(D_r), donde solo valen los cambios de signo en numero par."""
    r = len(x)
    a = [abs(int(v)) for v in x]
    if len(set(a)) != r:
        return None
    idx = sorted(range(r), key=lambda i: -a[i])
    s, visto = 1, [False] * r
    for i in range(r):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = idx[j]
            L += 1
        if L % 2 == 0:
            s = -s
    y = [int(x[i]) for i in idx]
    neg = sum(1 for v in y if v < 0)
    cero = any(v == 0 for v in y)
    y = [abs(v) for v in y]
    if neg % 2 == 1 and not cero:
        y[-1] = -y[-1]
    return (tuple(y), int(s))


def desplazar(x, exps, coefs, ender):
    """multiplica el alternante A_x por prod_j g(z_j), con g = sum coefs[k] z^{exps[k]}."""
    r = len(x)
    out = defaultdict(lambda: 0)
    for pick in itertools.product(range(len(exps)), repeat=r):
        c = 1
        for k in pick:
            c *= coefs[k]
        y = tuple(int(x[j]) + exps[pick[j]] for j in range(r))
        e = ender(y)
        if e is None:
            continue
        out[e[0]] += c * e[1]
    return {k: v for k, v in out.items() if v != 0}


# ------------------------------------------------------------------ caracteres
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


_TAU = {}
def tauC(eta, t, m):
    key = (tuple(int(v) for v in eta), t, m)
    if key not in _TAU:
        K = CyclotomicField(t)
        z = K.gen()
        s = K(0)
        for wt, mult in car("C", m, eta).items():
            s += mult * z ** (sum((i + 1) * int(wt[i]) for i in range(m)) % t)
        _TAU[key] = int(QQ(s)) if s in QQ else None
    return _TAU[key]


_BR = {}
def branch(Rr, m, r, Lam):
    key = (Rr, m, r, tuple(int(v) for v in Lam))
    if key not in _BR:
        W = WeylCharacterRing("C%d" % Rr)
        X = WeylCharacterRing("C%dxC%d" % (m, r))
        br = branching_rule("C%d" % Rr, "C%dxC%d" % (m, r), "orthogonal_sum")
        el = W(W.space().from_vector(vector([Integer(v) for v in Lam])))
        d = {}
        for wt, c in el.branch(X, rule=br).monomial_coefficients().items():
            v = [int(u) for u in wt.to_vector()]
            d[(tuple(v[:m]), tuple(v[m:]))] = int(c)
        _BR[key] = d
    return _BR[key]


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


def dominantes(rk, cota):
    def rec(k, tope):
        if k == 0:
            yield ()
            return
        for a in range(tope, -1, -1):
            for resto in rec(k - 1, a):
                yield (a,) + resto
    return list(rec(rk, cota))


print("=" * 112)
print("E6: LA IDENTIDAD GKRS PAR, COMPLETA")
print("=" * 112)
sys.stdout.flush()

RES = []
for (t, r, cota) in CASOS:
    m = (t - 2) // 2
    Rr = m + r
    if m < 1:
        continue
    h = t // 2
    EXPS = [h - 1 - 2 * k for k in range(h)]          # t/2-1, t/2-3, ..., -(t/2-1)
    COEFS = [1] * h
    rhoC = [r - j for j in range(r)]                  # rho_{C_r} = (r, ..., 1) -> +1 abajo
    n = e6 = d1 = d2 = d3 = ntriv = 0
    signos = Counter()
    fallo = None
    for Lam in dominantes(Rr, cota):
        # ---- c(Lambda, mu) por branching + filtro
        c = defaultdict(lambda: 0)
        for (eta, mu), b in branch(Rr, m, r, Lam).items():
            tv = tauC(eta, t, m)
            if not tv:
                continue
            c[tuple(int(u) for u in mu)] += int(b) * int(tv)
        c = {k: v for k, v in c.items() if v != 0}

        # ---- nu por W^1 : subconjunto S de tamano m, todo positivo (los dos bloques son tipo C)
        v = [int(Lam[i]) + Rr - i for i in range(Rr)]
        nu = defaultdict(lambda: 0)
        for S in itertools.combinations(range(Rr), m):
            Sc = [i for i in range(Rr) if i not in S]
            A = sorted([v[i] for i in S], reverse=True)
            eta = tuple(A[j] - (m - j) for j in range(m))
            if min(eta) < 0:
                continue
            tv = tauC(eta, t, m)
            if not tv:
                continue
            orden = sorted(S, key=lambda i: -v[i]) + sorted(Sc, key=lambda i: -v[i])
            sg = sgn_perm([orden.index(i) for i in range(Rr)])
            x = tuple(sorted([v[i] for i in Sc], reverse=True))
            nu[x] += sg * int(tv)
        nu = {k: val for k, val in nu.items() if val != 0}

        n += 1
        if not c and not nu:
            ntriv += 1

        def lado(exps, coefs, ender):
            out = defaultdict(lambda: 0)
            for mu, cv in c.items():
                x = tuple(int(mu[j]) + rhoC[j] for j in range(r))   # x = mu + rho_{C_r} = mu + (r,...,1)
                for y, sg in desplazar(x, exps, coefs, ender).items():
                    out[y] += cv * sg
            return {k: val for k, val in out.items() if val != 0}

        lhs = lado(EXPS, COEFS, enderezar_C)
        neg = {k: -val for k, val in lhs.items()}
        ok = (lhs == nu) or (neg == nu)
        e6 += 1 if ok else 0
        if lhs and lhs == nu:
            signos[+1] += 1
        elif lhs and neg == nu:
            signos[-1] += 1
        if not ok and fallo is None:
            fallo = {"Lambda": [int(a) for a in Lam],
                     "c": {str(k): int(val) for k, val in c.items()},
                     "nu": {str(k): int(val) for k, val in nu.items()},
                     "lhs": {str(k): int(val) for k, val in lhs.items()}}

        if c:
            # D1: el denominador del impar (dos terminos con signo alterno)
            l1 = lado([h, -h], [1, -1], enderezar_C)
            d1 += 1 if (l1 == nu or {k: -val for k, val in l1.items()} == nu) else 0
            # D2: los exponentes de t+2
            h2 = (t + 2) // 2
            l2 = lado([h2 - 1 - 2 * k for k in range(h2)], [1] * h2, enderezar_C)
            d2 += 1 if (l2 == nu or {k: -val for k, val in l2.items()} == nu) else 0
            # D3: exponentes correctos, enderezado por W(D_r)
            l3 = lado(EXPS, COEFS, enderezar_D)
            d3 += 1 if (l3 == nu or {k: -val for k, val in l3.items()} == nu) else 0
    nt = n - ntriv
    print("")
    print("  t=%d r=%d (m=%d R=%d)  %d pesos  (%d con c y nu vacios)" % (t, r, m, Rr, n, ntriv))
    print("     E6  ( sum c A_mu ) . Delta^par == +- nu   : %3d de %3d   <== FATAL" % (e6, n))
    print("     E6b signo global                          : %s" % dict(sorted(signos.items())))
    print("     D1  SENUELO denominador del impar         : %3d de %3d  (debe ser 0)" % (d1, nt))
    print("     D2  SENUELO exponentes de t+2             : %3d de %3d  (debe ser 0)" % (d2, nt))
    print("     D3  SENUELO enderezar por W(D_r)          : %3d de %3d  (debe ser 0)" % (d3, nt))
    if fallo:
        print("     !! primer fallo: %s" % json.dumps(fallo)[:500])
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "m": int(m), "R": int(Rr), "n": int(n),
                "no_triviales": int(nt), "E6": int(e6),
                "signo": {str(k): int(val) for k, val in signos.items()},
                "D1": int(d1), "D2": int(d2), "D3": int(d3), "fallo": fallo})

json.dump(RES, open("even_gkrs_identity_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr:")
print("   * si E6 acierta y los tres senuelos no, la identidad GKRS especializada vale TAMBIEN en el")
print("     par, con su propio denominador (1-z^t)/(1-z^2) -- y entonces la unica asimetria que")
print("     queda entre las paridades es la virtualidad de los a_Lambda, no el aparato.")
print("   * si E6 falla, el denominador par no es el que la cuenta de raices predice, y hay que")
print("     mirar que raiz complementaria se ha contado mal.")
print("=" * 112)
print("DONE")
