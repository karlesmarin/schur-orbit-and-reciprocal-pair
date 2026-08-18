# -*- coding: utf-8 -*-
# .ES EL TRANSVERSAL DE mu_max EL DE COORDENADAS MAS PEQUENAS?   16 de agosto de 2026.
#
# DE DONDE SALE.  prop:transversal deja el numerador nu indexado por transversales S de las clases
# plegadas de  v = 2 Lambda + 2 rho_{B_{R'}},  y el indice alternante del lado libre es  v|_{S^c}.
# Como v es estrictamente decreciente, "coordenadas mas pequenas" = "indices mas grandes".  Y
# elegir en cada clase el indice MAS GRANDE minimiza v|_S coordenada a coordenada, luego deja en el
# complemento los valores mas grandes.  La pregunta es si ESE es el que da el peso maximo.
#
# Es la pregunta que cruza prop:transversal con (L2) -- "en mu_max hay exactamente un Lambda con
# c != 0" -- porque si el maximo lo da un transversal CANONICO, el peso maximo se lee de Lambda sin
# recorrer W^1.
#
# LO QUE SE MIDE, por Lambda
#   M1  el maximo de supp(nu) (por (suma, lex)) .lo da  S_min = {indice mayor de cada clase}
#       con quiralidad +1?
#   M2  .es unico ese maximo, o hay empates?
#   M3  el maximo de supp(c) .es  enderezar(max supp(nu) - t.1) ?   (el primer paso de la division)
#   M4  y .es  S_min  tambien el que minimiza  sum(v|_S)  sobre todos los transversales?  (esto es
#       automatico coordenada a coordenada; se comprueba porque si fallara, el enunciado esta mal)
#
# CONTROLES
#   C0  M1 es fatal para la hipotesis.
#   D1  SENUELO: S_max = {indice MENOR de cada clase}.  Si tambien diera el maximo, la pregunta no
#       distingue nada.
#   D2  SENUELO: quiralidad -1 en S_min.
#   C2  n impreso siempre; y se guarda el primer contraejemplo entero.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python mumax_transversal.py

import itertools
import json


def plegar(v, t):
    v %= t
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
    return s


def jacobi(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def delta_de(A, t, mp):
    cl, ep = [], []
    for v in A:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, mp + 1)):
        return 0
    s = sgn_perm([c - 1 for c in cl])
    for e in ep:
        s *= e
    return int(s)


def enderezar(x):
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
    return (tuple(y), s)


def desplazar(x, paso):
    r = len(x)
    out = {}
    for eps in itertools.product((1, -1), repeat=r):
        sg = 1
        for e in eps:
            sg *= e
        e2 = enderezar(tuple(int(x[j]) + paso * eps[j] for j in range(r)))
        if e2 is None:
            continue
        out[e2[0]] = out.get(e2[0], 0) + sg * e2[1]
    return {k: v for k, v in out.items() if v != 0}


def cabeza(d):
    return max(d, key=lambda k: (sum(k), k))


def dividir(nu, t, r, tope=20000):
    P = dict(nu)
    c = {}
    for _ in range(tope):
        P = {k: v for k, v in P.items() if v != 0}
        if not P:
            return c, {}
        y = cabeza(P)
        cand = None
        for eps in itertools.product((1, -1), repeat=r):
            e = enderezar(tuple(int(y[j]) - t * eps[j] for j in range(r)))
            if e is None:
                continue
            D = desplazar(e[0], t)
            if D and cabeza(D) == y:
                cand = (e[0], D)
                break
        if cand is None:
            return c, P
        x, D = cand
        if P[y] % D[y] != 0:
            return c, P
        cv = P[y] // D[y]
        c[x] = c.get(x, 0) + cv
        for k, v in D.items():
            nv = P.get(k, 0) - cv * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return c, P


def clases(v, t):
    d = {}
    for i, x in enumerate(v):
        c_, _ = plegar(x, t)
        d.setdefault(c_, []).append(i)
    return d


def x_de(S, v, Rp, qui):
    """el indice alternante del bloque libre para el transversal S y la quiralidad qui."""
    Sc = [i for i in range(Rp) if i not in S]
    libre = sorted([v[i] for i in Sc], reverse=True)
    libre[-1] *= qui
    return tuple(libre)


def nu_de(Lam, t, r):
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    v = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = clases(v, t)
    if any(len(d.get(j, [])) == 0 for j in range(1, mp + 1)):
        return {}, v, d
    out = {}
    eps_t = jacobi((-2) % t, t) ** ((t + 3) // 2)
    for pick in itertools.product(*[d[j] for j in range(1, mp + 1)]):
        S = frozenset(pick)
        if len(S) != mp:
            continue
        A = sorted([v[i] for i in S], reverse=True)
        Sc = [i for i in range(Rp) if i not in S]
        for qui in (1, -1):
            orden = sorted(S, key=lambda i: -v[i]) + sorted(Sc, key=lambda i: -v[i])
            sg = sgn_perm([orden.index(i) for i in range(Rp)])
            if qui == -1:
                sg = -sg
            val = sg * eps_t * delta_de(A, t, mp)
            if val:
                out[x_de(S, v, Rp, qui)] = val
    return out, v, d


def dominantes(rk, cota):
    def rec(k, tope):
        if k == 0:
            yield ()
            return
        for a in range(tope, -1, -1):
            for resto in rec(k - 1, a):
                yield (a,) + resto
    return list(rec(rk, cota))


CASOS = [(3, 2, 6), (5, 2, 4), (7, 2, 3), (3, 3, 4), (5, 3, 2)]

print("=" * 112)
print(".ES EL TRANSVERSAL DE mu_max EL DE COORDENADAS MAS PEQUENAS?")
print("=" * 112)

RES = []
for (t, r, cota) in CASOS:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    n = m1 = m2 = m3 = m4 = d1 = d2 = 0
    contra = None
    for Lam in dominantes(Rp, cota):
        nu, v, d = nu_de(Lam, t, r)
        if not nu:
            continue
        n += 1
        top = cabeza(nu)
        # S_min: en cada clase no nula, el indice MAS GRANDE (valor de v mas pequeno)
        S_min = frozenset(max(d[j]) for j in range(1, mp + 1))
        S_max = frozenset(min(d[j]) for j in range(1, mp + 1))
        ok_min = (len(S_min) == mp) and (x_de(S_min, v, Rp, 1) == top)
        m1 += 1 if ok_min else 0
        if not ok_min and contra is None:
            contra = {"Lambda": list(Lam), "v": v,
                      "clases": {str(k): val for k, val in d.items()},
                      "top_real": list(top),
                      "top_S_min": list(x_de(S_min, v, Rp, 1)) if len(S_min) == mp else None}
        # M2 unicidad del maximo
        alturas = sorted(((sum(k), k) for k in nu), reverse=True)
        m2 += 1 if len(alturas) == 1 or alturas[0] != alturas[1] else 0
        # M3 el primer paso de la division
        c, resto = dividir(nu, t, r)
        if c and not resto:
            e = enderezar(tuple(top[j] - t for j in range(r)))
            m3 += 1 if (e is not None and cabeza(c) == e[0]) else 0
        # M4 S_min minimiza la suma
        sumas = []
        for pick in itertools.product(*[d[j] for j in range(1, mp + 1)]):
            S = frozenset(pick)
            if len(S) == mp:
                sumas.append((sum(v[i] for i in S), S))
        if sumas:
            m4 += 1 if min(sumas)[1] == S_min else 0
        # senuelos
        if len(S_max) == mp and x_de(S_max, v, Rp, 1) == top:
            d1 += 1
        if len(S_min) == mp and x_de(S_min, v, Rp, -1) == top:
            d2 += 1
    print("")
    print("  t=%d r=%d (m'=%d R'=%d)  %d formas con nu != 0" % (t, r, mp, Rp, n))
    print("     M1  top(nu) lo da S_min con quiralidad +1 : %3d de %3d   <== la hipotesis" % (m1, n))
    print("     M2  ese maximo es unico (sin empate)      : %3d de %3d" % (m2, n))
    print("     M3  top(c) == enderezar(top(nu) - t.1)    : %3d de %3d" % (m3, n))
    print("     M4  S_min minimiza sum(v|_S)              : %3d de %3d" % (m4, n))
    print("     D1  SENUELO S_max (indice menor)          : %3d de %3d  (debe ser bajo)" % (d1, n))
    print("     D2  SENUELO quiralidad -1                 : %3d de %3d  (debe ser 0)" % (d2, n))
    if contra:
        print("     !! primer contraejemplo de M1: %s" % json.dumps(contra))
    RES.append({"t": t, "r": r, "n": n, "M1": m1, "M2": m2, "M3": m3, "M4": m4,
                "D1": d1, "D2": d2, "contraejemplo": contra})

json.dump(RES, open("mumax_transversal_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr:")
print("   * si M1 sale limpia, el peso maximo del numerador se lee de Lambda SIN recorrer W^1:")
print("     basta tomar en cada clase plegada el indice mayor.  Eso es la mitad de (L2) hecha")
print("     computo, no busqueda.")
print("   * si M1 falla, el contraejemplo dice que clase rompe el orden, y ahi esta el contenido.")
print("   * si el senuelo D1 tambien acierta, la pregunta no distingue y hay que afinarla.")
print("=" * 112)
print("DONE")
