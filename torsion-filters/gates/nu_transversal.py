# -*- coding: utf-8 -*-
# EL NUMERADOR GKRS ES UN RECUENTO DE TRANSVERSALES.   16 de agosto de 2026.
#
# DE DONDE SALE.  Cruce de dos cosas de hoy que estaban en secciones distintas:
#
#   (a)  lem:muinj  --  w en W^1 ES el dato de un subconjunto S de tamano m' mas una quiralidad,
#        y  eta_w  es el bloque S de  v = 2 Lambda + 2 rho_{B_{R'}}.
#   (b)  cor:oddsign --  tau^B_t(eta) = eps_t . delta(A),  y delta != 0 exactamente cuando las
#        clases plegadas de A son una permutacion de {1,...,m'}.
#
# Y AHORA EL PASO QUE LO CIERRA: como  t = 2m'+1,  las clases plegadas NO NULAS son exactamente
# {1,...,m'}, ni una mas.  Luego pedir que las m' clases de v|_S sean una permutacion de {1..m'}
# es pedir que S tome UN elemento de CADA clase no nula y NINGUNO de la clase 0.  Es decir:
#
#            S admisible  <=>  S es un TRANSVERSAL de las clases plegadas no nulas de v.
#
# Consecuencias inmediatas, y son las que se miden:
#
#      |supp nu|  =  2 . prod_{j=1}^{m'} n_j ,     n_j = #{ i : v_i pliega a la clase j }
#      nu = 0     <=>  algun n_j = 0   (alguna clase no nula queda sin cubrir)
#      y el signo de cada termino es  eps_t . sgn(w) . delta(v|_S),  todo explicito.
#
# POR QUE IMPORTA.  Deja el numerador de (L1) SIN representaciones: es un recuento con signo de
# transversales de una particion de R' residuos en m'+1 clases.  Y (L1) pasa a ser un enunciado
# combinatorio sobre transversales con signo, dividido por Delta_t.
#
# LO QUE SE MIDE
#   G1  |supp nu| == 2 prod n_j                          FATAL
#   G2  supp nu == { transversales } x {+-}  PUNTO A PUNTO, no solo el cardinal
#   G3  el valor de nu en cada punto == eps_t . sgn(w) . delta(v|_S)
#   G4  nu vacio  <=>  algun n_j = 0
#
# CONTROLES
#   C0  G1 y G2 son fatales.
#   D1  SENUELO: contar las clases modulo t+2 en vez de t.  Debe fallar.
#   D2  SENUELO: prod n_j sin el factor 2 de la quiralidad.  Debe fallar.
#   C2  n impreso siempre, y el reparto de los n_j.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python nu_transversal.py

import itertools
import json
from collections import Counter


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


def eps_de(t):
    return jacobi((-2) % t, t) ** ((t + 3) // 2)


def nu_directo(Lam, t, r):
    """nu por la enumeracion honesta de W^1, sin usar la lectura de transversales."""
    mp = (t - 1) // 2
    Rp = mp + r
    v = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
    out = {}
    for perm in itertools.permutations(range(Rp)):
        s = sgn_perm(list(perm))
        base = [v[perm[i]] for i in range(Rp)]
        for eps in itertools.product((1, -1), repeat=Rp):
            u = [base[i] * eps[i] for i in range(Rp)]
            if not (all(u[i] > u[i + 1] for i in range(mp - 1)) and u[mp - 1] > 0):
                continue
            f = u[mp:]
            if not (all(f[i] > f[i + 1] for i in range(r - 1)) and f[r - 2] > abs(f[r - 1])):
                continue
            sg = s
            for e in eps:
                sg *= e
            tv = eps_de(t) * delta_de(u[:mp], t, mp)
            if tv == 0:
                continue
            out[tuple(u[mp:])] = out.get(tuple(u[mp:]), 0) + sg * tv
    return {k: v2 for k, v2 in out.items() if v2 != 0}, v


def clases(v, t):
    """particion de los indices por clase plegada de v_i mod t."""
    d = {}
    for i, x in enumerate(v):
        c_, _ = plegar(x, t)
        d.setdefault(c_, []).append(i)
    return d


def transversales(v, t, mp):
    """los S que toman un elemento de cada clase no nula 1..mp y ninguno de la clase 0."""
    d = clases(v, t)
    if any(len(d.get(j, [])) == 0 for j in range(1, mp + 1)):
        return []
    out = []
    for pick in itertools.product(*[d[j] for j in range(1, mp + 1)]):
        out.append(frozenset(pick))
    return out


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
print("EL NUMERADOR GKRS COMO RECUENTO DE TRANSVERSALES")
print("=" * 112)

RES = []
for (t, r, cota) in CASOS:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    n = g1 = g2 = g3 = g4 = d1 = d2 = 0
    rep_nj = Counter()
    fallo = None
    for Lam in dominantes(Rp, cota):
        nu, v = nu_directo(Lam, t, r)
        d = clases(v, t)
        nj = [len(d.get(j, [])) for j in range(1, mp + 1)]
        for x in nj:
            rep_nj[x] += 1
        pred = 2
        for x in nj:
            pred *= x
        n += 1
        # G1
        if len(nu) == pred:
            g1 += 1
        elif fallo is None:
            fallo = {"Lambda": list(Lam), "v": v, "n_j": nj, "pred": pred, "real": len(nu)}
        # G4
        vacio = any(x == 0 for x in nj)
        if (len(nu) == 0) == vacio:
            g4 += 1
        # G2 y G3: punto a punto
        Ts = transversales(v, t, mp)
        prev = {}
        for S in Ts:
            Sc = [i for i in range(Rp) if i not in S]
            # el bloque congelado, ordenado decreciente y positivo; el libre igual con quiralidad
            A = sorted([v[i] for i in S], reverse=True)
            for qui in (1, -1):
                libre = sorted([v[i] for i in Sc], reverse=True)
                libre[-1] *= qui
                # signo de la permutacion con signo que lleva v a (A | libre)
                orden = sorted(S, key=lambda i: -v[i]) + sorted(Sc, key=lambda i: -v[i])
                sg = sgn_perm([orden.index(i) for i in range(Rp)])
                if qui == -1:
                    sg = -sg
                prev[tuple(libre)] = sg * eps_de(t) * delta_de(A, t, mp)
        prev = {k: v2 for k, v2 in prev.items() if v2 != 0}
        if set(prev) == set(nu):
            g2 += 1
        if prev == nu:
            g3 += 1
        # senuelos
        d_mal = clases(v, t + 2)
        nj_mal = [len(d_mal.get(j, [])) for j in range(1, mp + 1)]
        p_mal = 2
        for x in nj_mal:
            p_mal *= x
        if len(nu) == p_mal:
            d1 += 1
        p_sin2 = 1
        for x in nj:
            p_sin2 *= x
        if len(nu) == p_sin2:
            d2 += 1
    print("")
    print("  t=%d r=%d (m'=%d R'=%d)  %d pesos" % (t, r, mp, Rp, n))
    print("     G1  |supp nu| == 2 prod n_j        : %3d de %3d" % (g1, n))
    print("     G2  supp nu == transversales x +-  : %3d de %3d  (punto a punto)" % (g2, n))
    print("     G3  y los VALORES tambien          : %3d de %3d" % (g3, n))
    print("     G4  nu vacio <=> algun n_j = 0     : %3d de %3d" % (g4, n))
    print("     D1  senuelo: clases mod t+2        : %3d de %3d  (debe ser bajo)" % (d1, n))
    print("     D2  senuelo: sin el factor 2       : %3d de %3d  (debe ser bajo)" % (d2, n))
    print("     C2  reparto de los n_j             : %s" % dict(sorted(rep_nj.items())))
    if fallo:
        print("     !! primer fallo de G1: %s" % json.dumps(fallo))
    RES.append({"t": t, "r": r, "n": n, "G1": g1, "G2": g2, "G3": g3, "G4": g4,
                "D1": d1, "D2": d2, "n_j": {str(k): v2 for k, v2 in rep_nj.items()},
                "fallo": fallo})

json.dump(RES, open("nu_transversal_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr:")
print("   * si G1-G3 salen limpias, el numerador de (L1) queda SIN representaciones: es un recuento")
print("     con signo de transversales de una particion de R' residuos en m'+1 clases.")
print("   * y G4 da un criterio de anulacion que no necesita calcular nada: basta una clase vacia.")
print("   * si G2 sale y G3 no, la lectura acierta el SOPORTE y falla el signo -- que es justo la")
print("     distincion que nos costo la vuelta 28.")
print("=" * 112)
print("DONE")
