# -*- coding: utf-8 -*-
# ¿LA CANCELACION OCURRE DENTRO DE CADA ORDER TYPE O ENTRE ELLOS?   15 de agosto de 2026.
#
# DE DONDE SALE.  Su sospecha de la vuelta 18: el emparejamiento, si existe, no vive en Lambda
# crudo -- lo mato el argumento de paridad con 10 supervivientes -- sino en el ORDER TYPE de
# Yacobi, la celda que dice que rama toma el min y que rama el max en cada posicion de la cadena.
# Dentro de un order type fijo los min/max dejan de serlo, cada r_i se vuelve una FORMA LINEAL en
# (Lambda,mu), y "r_i par" pasa a ser un sistema lineal modulo 2.  Es decir: dentro de una celda hay
# mucha mas estructura que fuera.
#
# LA PREGUNTA, tal como el la formula, y es binaria:
#
#     ¿la cancelacion a +-1 ocurre DENTRO de cada order type, o ENTRE order types?
#
#   * DENTRO  -> cada celda suma 0 salvo una, que suma +-1.  La accion local probablemente existe y
#                la involucion hay que buscarla dentro de la celda.
#   * ENTRE   -> las celdas suman cosas distintas de 0 y solo se cancelan al juntarlas.  Entonces lo
#                que organiza la cancelacion es una operacion sobre el POSET de order types, no una
#                involucion local.
#
# CONTROLES
#   C0  FATAL.  La suma de todas las celdas tiene que dar el A_{mu_max} que ya conocemos.
#   C1  se imprimen TODAS las celdas con su suma, no un resumen: con 9 celdas caben.
#   C2  SEÑUELO.  El mismo reparto agrupando por una particion FALSA del mismo tamaño -- el resto de
#       |Lambda| modulo el numero de celdas.  Si el azar concentra igual, el order type no explica
#       nada y la lectura no vale.
#   C3  no vacuidad: n por celda impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage ordertype_split.sage

import itertools, json, sys
from collections import defaultdict

t, r, R, N = 4, 2, 3, 8
INF = 10 ** 9


def phi_bialternante(beta, tt, nvar):
    Nn = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    Lr = LaurentPolynomialRing(K, nvar, 'z')
    zs = Lr.gens()
    x = [Lr(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(Nn - 1, -1, -1))
    def alt(expo):
        return matrix(Lr, Nn, Nn, lambda i, j: x[i] ** expo[j]).determinant()
    den = alt(delta)
    if den == 0:
        return None
    q = alt(list(beta)) / den
    try:
        q = Lr(q)
    except Exception:
        return "NO-POL"
    return {tuple(e) if hasattr(e, '__iter__') else (e,): c
            for e, c in zip(q.exponents(), q.coefficients()) if c != 0}


_SP = {}
def sp_char(mu, rr):
    key = (tuple(mu), rr)
    if key not in _SP:
        W = WeylCharacterRing("C%d" % rr)
        el = W(W.space().from_vector(vector(list(mu))))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + mult
        _SP[key] = d
    return _SP[key]


def pelar(P, rr, tope=6000):
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
        for k, v in sp_char(mu, rr).items():
            nv = P.get(k, 0) - c * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P


def entrelaza(mu, Lam):
    Lp = list(Lam) + [0]
    return all(Lp[i] >= mu[i] >= Lp[i + 2] for i in range(len(Lam) - 1))


def rs_locales(mu, Lam):
    L = list(Lam) + [0]
    M = [INF] + list(mu) + [0] * (len(Lam) + 2)
    return [min(L[i], M[i]) - max(L[i + 1], M[i + 1]) for i in range(len(Lam))]


def order_type(mu, Lam):
    """que rama toma el min y cual el max en cada posicion.  La celda de Yacobi."""
    L = list(Lam) + [0]
    M = [INF] + list(mu) + [0] * (len(Lam) + 2)
    return tuple((0 if L[i] <= M[i] else 1, 0 if L[i + 1] >= M[i + 1] else 1)
                 for i in range(len(Lam)))


def eps_de(rs):
    return 0 if any(x % 2 for x in rs) else (-1) ** (sum(rs) // 2)


CASOS = [(10, 9, 7, 4, 3, 2, 1, 0), (12, 11, 10, 5, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
         (13, 9, 8, 7, 5, 4, 2, 0), (14, 13, 11, 4, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0),
         (18, 17, 11, 8, 7, 6, 1, 0)]

print("=" * 122)
print("LA CANCELACION, SEGMENTADA POR ORDER TYPE   --   t=4, r=2, R=3")
print("=" * 122)
RES = []
for b in CASOS:
    Psi = phi_bialternante(b, 2, R)
    if Psi in (None, "NO-POL"):
        continue
    aL, _ = pelar(Psi, R)
    aL = {k: int(v) for k, v in aL.items() if v != 0}
    MUS = set()
    for Lam in aL:
        for u1 in range(Lam[0] + 1):
            for u2 in range(u1 + 1):
                if entrelaza((u1, u2), Lam):
                    MUS.add((u1, u2))
    A, sup = {}, {}
    for mu in MUS:
        tot, L = 0, []
        for Lam, a in aL.items():
            if not entrelaza(mu, Lam):
                continue
            e = eps_de(rs_locales(mu, Lam))
            if e:
                tot += a * e
                L.append((Lam, a * e))
        if tot:
            A[mu] = tot; sup[mu] = L
    S = list(A)
    maxi = [mu for mu in S if not any(nu != mu and all(sum(nu[:k + 1]) >= sum(mu[:k + 1])
                                                       for k in range(r)) for nu in S)]
    if len(maxi) != 1:
        continue
    mm = maxi[0]
    celdas = defaultdict(list)
    for Lam, c in sup[mm]:
        celdas[order_type(mm, Lam)].append((Lam, c))
    sumas = {tp: sum(c for _, c in L) for tp, L in celdas.items()}
    total = sum(sumas.values())
    nulas = sum(1 for v in sumas.values() if v == 0)
    # C2 SEÑUELO: agrupar por |Lambda| mod (numero de celdas), una particion falsa del mismo tamaño
    k = max(1, len(celdas))
    falsas = defaultdict(int)
    for Lam, c in sup[mm]:
        falsas[sum(Lam) % k] += c
    nulas_f = sum(1 for v in falsas.values() if v == 0)
    print("")
    print("  beta=%s   mu_max=%s   A=%d   %d supervivientes en %d celdas"
          % (str(b), str(mm), A[mm], len(sup[mm]), len(celdas)))
    for tp, L in sorted(celdas.items(), key=lambda kv: -len(kv[1])):
        print("      %-30s n=%2d  suma=%+d   %s"
              % (str(tp), len(L), sumas[tp], str([c for _, c in L])[:46]))
    print("      C0 total = %d, y A = %d : %s" % (total, A[mm], "ok" if total == A[mm] else "*** FALLA ***"))
    print("      celdas que suman 0: %d de %d   |   SEÑUELO (particion falsa): %d de %d"
          % (nulas, len(celdas), nulas_f, len(falsas)))
    sys.stdout.flush()
    RES.append({"beta": [int(x) for x in b], "A": int(A[mm]), "n_sup": len(sup[mm]),
                "n_celdas": len(celdas), "sumas": [int(v) for v in sumas.values()],
                "celdas_cero": int(nulas), "senuelo_cero": int(nulas_f),
                "senuelo_celdas": int(len(falsas))})

print("")
print("=" * 122)
print("  LECTURA, escrita ANTES de correr:")
print("   * si en cada forma TODAS las celdas suman 0 salvo una que suma +-1  ->  la cancelacion es")
print("     DENTRO del order type, y la involucion hay que buscarla dentro de la celda.")
print("   * si las celdas suman valores variados que solo se cancelan al juntarlas  ->  lo que")
print("     organiza la cancelacion vive en el POSET de order types, no en una accion local.")
print("   * y si el SEÑUELO concentra igual, el order type no explica nada.")
print("")
tot_c = sum(x["celdas_cero"] for x in RES); tot_n = sum(x["n_celdas"] for x in RES)
tot_f = sum(x["senuelo_cero"] for x in RES); tot_fn = sum(x["senuelo_celdas"] for x in RES)
print("  AGREGADO:  celdas que suman 0: %d de %d (%.1f %%)   |   señuelo: %d de %d (%.1f %%)"
      % (tot_c, tot_n, 100.0 * tot_c / max(1, tot_n), tot_f, tot_fn, 100.0 * tot_f / max(1, tot_fn)))
json.dump(RES, open("ordertype_split_DUMP.json", "w"), indent=1)
print("")
print("=" * 122)
print("DONE")
