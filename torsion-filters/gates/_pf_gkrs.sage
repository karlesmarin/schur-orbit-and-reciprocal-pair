# -*- coding: utf-8 -*-
# GKRS PARA (L1).   16 de agosto de 2026.
#
# DE DONDE SALE.  El factor cruzado  z_j^t - 1  que mato la ruta de Laplace (rem:L1route) es, por
# prop:crossden, el denominador relativo del par de RANGO IGUAL  B_{R'} > B_{m'} x D_r.  Los pares de
# rango igual tienen su formula de caracteres: Gross-Kostant-Ramond-Sternberg.  La ruta que se sigue
# aqui es aplicar el filtro  tau^B_t  al factor  B_{m'}  de la identidad GKRS y preguntar por
# que dividir por el denominador especializado deja  0, +-1  en la base de  D_r.
#
# LA IDENTIDAD, especializada.   Con  v = Lambda + rho_{B_{R'}}  y  W^1  los  w  tales que  w(v)  es
# estrictamente H-dominante,
#
#     ch V_Lambda . (D_G / D_H)  =  sum_{w in W^1} (-1)^{l(w)} ch_{B_{m'}}(eta_w) (x) ch_{D_r}(mu_w)
#
# Congelando el bloque  y_i = zeta^i  el lado izquierdo es  ch V_Lambda(xi_t, z) . prod_j (z_j^{t/2} -
# z_j^{-t/2})  y el derecho pierde todo  w  con  tau^B_t(eta_w) = 0.  Escribiendo
# c(Lambda,mu) = sum_eta B_{Lambda;eta,mu} tau^B(eta)  queda
#
#     ( sum_mu c(Lambda,mu) chi^{D_r}_mu ) . Delta_t  =  sum_{w} sgn(w) tau^B_t(eta_w) chi^{D_r}_{mu_w}
#                                                        =: sum_mu nu(Lambda,mu) chi^{D_r}_mu
#
# y  nu  esta en  {0,+-1}  POR TEOREMA (GKRS multiplicity-free + [NPP25] Thm 4.1) salvo colisiones
# de mu_w.  Luego (L1) es exactamente la pregunta de si DIVIDIR por Delta_t conserva  {0,+-1}.
#
# LA CUENTA SE VUELVE ENTERA.  En la base de ALTERNANTES la division es explicita.  Coordenadas
# DOBLADAS  x := 2(mu + rho_{D_r}):
#
#     multiplicar por Delta_t  =  x -> x + t.epsilon,  epsilon in {+-1}^r,  con signo  prod epsilon_j,
#                                  seguido de enderezar por W(D_r) con su det.
#
# y el indice alternante del lado GKRS es LITERALMENTE  x_w = w(v) restringido al bloque libre (doblado).
#
# LO QUE SE MIDE
#   T1  CONTROL FATAL:  sum_mu c(L,mu) . [las 2^r copias desplazadas y enderezadas]  ==  nu.
#       Si esto falla, o GKRS no es lo que creemos o prop:crossden no es el denominador.
#   T2  el reparto de nu.  Y si los mu_w son dos a dos distintos (si lo son, nu en {0,+-1} sin medir).
#   T3  COLISIONES: las 2^r imagenes desplazadas de supp(c), .se pisan?  Si NO se pisan, entonces
#       c(L,mu) = +- nu(enderezar(2mu+2rho+t.1))  y (L1) queda PROBADA a partir de nu en {0,+-1}.
#   T4  esa reconstruccion, medida: c leido desde nu por la copia todo-+.
#
# CONTROLES / SENUELOS
#   C0  T1 es fatal, forma a forma.
#   D1  senuelo del denominador: desplazar por  t+2  en vez de  t  (mismo paridad, otro numero).
#   D2  senuelo del signo: olvidar  prod epsilon_j.
#   C2  se imprimen n y los peores casos, no solo agregados.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage gkrs_L1.sage

import json
import sys
import itertools
from collections import Counter, defaultdict

# (t, r, cota de Lambda_1)
CASOS = [(3, 2, 2)]


# ----------------------------------------------------------------------------- caracteres y branching
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
        _TAU[key] = int(QQ(s)) if s in QQ else None
    return _TAU[key]


# ----------------------------------------------------------------------------- W(D_r): enderezar
def enderezar(x):
    """x en coordenadas DOBLADAS.  Devuelve (x_dominante_regular, det w) o None si es singular."""
    r = len(x)
    a = [abs(int(v) )for v in x]
    if len(set(a)) != r:                      # x_i = +- x_j  ->  alternante nulo
        return None
    idx = sorted(range(r), key=lambda i: -a[i])
    # signo de la permutacion
    perm = list(idx)
    s = 1
    seen = [False] * r
    for i in range(r):
        if seen[i]:
            continue
        j, L = i, 0
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            s = -s
    y = [int(x[i]) for i in idx]
    neg = sum(1 for v in y if v < 0)
    hay_cero = any(v == 0 for v in y)
    y = [abs(v) for v in y]
    if neg % 2 == 1 and not hay_cero:
        y[-1] = -y[-1]
    # los cambios de signo usados son SIEMPRE en numero par -> det de esa parte = +1
    return (tuple(y), s)


def desplazar(x, paso, con_signo=True):
    """las 2^r copias  x + paso.epsilon  enderezadas, como dict {x_dom: coef}."""
    r = len(x)
    out = defaultdict(lambda: 0)
    for eps in itertools.product((1, -1), repeat=r):
        sg = 1
        for e in eps:
            sg *= e
        if not con_signo:
            sg = 1
        y = tuple(int(x[j]) + paso * eps[j] for j in range(r))
        e = enderezar(y)
        if e is None:
            continue
        out[e[0]] += sg * e[1]
    return {k: v for k, v in out.items() if v != 0}


# ----------------------------------------------------------------------------- W^1 del par de rango igual
def coset_reps(v, mp, rr):
    """v = 2(Lambda + rho_{B_{R'}}), enteros impares decrecientes > 0.
       Devuelve [(u, sgn)] con u = w(v) estrictamente H-dominante."""
    Rp = mp + rr
    out = []
    for perm in itertools.permutations(range(Rp)):
        # signo de la permutacion
        s = 1
        seen = [False] * Rp
        for i in range(Rp):
            if seen[i]:
                continue
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                L += 1
            if L % 2 == 0:
                s = -s
        base = [v[perm[i]] for i in range(Rp)]
        for eps in itertools.product((1, -1), repeat=Rp):
            u = [base[i] * eps[i] for i in range(Rp)]
            # B_{m'} estrictamente dominante:  u_1 > ... > u_{m'} > 0
            ok = all(u[i] > u[i + 1] for i in range(mp - 1)) and (mp == 0 or u[mp - 1] > 0)
            if not ok:
                continue
            # D_r estrictamente dominante:  u_{m'+1} > ... > u_{R'-1} > |u_{R'}|
            f = u[mp:]
            if rr >= 2:
                ok = all(f[i] > f[i + 1] for i in range(rr - 1)) and f[rr - 2] > abs(f[rr - 1])
                if not ok:
                    continue
            sg = s
            for e in eps:
                sg *= e
            out.append((tuple(u), sg))
    return out


def dominantes(rk, cota):
    """Lambda dominante integral para B_rk con Lambda_1 <= cota."""
    def rec(k, tope):
        if k == 0:
            yield ()
            return
        for a in range(tope, -1, -1):
            for resto in rec(k - 1, a):
                yield (a,) + resto
    return list(rec(rk, cota))


# ----------------------------------------------------------------------------- el barrido
print("=" * 118)
print("GKRS PARA (L1):  la identidad de rango igual  B_{R'} > B_{m'} x D_r  filtrada por tau^B_t")
print("=" * 118)
sys.stdout.flush()

RES = []
for (t, r, cota) in CASOS:
    mp = (t - 1) // 2
    Rp = mp + r
    rho2 = [2 * (Rp - i) - 1 for i in range(Rp)]        # 2 rho_{B_{R'}} = (2R'-1, ..., 1)
    rhoD2 = [2 * (r - j) - 2 for j in range(r)]          # 2 rho_{D_r} = (2r-2, ..., 0)

    agg = {"n": 0, "T1": 0, "nu_rep": Counter(), "mu_w_distintos": 0, "sin_colision": 0,
           "T4": 0, "c_rep": Counter(), "D1": 0, "D2": 0, "fallos": [], "colisiones": []}

    for Lam in dominantes(Rp, cota):
        v = [2 * Lam[i] + rho2[i] for i in range(Rp)]

        # ---- c(Lambda, mu) por branching + filtro    (mu SIN plegar el signo quiral)
        c = defaultdict(lambda: 0)
        for (eta, mu), b in branch(Rp, mp, r, Lam).items():
            tv = tauB(mp, eta, t)
            if not tv:
                continue
            c[tuple(int(u) for u in mu)] += int(b) * int(tv)
        c = {k: v2 for k, v2 in c.items() if v2 != 0}

        # ---- nu(Lambda, mu) por GKRS
        reps = coset_reps(v, mp, r)
        nu = defaultdict(lambda: 0)
        mu_w = []
        for u, sg in reps:
            eta = tuple((u[i] - (2 * (mp - i) - 1)) // 2 for i in range(mp))
            tv = tauB(mp, eta, t)
            if not tv:
                continue
            x = tuple(u[mp:])                      # = 2(mu_w + rho_{D_r}), ya dominante regular
            mu_w.append(x)
            nu[x] += sg * int(tv)
        nu = {k: v2 for k, v2 in nu.items() if v2 != 0}

        # ---- T1: c . Delta_t  ==  nu     (en la base de alternantes, coordenadas dobladas)
        lhs = defaultdict(lambda: 0)
        for mu, cv in c.items():
            x = tuple(2 * int(mu[j]) + rhoD2[j] for j in range(r))
            for y, sg in desplazar(x, t).items():
                lhs[y] += cv * sg
        lhs = {k: v2 for k, v2 in lhs.items() if v2 != 0}
        ok1 = (lhs == nu)

        # ---- senuelos
        d1 = defaultdict(lambda: 0)
        for mu, cv in c.items():
            x = tuple(2 * int(mu[j]) + rhoD2[j] for j in range(r))
            for y, sg in desplazar(x, t + 2).items():
                d1[y] += cv * sg
        okD1 = ({k: v2 for k, v2 in d1.items() if v2 != 0} == nu)
        d2 = defaultdict(lambda: 0)
        for mu, cv in c.items():
            x = tuple(2 * int(mu[j]) + rhoD2[j] for j in range(r))
            for y, sg in desplazar(x, t, con_signo=False).items():
                d2[y] += cv * sg
        okD2 = ({k: v2 for k, v2 in d2.items() if v2 != 0} == nu)

        # ---- T3: colisiones entre las 2^r copias desplazadas de supp(c)
        golpes = Counter()
        for mu in c:
            x = tuple(2 * int(mu[j]) + rhoD2[j] for j in range(r))
            for eps in itertools.product((1, -1), repeat=r):
                y = tuple(x[j] + t * eps[j] for j in range(r))
                e = enderezar(y)
                if e is not None:
                    golpes[e[0]] += 1
        sin_col = all(g == 1 for g in golpes.values())

        # ---- T4: leer c desde nu por la copia todo-+
        ok4 = True
        for mu, cv in c.items():
            x = tuple(2 * int(mu[j]) + rhoD2[j] + t for j in range(r))
            e = enderezar(x)
            if e is None or nu.get(e[0], 0) != cv * e[1]:
                ok4 = False
                break

        agg["n"] += 1
        agg["T1"] += 1 if ok1 else 0
        agg["D1"] += 1 if okD1 else 0
        agg["D2"] += 1 if okD2 else 0
        agg["sin_colision"] += 1 if sin_col else 0
        agg["T4"] += 1 if ok4 else 0
        agg["mu_w_distintos"] += 1 if len(mu_w) == len(set(mu_w)) else 0
        for v2 in nu.values():
            agg["nu_rep"][int(v2)] += 1
        for v2 in c.values():
            agg["c_rep"][int(v2)] += 1
        if not ok1 and len(agg["fallos"]) < 4:
            agg["fallos"].append({"Lambda": [int(a) for a in Lam],
                                  "c": {str(k): int(v2) for k, v2 in c.items()},
                                  "nu": {str(k): int(v2) for k, v2 in nu.items()},
                                  "lhs": {str(k): int(v2) for k, v2 in lhs.items()}})
        if not sin_col and len(agg["colisiones"]) < 4:
            agg["colisiones"].append({"Lambda": [int(a) for a in Lam],
                                      "golpes": {str(k): int(v2) for k, v2 in golpes.items() if v2 > 1},
                                      "c": {str(k): int(v2) for k, v2 in c.items()}})

    print("")
    print("  t=%d  r=%d  (m'=%d, R'=%d)   Lambda_1 <= %d   ->  %d pesos" % (t, r, mp, Rp, cota, agg["n"]))
    print("     T1  CONTROL FATAL   c . Delta_t == nu            : %d de %d" % (agg["T1"], agg["n"]))
    print("     D1  senuelo  paso t+2                            : %d de %d  (debe ser ~0)" % (agg["D1"], agg["n"]))
    print("     D2  senuelo  sin el signo prod eps               : %d de %d  (debe ser ~0)" % (agg["D2"], agg["n"]))
    print("     T2  reparto de nu                                : %s" % dict(sorted(agg["nu_rep"].items())))
    print("         mu_w dos a dos distintos                     : %d de %d" % (agg["mu_w_distintos"], agg["n"]))
    print("     T3  las 2^r copias de supp(c) NO se pisan        : %d de %d" % (agg["sin_colision"], agg["n"]))
    print("     T4  c(L,mu) = +- nu(mu + t/2 . 1)                : %d de %d" % (agg["T4"], agg["n"]))
    print("         reparto de c                                 : %s" % dict(sorted(agg["c_rep"].items())))
    if agg["fallos"]:
        print("     !! fallos de T1: %s" % json.dumps(agg["fallos"][:2]))
    if agg["colisiones"]:
        print("     ** colisiones: %s" % json.dumps(agg["colisiones"][:2]))
    sys.stdout.flush()

    RES.append({"t": int(t), "r": int(r), "mp": int(mp), "Rp": int(Rp), "cota": int(cota),
                "n": int(agg["n"]), "T1": int(agg["T1"]), "D1": int(agg["D1"]), "D2": int(agg["D2"]),
                "nu_rep": {str(k): int(v2) for k, v2 in agg["nu_rep"].items()},
                "c_rep": {str(k): int(v2) for k, v2 in agg["c_rep"].items()},
                "mu_w_distintos": int(agg["mu_w_distintos"]),
                "sin_colision": int(agg["sin_colision"]), "T4": int(agg["T4"]),
                "fallos": agg["fallos"], "colisiones": agg["colisiones"]})

print("")
print("=" * 118)
print("  LECTURA, escrita ANTES de correr:")
print("   * si T1 acierta y los senuelos no, la identidad GKRS especializada ES prop:crossden y la")
print("     reformulacion de (L1) esta justificada como CUENTA, no solo como analogia.")
print("   * si ademas T3 sale limpio (sin colisiones) y nu vive en {0,+-1}, entonces")
print("         c(Lambda,mu) = +- nu(mu + (t/2).1)  en  {0,+-1}   -- eso es (L1) PROBADA.")
print("   * si hay colisiones, esos Lambda son EXACTAMENTE la obstruccion, y hay que mirar si las")
print("     copias que se pisan se cancelan (entonces c=0) o se suman (entonces (L1) peligra).")
json.dump(RES, open("_pf_gkrs_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
