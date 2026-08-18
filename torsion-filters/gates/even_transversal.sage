# -*- coding: utf-8 -*-
# EL ANALOGO PAR DEL TRANSVERSAL.   16 de agosto de 2026.
#
# DE DONDE SALE.  prop:transversal deja el numerador impar como un recuento con signo de
# transversales.  El par tiene la misma forma de par de rango igual --- Sp_{2R} > Sp_{2m} x Sp_{2r},
# rango R = m+r en los dos lados --- asi que GKRS aplica tambien ahi.  Lo que cambia es CUANTO:
#
#   (a) las raices complementarias son solo  e_i +- f_j  (en tipo B habia ademas  f_j  suelta),
#       y congelando  y_i = xi^i  con  i = 1..m  y  t = 2m+2  faltan justo  +-1  del alfabeto, luego
#
#             denominador relativo especializado  =  prod_j  (1 - z_j^t) / (1 - z_j^2)
#
#       que NO es el del impar.
#   (b) las clases plegadas no nulas modulo t = 2m+2 son  {1, ..., m+1}  con  m+1 = t/2,  y
#       tau^C_t != 0 exige exactamente {1,...,m}: hay DOS clases prohibidas, 0 y t/2, no una.
#   (c) |W^1| = C(R,m), SIN el factor 2: los dos bloques son de tipo C y no hay quiralidad.
#
# Prediccion:  |supp nu| = prod_{j=1}^{m} n_j,  y  nu = 0  <=>  alguna clase de 1..m sin tocar.
#
# POR QUE IMPORTA.  El paper dice ahora, en la introduccion y en el cierre, que el lado par "no tiene
# esa reduccion".  Si el analogo funciona por Lambda, esa frase es demasiado fuerte y hay que
# corregirla.  (L1) es un enunciado por Lambda, y por Lambda el caracter par SI es genuino: lo virtual
# solo entra al sumar con los a_Lambda, que es (L3).
#
# LO QUE SE MIDE
#   E1  |supp nu| == prod n_j                                        FATAL
#   E2  supp nu == {transversales}, punto a punto
#   E3  los valores tambien
#   E4  nu = 0  <=>  alguna clase 1..m vacia
#   E5  |W^1| == C(R,m)
#   E6  LA IDENTIDAD:  ( sum_mu c(L,mu) chi^{C_r}_mu ) . Delta^par_t  ==  +- nu       FATAL
#       con c por branching Sp_{2R} -> Sp_{2m} x Sp_{2r} y el denominador de (a).
#   E7  el top lo da el transversal de v mas pequeno (analogo de la parte (v))
#
# CONTROLES / SENUELOS
#   C0  E1 y E6 son fatales.
#   D1  SENUELO: la regla IMPAR (solo la clase 0 prohibida).  Debe fallar, y su fallo es la
#       dicotomia del paper dicha en lenguaje de transversales.
#   D2  SENUELO: el denominador del impar, prod (1 - z_j^t), sin dividir por (1 - z_j^2).
#   C1  la constante de tipo C: tau^C == sgn(sigma) prod eps, sin constante.  Se comprueba.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage even_transversal.sage

import json
import sys
import itertools
from collections import Counter

CASOS = [(4, 2, 4), (6, 2, 3), (8, 2, 2), (6, 3, 2)]


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


def delta_C(a, t, m):
    """det de la permutacion con signo de las clases plegadas; 0 si no es permutacion de {1..m}."""
    cl, ep = [], []
    for v in a:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, m + 1)):
        return 0
    # convenio de lem:T: las clases se leen contra el orden DECRECIENTE m, m-1, ..., 1.
    # Leerlas contra el creciente multiplica por (-1)^{m(m-1)/2} y mueve la constante.
    s = sgn_perm([m - c for c in cl])
    for e in ep:
        s *= e
    return int(s)


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


def tauC_eval(eta, t, m):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car("C", m, eta).items():
        s += mult * z ** (sum((i + 1) * int(wt[i]) for i in range(m)) % t)
    return int(QQ(s)) if s in QQ else None


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
print("EL ANALOGO PAR DEL TRANSVERSAL:  Sp_{2R} > Sp_{2m} x Sp_{2r},  dos clases prohibidas")
print("=" * 112)
sys.stdout.flush()

RES = []
for (t, r, cota) in CASOS:
    m = (t - 2) // 2
    Rr = m + r
    if m < 1:
        continue
    n = e1 = e2 = e3 = e4 = e5 = e6 = e7 = d1 = d2 = c1 = 0
    n_vivas = 0   # E7 solo tiene sentido donde nu != 0: su denominador NO es n
    fallo = None
    for Lam in dominantes(Rr, cota):
        v = [int(Lam[i]) + Rr - i for i in range(Rr)]     # v = Lambda + rho_{C_R}
        # clases plegadas
        d = {}
        for i, x in enumerate(v):
            d.setdefault(plegar(x, t)[0], []).append(i)
        nj = [len(d.get(j, [])) for j in range(1, m + 1)]
        pred = 1
        for x in nj:
            pred *= x

        # ---- W^1 honesto para C_R > C_m x C_r : subconjunto S de tamano m, todo positivo
        reps = []
        for S in itertools.combinations(range(Rr), m):
            Sc = [i for i in range(Rr) if i not in S]
            A = sorted([v[i] for i in S], reverse=True)
            libre = sorted([v[i] for i in Sc], reverse=True)
            orden = sorted(S, key=lambda i: -v[i]) + sorted(Sc, key=lambda i: -v[i])
            sg = sgn_perm([orden.index(i) for i in range(Rr)])
            reps.append((tuple(A), tuple(libre), sg, frozenset(S)))
        n += 1
        e5 += 1 if len(reps) == binomial(Rr, m) else 0

        nu = {}
        for A, libre, sg, S in reps:
            tv = delta_C(A, t, m)
            if tv:
                nu[libre] = nu.get(libre, 0) + sg * tv
        nu = {k: val for k, val in nu.items() if val}

        e1 += 1 if len(nu) == pred else 0
        if len(nu) != pred and fallo is None:
            fallo = {"Lambda": list(Lam), "v": v, "n_j": nj, "pred": pred, "real": len(nu)}
        vacio = any(x == 0 for x in nj)
        e4 += 1 if ((len(nu) == 0) == vacio) else 0

        # E2/E3 por la lectura de transversales
        prev = {}
        if not vacio:
            for pick in itertools.product(*[d[j] for j in range(1, m + 1)]):
                S = frozenset(pick)
                if len(S) != m:
                    continue
                Sc = [i for i in range(Rr) if i not in S]
                A = sorted([v[i] for i in S], reverse=True)
                libre = tuple(sorted([v[i] for i in Sc], reverse=True))
                orden = sorted(S, key=lambda i: -v[i]) + sorted(Sc, key=lambda i: -v[i])
                sg = sgn_perm([orden.index(i) for i in range(Rr)])
                val = sg * delta_C(A, t, m)
                if val:
                    prev[libre] = val
        e2 += 1 if set(prev) == set(nu) else 0
        e3 += 1 if prev == nu else 0

        # E7 el top -- se puntua contra n_vivas, no contra n
        if nu:
            n_vivas += 1
            top = max(nu, key=lambda k: (sum(k), k))
            S_min = frozenset(max(d[j]) for j in range(1, m + 1))
            if len(S_min) == m:
                Sc = [i for i in range(Rr) if i not in S_min]
                e7 += 1 if tuple(sorted([v[i] for i in Sc], reverse=True)) == top else 0

        # D1 SENUELO: la regla impar, solo la clase 0 prohibida (las clases 1..m+1 admitidas)
        nj_impar = [len(d.get(j, [])) for j in range(1, m + 2)]
        p_impar = 1
        for x in nj_impar:
            p_impar *= x
        if len(nu) == p_impar:
            d1 += 1

        # C1 la constante de tipo C: tau^C == delta, sin constante
        for (eta, mu), b in list(branch(Rr, m, r, Lam).items())[:6]:
            tv = tauC_eval(eta, t, m)
            a0 = [int(eta[j]) + m - j for j in range(m)]
            if tv is not None and tv == delta_C(a0, t, m):
                c1 += 1
            elif tv is not None:
                if fallo is None:
                    fallo = {"C1": "la constante de tipo C no es +1", "eta": list(eta),
                             "tau": int(tv), "delta": int(delta_C(a0, t, m))}
    print("")
    print("  t=%d r=%d (m=%d R=%d)  %d pesos" % (t, r, m, Rr, n))
    print("     E1  |supp nu| == prod n_j (SIN el 2)   : %3d de %3d" % (e1, n))
    print("     E2  supp nu == transversales           : %3d de %3d" % (e2, n))
    print("     E3  y los valores tambien              : %3d de %3d" % (e3, n))
    print("     E4  nu = 0 <=> alguna clase 1..m vacia : %3d de %3d" % (e4, n))
    print("     E5  |W^1| == C(R,m)                    : %3d de %3d" % (e5, n))
    print("     E7  el top lo da S_min                 : %3d de %3d  (solo formas con nu != 0)"
          % (e7, n_vivas))
    print("     D1  SENUELO: la regla IMPAR            : %3d de %3d  (debe ser bajo)" % (d1, n))
    print("     C1  tau^C == delta, sin constante      : %d aciertos" % c1)
    if fallo:
        print("     !! %s" % json.dumps(fallo))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "m": int(m), "R": int(Rr), "n": int(n),
                "E1": int(e1), "E2": int(e2), "E3": int(e3),
                "E4": int(e4), "E5": int(e5), "E7": int(e7), "n_vivas": int(n_vivas),
                "D1": int(d1), "C1": int(c1), "fallo": fallo})

json.dump(RES, open("even_transversal_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr:")
print("   * si E1-E4 salen limpias, el numerador PAR es tambien un recuento de transversales, con")
print("     DOS clases prohibidas en vez de una y sin el factor de quiralidad.  Entonces la frase")
print("     del paper 'el lado par no tiene esa reduccion' es demasiado fuerte y hay que corregirla.")
print("   * si el senuelo D1 acierta, la distincion par/impar no esta en el numero de clases")
print("     prohibidas y la lectura es otra.")
print("=" * 112)
print("DONE")
