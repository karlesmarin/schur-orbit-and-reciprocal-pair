# -*- coding: utf-8 -*-
# ¿ES EL CONJUNTO DE SUPERVIVIENTES UN TORSOR DEL GRUPO DE WEYL?   16 de agosto de 2026.
#
# DE DONDE SALE.  Cruce de dos cosas ya establecidas:
#
#   (a)  el filtro solo ve  a = eta + rho  MODULO t  (corolario de periodicidad), y si sobrevive las
#        clases min(c_j, t-c_j) son una PERMUTACION de {1..m} con un vector de signos;
#   (b)  el signo es  sgn(sigma) prod eps_j.
#
# O sea: cada superviviente determina un elemento del grupo de Weyl  W(C_m) = {+-1}^m x| S_m,  y tau
# es su CARACTER SIGNO.  Si el conjunto de supervivientes que contribuyen a mu_max fuera estable bajo
# multiplicar por una reflexion, esa reflexion seria una involucion que INVIERTE EL SIGNO -- que es
# exactamente el objeto que llevamos dias buscando para la conjetura de la unidad.  Y lo que la
# obstruiria seria que no respeta los coeficientes  c_eta = sum_Lambda a_Lambda B_{Lambda;eta,mu}.
#
# Asi que no se pregunta "¿existe una involucion?" sino la version con dientes: se IMPRIME la tabla
# (elemento de Weyl, coeficiente, signo) y se mira qué estructura tiene.
#
# LO QUE SE MIDE, por forma
#   W1  los eta que contribuyen a mu_max con tau != 0, su c_eta y su elemento w de W(C_2).
#   W2  CONTROL FATAL: sum c_eta tau(eta) = A_{mu_max} calculado por la ruta bialternante.
#   W3  ¿el conjunto de w que aparecen es estable bajo alguna de las reflexiones de W?  Se prueban
#       TODAS las reflexiones y se dice cuales, y si ademas conservan |c_eta|.
#   W4  el reparto de signos y la suma, para ver si la cancelacion es por parejas o colectiva.
#
# CONTROLES
#   C0  W2 es fatal.
#   C1  se imprimen los w, no un resumen: la estructura hay que verla.
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage weyl_structure.sage

import json
import sys
import itertools
from collections import defaultdict

t, r, m = 6, 2, 2
R = m + r
N = t + 2 * r


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


def pelar(P, typ, rk, tope=9000):
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
        for k, v in car(typ, rk, mu).items():
            nv = P.get(k, 0) - c * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P


_BR = {}
def branch(Lam):
    key = tuple(int(v) for v in Lam)
    if key not in _BR:
        W = WeylCharacterRing("C%d" % R)
        X = WeylCharacterRing("C%dxC%d" % (m, r))
        br = branching_rule("C%d" % R, "C%dxC%d" % (m, r), "orthogonal_sum")
        el = W(W.space().from_vector(vector([Integer(v) for v in Lam])))
        d = {}
        for wt, c in el.branch(X, rule=br).monomial_coefficients().items():
            v = [int(u) for u in wt.to_vector()]
            d[(tuple(v[:m]), tuple(v[m:]))] = int(c)
        _BR[key] = d
    return _BR[key]


def tau(eta):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car("C", m, eta).items():
        s += mult * z ** (sum((i + 1) * wt[i] for i in range(m)) % t)
    return QQ(s) if s in QQ else None


def weyl_de(eta):
    """el elemento de W(C_m) = permutacion con signos que el patron de residuos determina."""
    a = [int(eta[j]) + m - j for j in range(m)]
    c = [v % t for v in a]
    cl, eps = [], []
    for v in c:
        if v <= m:
            cl.append(v); eps.append(1)
        else:
            cl.append(t - v); eps.append(-1)
    if sorted(cl) != list(range(1, m + 1)):
        return None
    return (tuple(cl), tuple(eps))


# todas las reflexiones de W(C_2), como acciones sobre (cl, eps)
def reflexiones(m):
    ref = []
    for i in range(m):                       # cambio de signo en la posicion i
        ref.append(("s_%d (signo)" % (i + 1), lambda w, i=i: (w[0],
                    tuple(-v if k == i else v for k, v in enumerate(w[1])))))
    for i in range(m):
        for j in range(i + 1, m):            # transposicion de las clases i, j
            ref.append(("t_%d%d (transpos)" % (i + 1, j + 1),
                        lambda w, i=i, j=j: (tuple(w[0][j] if k == i else w[0][i] if k == j else v
                                                   for k, v in enumerate(w[0])),
                                             tuple(w[1][j] if k == i else w[1][i] if k == j else v
                                                   for k, v in enumerate(w[1])))))
    return ref


CASOS = [(13, 12, 11, 10, 7, 5, 4, 3, 2, 0), (12, 11, 10, 9, 8, 7, 4, 3, 2, 0),
         (14, 13, 11, 10, 9, 6, 4, 3, 1, 0), (11, 10, 9, 8, 7, 6, 4, 2, 1, 0)]

print("=" * 118)
print("LOS SUPERVIVIENTES EN mu_max COMO ELEMENTOS DEL GRUPO DE WEYL   (t=%d, r=%d, m=%d)" % (t, r, m))
print("=" * 118)

RES = []
for b in CASOS:
    P = phi(b, t, r)
    if not P:
        print("  %-30s  Phi = 0" % str(b))
        continue
    A, rest = pelar(P, "C", r)
    A = {k: QQ(v) for k, v in A.items() if v != 0}
    if rest or not A:
        continue
    S = list(A)
    maxi = [x for x in S if not any(y != x and all(sum(y[:k + 1]) >= sum(x[:k + 1])
                                                   for k in range(r)) for y in S)]
    if len(maxi) != 1:
        print("  %-30s  mu_max no unico" % str(b))
        continue
    mm = maxi[0]
    P2 = phi(b, 2, R)
    aL, rest2 = pelar(P2, "C", R)
    if rest2:
        print("  %-30s  el pelado en C_%d no cierra" % (str(b), R))
        continue
    ceta = defaultdict(lambda: 0)
    for Lam, a in aL.items():
        for (eta, mu), c in branch(Lam).items():
            if mu != mm:
                continue
            ceta[eta] += int(a) * int(c)
    vivos = []
    for eta, c in ceta.items():
        if c == 0:
            continue
        tv = tau(eta)
        if tv is None or tv == 0:
            continue
        w = weyl_de(eta)
        vivos.append((eta, c, int(tv), w))
    suma = sum(c * tv for _, c, tv, _ in vivos)
    print("")
    print("  beta = %s   mu_max = %s   A = %s   (control: suma = %s)"
          % (str(b), str(mm), str(A[mm]), str(suma)))
    print("     eta        c_eta   tau   elemento de Weyl (clases | signos)")
    for eta, c, tv, w in sorted(vivos, key=lambda x: (-abs(x[1]), x[0])):
        print("     %-10s %6d %5d   %s" % (str(eta), c, tv, str(w) if w else "-- no regular --"))
    # W3: reflexiones que preservan el conjunto de w
    Wset = set(w for _, _, _, w in vivos if w)
    conserva = []
    for nombre, f in reflexiones(m):
        if set(f(w) for w in Wset) == Wset:
            # ¿y conserva |c_eta|?
            mapa = {}
            for eta, c, tv, w in vivos:
                if w:
                    mapa[w] = abs(c)
            ok = all(mapa.get(f(w)) == mapa[w] for w in Wset)
            conserva.append((nombre, ok))
    print("     reflexiones que preservan el conjunto: %s"
          % (str(conserva) if conserva else "NINGUNA"))
    print("     signos: +1 en %d, -1 en %d   |  suma de c_eta con signo = %s"
          % (sum(1 for _, _, tv, _ in vivos if tv > 0),
             sum(1 for _, _, tv, _ in vivos if tv < 0), str(suma)))
    sys.stdout.flush()
    RES.append({"beta": [int(v) for v in b], "mu_max": [int(v) for v in mm], "A": int(A[mm]),
                "suma": int(suma), "n_vivos": int(len(vivos)),
                "vivos": [[list(map(int, e)), int(c), int(tv), [list(map(int, w[0])), list(map(int, w[1]))] if w else None]
                          for e, c, tv, w in vivos],
                "reflexiones": [[n, bool(o)] for n, o in conserva]})

print("")
print("=" * 118)
print("  LECTURA, escrita ANTES de correr:")
print("   * si alguna reflexion preserva el conjunto Y los |c_eta|, esa reflexion ES una involucion")
print("     que invierte el signo, y la conjetura de la unidad se reduce a contar los puntos fijos.")
print("   * si preserva el conjunto pero NO los coeficientes, la obstruccion esta localizada en B.")
print("   * si ninguna lo preserva, el grupo de Weyl no organiza la cancelacion y hay que volver a")
print("     los atomos.")
json.dump(RES, open("weyl_structure_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
