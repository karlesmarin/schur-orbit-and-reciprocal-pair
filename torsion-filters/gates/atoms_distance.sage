# -*- coding: utf-8 -*-
# LOS ATOMOS (nu,Q) Y SUS DISTANCIAS.   15 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 20: atomizar por el lado que NO es el branching.  Con
#
#     A_mu = sum_nu eps_{lambda,nu}^{(t)} * m_{nu,mu},     eps en {0,+-1},
#
# el signed set es  X = {(nu,Q) : eps_nu != 0, Q cuenta la multiplicidad m_{nu,mu}}, y cada atomo
# pesa +-1.  El experimento que pide: construirlo y medir la distancia combinatoria minima entre
# atomos de signo OPUESTO, sin preguntar todavia si hay involucion.
#
# LA SIMPLIFICACION QUE LO HACE MEDIBLE HOY.  eps depende SOLO de nu, no de Q: todos los atomos con
# el mismo nu llevan el mismo signo y hay m_{nu,mu} de ellos.  Luego un emparejamiento que invierta
# el signo tiene que casar un nu POSITIVO con uno NEGATIVO, y la geometria que importa primero es la
# de los nu, no la de los Q.  Con eso no hace falta implementar los Rec de Watanabe para la primera
# medida -- solo para refinarla despues.
#
# Y SU PISTA DEL RIBBON NO ES ANALOGIA.  eps_nu = s_{lambda/nu}(mu_t) y, por Littlewood, eso es el
# signo de un teselado de lambda/nu por t-RIBBONS: vale 0 si no se puede teselar y +-1 con el signo
# de las alturas si se puede.  Asi que el movimiento natural entre nu supervivientes ES quitar o
# poner un t-ribbon.  El tamaño del paso lo dicta la estructura, no un ajuste.
#
# LO QUE SE MIDE
#   D1  para cada forma: los nu con eps != 0, su signo, y su multiplicidad m_{nu,mu}.
#   D2  la distancia minima |nu \ nu'| + |nu' \ nu| entre nu de signo opuesto.
#   D3  la pregunta afilada: ¿los pares opuestos mas cercanos difieren por UN t-ribbon?  Se mide si
#       nu/nu' es un ribbon conexo de t cajas.
#   D4  si esas distancias son todas iguales, hay una operacion local repetible; si son caoticas,
#       se abandona la involucion y se pasa a la interpretacion de caracteristica de Euler.
#
# CONTROLES
#   C0  FATAL.  sum_nu eps_nu m_{nu,mu} tiene que dar el A_mu conocido por la otra via.
#   C1  eps calculado por el determinante 0/1 se contrasta con el signo de ribbon, que es
#       independiente: si discrepan, uno de los dos esta mal.
#   C2  SEÑUELO para D3: se mide tambien cuantos pares opuestos difieren por un ribbon de tamaño
#       t' != t.  Si sale lo mismo, el tamaño t no esta diciendo nada.
#   C3  no vacuidad: n de atomos y de nu impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage atoms_distance.sage

import itertools, json, sys
from collections import defaultdict

Sym = SymmetricFunctions(QQ)
sch = Sym.schur()
t, r = 4, 2
m, R, N = 1, 3, 8
DELTA = list(range(N - 1, -1, -1))


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


def pelar_sp(P, rr, tope=4000):
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


def eps_det(lam, nu, tt):
    """eps = det([ tt divide lam_i - nu_j - i + j ]), el determinante 0/1."""
    n = len(lam)
    nn = list(nu) + [0] * (n - len(nu))
    M = matrix(ZZ, n, n, lambda i, j:
               1 if (lam[i] - nn[j] - i + j) >= 0 and (lam[i] - nn[j] - i + j) % tt == 0 else 0)
    return M.determinant()


_REC = {}
def s_nu_rec(nu, rr):
    key = (tuple(nu), rr)
    if key not in _REC:
        L = LaurentPolynomialRing(QQ, rr, 'z')
        zs = L.gens()
        alf = [g ** e for g in zs for e in (1, -1)]
        pol = sch[Partition(list(nu))].expand(2 * rr)
        q = L(pol(*alf))
        _REC[key] = {tuple(e) if hasattr(e, '__iter__') else (e,): c
                     for e, c in zip(q.exponents(), q.coefficients()) if c != 0}
    return _REC[key]


def celdas(p):
    return {(i, j) for i, x in enumerate(p) for j in range(x)}


def es_ribbon(a, b, k):
    """b subconjunto de a; ¿a/b es un ribbon CONEXO de k cajas (sin cuadrado 2x2)?"""
    A, B = celdas(a), celdas(b)
    if not B <= A:
        return False
    D = A - B
    if len(D) != k:
        return False
    if any((i, j) in D and (i + 1, j) in D and (i, j + 1) in D and (i + 1, j + 1) in D for i, j in D):
        return False
    # conexo por lados
    v = {next(iter(D))}
    cam = True
    while cam:
        cam = False
        for (i, j) in list(D - v):
            if any(((i + di, j + dj) in v) for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1))):
                v.add((i, j)); cam = True
    return v == D


def dist(a, b):
    A, B = celdas(a), celdas(b)
    return len(A - B) + len(B - A)


CASOS = [(10, 9, 7, 4, 3, 2, 1, 0), (12, 11, 10, 5, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
         (13, 9, 8, 7, 5, 4, 2, 0), (14, 13, 11, 4, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0),
         (18, 17, 11, 8, 7, 6, 1, 0)]

print("=" * 124)
print("LOS ATOMOS (nu,Q) Y SUS DISTANCIAS   --   t=%d, r=%d" % (t, r))
print("=" * 124)
RES = []
for b in CASOS:
    lam = tuple(x for x in (b[i] - DELTA[i] for i in range(N)) if x > 0)
    Phi = phi_bialternante(b, t, r)
    if Phi in (None, "NO-POL"):
        continue
    # los nu candidatos: contenidos en lambda, con a lo sumo 2r partes
    NUS = []
    for k in range(sum(lam) + 1):
        for nu in Partitions(k, max_length=2 * r).list():
            nu = tuple(nu)
            if len(nu) <= len(lam) and all(nu[i] <= lam[i] for i in range(len(nu))):
                NUS.append(nu)
    eps = {nu: int(eps_det(lam, nu, t)) for nu in NUS}
    eps = {k: v for k, v in eps.items() if v != 0}
    # A_mu por este lado
    A = defaultdict(lambda: 0)
    B = {}
    for nu, e in eps.items():
        bb, _ = pelar_sp(s_nu_rec(nu, r), r)
        B[nu] = {k: int(v) for k, v in bb.items() if v != 0}
        for mu, c in B[nu].items():
            A[mu] += e * c
    A = {mu: int(c) for mu, c in A.items() if c != 0}
    S = list(A)
    maxi = [mu for mu in S if not any(nu2 != mu and all(sum(nu2[:k + 1]) >= sum(mu[:k + 1])
                                                        for k in range(r)) for nu2 in S)]
    if len(maxi) != 1:
        continue
    mm = maxi[0]
    # los atomos: (nu, indice 1..m) con signo eps_nu
    at = [(nu, eps[nu], B[nu].get(mm, 0)) for nu in eps if B[nu].get(mm, 0)]
    pos = [(nu, mult) for nu, e, mult in at if e > 0]
    neg = [(nu, mult) for nu, e, mult in at if e < 0]
    n_at = sum(mult for _, _, mult in at)
    # D2/D3: distancias entre nu de signo opuesto
    dists = defaultdict(int)
    rib_t = rib_o = 0
    mind = None
    for nu1, _ in pos:
        for nu2, _ in neg:
            d = dist(nu1, nu2)
            dists[d] += 1
            mind = d if mind is None else min(mind, d)
            a1, a2 = (nu1, nu2) if sum(nu1) > sum(nu2) else (nu2, nu1)
            if es_ribbon(a1, a2, t):
                rib_t += 1
            for tp in (2, 3, 6):
                if tp != t and es_ribbon(a1, a2, tp):
                    rib_o += 1
    print("")
    print("  beta=%s  mu_max=%s  A=%d" % (str(b), str(mm), A[mm]))
    print("      nu con eps != 0 que tocan mu_max: %d  (%d positivos, %d negativos)"
          % (len(at), len(pos), len(neg)))
    print("      atomos totales (con multiplicidad): %d" % n_at)
    print("      distancias |nu triangulo nu'| entre signos opuestos: %s" % dict(sorted(dists.items())))
    print("      de esos pares, difieren por UN %d-ribbon: %d   |  SEÑUELO otros tamaños: %d"
          % (t, rib_t, rib_o))
    sys.stdout.flush()
    RES.append({"beta": [int(x) for x in b], "A": int(A[mm]), "n_nu": len(at),
                "n_atomos": int(n_at), "dists": {str(k): int(v) for k, v in dists.items()},
                "ribbon_t": int(rib_t), "ribbon_otros": int(rib_o)})

print("")
print("=" * 124)
print("  LECTURA, escrita ANTES de correr:")
print("   * si casi todos los pares opuestos cercanos difieren por UN t-ribbon -> hay operacion local")
print("     repetible, y la involucion tiene donde agarrarse.")
print("   * si las distancias son caoticas y el señuelo acierta igual -> se abandona la involucion y")
print("     se pasa a la interpretacion de caracteristica de Euler.")
json.dump(RES, open("atoms_distance_DUMP.json", "w"), indent=1)
print("")
print("=" * 124)
print("DONE")
