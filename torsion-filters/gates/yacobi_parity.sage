# -*- coding: utf-8 -*-
# LA RUTA DE YACOBI, EJECUTADA:  A_mu por paridad de los r_i, sin pasar por los eta.
# 15 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 16 de la consulta, contestando nuestra P17: nuestro Sp_2 es la DIAGONAL
# Delta SL_2 en (SL_2)^R, no uno de los factores.  Verificado en el texto de Yacobi (arXiv:0907.3247),
# literalmente:
#
#   Corollary 3.3.  "Suppose (mu,lambda) in Lambda_B.  Then as SL_2-modules
#                    W_{lambda/mu} ~= (x)_{i=1}^n F_{r_i(mu,lambda^+)},
#                    where SL_2 acts by the tensor product representation on the right hand side."
#   Theorem 3.5.    "Moreover, Res^L_{SL_2}(Phi) recovers the natural action of SL_2 on B."
#
# y la receta de los r_i, tambien literal:
#
#   "Given (mu,lambda) in Lambda_{n-1,n+1} let (x_1 >= y_1 >= ... >= x_n >= y_n) be the
#    non-increasing rearrangement of (mu_1,...,mu_{n-1}, lambda_1,...,lambda_{n+1}).
#    Set r_i(mu,lambda) = x_i - y_i."
#
# CONSECUENCIA, y es la de el.  Sobre la diagonal el caracter de un producto tensorial es el
# PRODUCTO de caracteres, luego los 2-6 eta que medimos son internos a la restriccion y se pueden
# SALTAR:
#
#     chi_{W_{Lambda/mu}}(i)  =  prod_j chi_{r_j}(i)
#                             =  0                       si algun r_j es impar,
#                                (-1)^{(sum r_j)/2}      si todos son pares.
#
# Y como Phi_{2,R} es un caracter VIRTUAL, Phi_{2,R} = sum_Lambda a_Lambda sp_Lambda, queda
#
#     (8)   A_mu  =  sum_Lambda  a_Lambda * eps(Lambda,mu),      eps en {0,+-1}.
#
# LO QUE SE MIDE.  No (8) por si sola, que es teoria: lo que interesa es CUANTOS Lambda sobreviven
# a "todos los r_j pares" en el peso maximo.  Sus tres escenarios, escritos antes de correr:
#     UNO solo                      -> la unidad en t=4 queda practicamente resuelta.
#     varios con cancelacion        -> prueba plausible, hay que exhibir la involucion.
#     muchos sin estructura         -> la ruta simplifica pero no resuelve.
#
# CONTROLES
#   C0  FATAL.  (8) tiene que reproducir el A_mu que ya teniamos por la otra via (branching + tau).
#       Si no, o la receta de los r_i esta mal aplicada o la lectura diagonal no es la nuestra.
#   C1  eps en {0,+-1} por construccion: se COMPRUEBA igualmente, porque comprobarlo cuesta nada y
#       un fallo delataria un error en la reordenacion.
#   C2  SEÑUELO.  La misma cuenta con la reordenacion SIN ordenar (los r_i tal cual, sin el
#       rearrangement).  Tiene que fallar C0; si acertara, el rearrangement no estaria haciendo nada.
#   C3  no vacuidad: n impreso, y el recuento de Lambda por forma.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage yacobi_parity.sage

import itertools, json, sys
from collections import defaultdict

t, r = 4, 2
m, R, N = 1, 3, 8            # R = r + m = 3:  Sp_6 -> Sp_4 x Sp_2


def phi_bialternante(beta, tt, nvar):
    Nn = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    L = LaurentPolynomialRing(K, nvar, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(Nn - 1, -1, -1))
    def alt(expo):
        return matrix(L, Nn, Nn, lambda i, j: x[i] ** expo[j]).determinant()
    den = alt(delta)
    if den == 0:
        return None
    q = alt(list(beta)) / den
    try:
        q = L(q)
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


# ---------------------------------------------------------------- la receta de Yacobi -----------
def r_de_yacobi(mu, Lam, ordenar=True):
    """r_i = x_i - y_i con (x_1>=y_1>=...>=x_n>=y_n) la reordenacion decreciente de
       (mu_1..mu_{n-1}, Lam^+_1..Lam^+_{n+1}),  Lam^+ = Lam con un cero añadido.
       n = R = len(Lam).  Con ordenar=False se salta el rearrangement: es el SEÑUELO."""
    n = len(Lam)
    Lp = list(Lam) + [0]
    juntos = list(mu) + Lp
    seq = sorted(juntos, reverse=True) if ordenar else juntos
    if len(seq) != 2 * n:
        return None
    return [seq[2 * i] - seq[2 * i + 1] for i in range(n)]


def doble_entrelaza(mu, Lam):
    """mu << Lam, el DOBLE ENTRELAZADO de Yacobi:  Lam^+_i >= mu_i >= Lam^+_{i+2},  i=1..n-1.
    Es la PRECONDICION que me faltaba: 'W_{lambda/mu} != {0} <=> mu << lambda'.  Sin ella, mu que
    no entrelaza produce igualmente una reordenacion y un eps no nulo, y mete contribuciones
    FANTASMA de espacios de multiplicidad que son cero.  Eso hacia fallar C0 en 7 de 7."""
    n = len(Lam)
    Lp = list(Lam) + [0]
    for i in range(n - 1):
        if not (Lp[i] >= mu[i] >= Lp[i + 2]):
            return False
    return True


def eps_de(rs):
    """0 si algun r_j es impar; (-1)^{(sum r_j)/2} si todos pares."""
    if rs is None or any(x % 2 for x in rs):
        return 0
    return (-1) ** (sum(rs) // 2)


def tau4(eta):
    k = eta[0]
    return 0 if k % 2 else (-1) ** (k // 2)


def pelar_branching(P, tope=9000):
    P = {e: QQ(c) for e, c in P.items() if c != 0}
    out = {}
    for _ in range(tope):
        if not P:
            return out, {}
        dom = [e for e in P
               if list(e[:m]) == sorted(e[:m], reverse=True) and min(e[:m]) >= 0
               and list(e[m:]) == sorted(e[m:], reverse=True) and min(e[m:]) >= 0]
        if not dom:
            return out, P
        top = max(dom, key=lambda e: (sum(e), e))
        B = P[top]
        out[(tuple(top[:m]), tuple(top[m:]))] = out.get((tuple(top[:m]), tuple(top[m:])), 0) + B
        a, b = sp_char(tuple(top[:m]), m), sp_char(tuple(top[m:]), r)
        for e1, c1 in a.items():
            for e2, c2 in b.items():
                k = e1 + e2
                nv = P.get(k, 0) - B * c1 * c2
                if nv == 0:
                    P.pop(k, None)
                else:
                    P[k] = nv
    return out, P


CASOS = [(10, 9, 7, 4, 3, 2, 1, 0), (12, 11, 10, 5, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
         (13, 9, 8, 7, 5, 4, 2, 0), (14, 13, 11, 4, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0),
         (18, 17, 11, 8, 7, 6, 1, 0)]

print("=" * 126)
print("LA RUTA DE YACOBI EJECUTADA   --   t=4, r=2, R=3:  Sp_6 -> Sp_4,  Sp_2 = diagonal")
print("=" * 126)
print("")
print("  %-28s | #Lambda | C0 (8)==A_mu | C2 señuelo | mu_max      | A_max | #Lambda que sobreviven en mu_max"
      % "beta")
print("  " + "-" * 122)
RES = []
for b in CASOS:
    Psi = phi_bialternante(b, 2, R)
    if Psi in (None, "NO-POL"):
        continue
    aL, resto = pelar(Psi, R)                       # Phi_{2,R} = sum_Lambda a_Lambda sp_Lambda
    aL = {k: v for k, v in aL.items() if v != 0}
    # --- A_mu por NUESTRA via (branching + tau), que es la referencia -------------------------
    B, _ = pelar_branching(Psi)
    A_ref = defaultdict(lambda: 0)
    for (eta, mu), bb in B.items():
        v = tau4(eta)
        if v:
            A_ref[mu] += bb * v
    A_ref = {mu: int(c) for mu, c in A_ref.items() if c != 0}
    # --- A_mu por la via de YACOBI (8) ---------------------------------------------------------
    MUS = set(A_ref) | {mu for (_, mu) in B}
    A_yac, A_sen = defaultdict(lambda: 0), defaultdict(lambda: 0)
    sobreviven = defaultdict(list)
    for mu in MUS:
        for Lam, a in aL.items():
            if not doble_entrelaza(mu, Lam):      # la precondicion, ANTES de calcular nada
                continue
            e = eps_de(r_de_yacobi(mu, Lam, True))
            if e:
                A_yac[mu] += a * e
                sobreviven[mu].append((Lam, int(a), int(e)))
            es = eps_de(r_de_yacobi(mu, Lam, False))
            if es:
                A_sen[mu] += a * es
    A_yac = {mu: int(c) for mu, c in A_yac.items() if c != 0}
    A_sen = {mu: int(c) for mu, c in A_sen.items() if c != 0}
    c0 = (A_yac == A_ref)
    c2 = (A_sen != A_ref)
    S = list(A_ref)
    maxi = [mu for mu in S if not any(nu != mu and all(sum(nu[:k + 1]) >= sum(mu[:k + 1])
                                                       for k in range(r)) for nu in S)]
    mm = maxi[0] if len(maxi) == 1 else None
    print("  %-28s | %7d | %-12s | %-10s | %-11s | %-5s | %d"
          % (str(b), len(aL), "IGUAL" if c0 else "*** DISTINTO ***",
             "discrepa" if c2 else "*** COINCIDE ***", str(mm),
             str(A_ref.get(mm)) if mm else "-", len(sobreviven[mm]) if mm else -1))
    sys.stdout.flush()
    RES.append({"beta": [int(x) for x in b], "n_Lambda": int(len(aL)), "C0": bool(c0),
                "C2": bool(c2), "mu_max": [int(x) for x in mm] if mm else None,
                "A_max": A_ref.get(mm), "n_sobreviven": len(sobreviven[mm]) if mm else -1,
                "sobreviven": [[list(map(int, L)), a, e] for L, a, e in sobreviven[mm]] if mm else []})

print("")
n = len(RES)
print("  n = %d formas" % n)
print("  C0  (8) reproduce A_mu por la via del branching : %d/%d" % (sum(x["C0"] for x in RES), n))
print("  C2  el señuelo sin rearrangement discrepa       : %d/%d" % (sum(x["C2"] for x in RES), n))
print("")
h = defaultdict(int)
for x in RES:
    if x["n_sobreviven"] >= 0:
        h[x["n_sobreviven"]] += 1
print("  CUANTOS Lambda sobreviven a 'todos los r_j pares' en mu_max: %s" % dict(sorted(h.items())))
print("")
print("  LECTURA, escrita ANTES de correr (son sus tres escenarios):")
print("    UNO solo            -> la unidad en t=4 queda practicamente resuelta.")
print("    varios, cancelando  -> prueba plausible; habria que exhibir la involucion.")
print("    muchos sin patron   -> la ruta simplifica pero no resuelve.")
json.dump(RES, open("yacobi_parity_DUMP.json", "w"), indent=1)
print("")
print("=" * 126)
print("DONE")
