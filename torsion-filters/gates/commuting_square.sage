# -*- coding: utf-8 -*-
# EL CUADRADO CONMUTATIVO:  las dos factorizaciones del mismo objeto.   15 de agosto de 2026.
#
# DE DONDE SALE.  Reseña del Paper II, propuesta 2.  El mismo objeto admite dos ordenes de proceso:
#
#   CAMINO A (el del companion) -- primero la RAIZ, despues el alfabeto libre:
#       s_lambda(mu_t u W) = sum_nu eps_nu s_nu(W),        eps_nu en {0,+-1}
#   y luego, restringiendo GL_{2r} -> Sp_{2r} sobre el alfabeto reciproco,
#       s_nu(z,1/z,...) = sum_mu b_{nu,mu} sp_mu(z).
#
#   CAMINO B (el de este paper) -- primero RAMIFICAR, despues filtrar:
#       A_mu = sum_eta B_{eta,mu} tau_t(eta).
#
# Los dos calculan el mismo A_mu, luego
#
#       (*)   sum_nu eps_nu b_{nu,mu}  =  sum_eta B_{eta,mu} tau_t(eta)   para todo mu.
#
# A la izquierda: Littlewood/raiz de la unidad + Littlewood restriction.  A la derecha: branching
# simplectico + plegado de torsion.  Si (*) vale, no es una composicion a posteriori: es un CUADRADO
# CONMUTATIVO entre las dos factorizaciones naturales, que es la respuesta fuerte al "operador que ve
# las dos mitades" que el companion dejo abierto.
#
# EL LADO IZQUIERDO TIENE FORMA CERRADA, y es lo que hace el test barato:
#       s_lambda(X u Y) = sum_nu s_{lambda/nu}(X) s_nu(Y)   =>   eps_nu = s_{lambda/nu}(mu_t),
# o sea un SCHUR SESGADO evaluado en la orbita entera de raices.  No hay que adivinarlo.
#
# CONTROLES
#   C0  FATAL.  Cada lado, por separado, tiene que reconstruir Phi_{t,r} monomio a monomio.  Si uno
#       de los dos no reconstruye, la igualdad de los dos no dice nada.
#   C1  el enunciado del companion:  eps_nu en {0,+-1}.  Es SUYO, no nuestro, y se verifica aqui.
#   C2  la igualdad (*) coeficiente a coeficiente, no solo la suma total.
#   C3  SEÑUELO.  El lado izquierdo con los eps del t EQUIVOCADO.  Tiene que DISCREPAR; si coincide,
#       el test no distingue nada.
#   C4  no vacuidad: n impreso siempre, y el numero de nu y de mu de cada forma.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage commuting_square.sage

import itertools, json, sys
from collections import defaultdict

Sym = SymmetricFunctions(QQ)
sch = Sym.schur()


def phi_bialternante(beta, tt, nvar):
    N = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    L = LaurentPolynomialRing(K, nvar, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(N - 1, -1, -1))
    def alt(expo):
        return matrix(L, N, N, lambda i, j: x[i] ** expo[j]).determinant()
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


def pelar_sp(P, rr, tope=4000):
    """P = sum_mu c_mu sp_mu(z).  Devuelve (coefs, resto)."""
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


def pelar_branching(P, m, rr, tope=9000):
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
        a, b = sp_char(tuple(top[:m]), m), sp_char(tuple(top[m:]), rr)
        for e1, c1 in a.items():
            for e2, c2 in b.items():
                k = e1 + e2
                nv = P.get(k, 0) - B * c1 * c2
                if nv == 0:
                    P.pop(k, None)
                else:
                    P[k] = nv
    return out, P


def tau(eta, tt, mm):
    a = [eta[j] + (mm - (j + 1) + 1) for j in range(mm)]
    cl, sg = [], 1
    for v in a:
        c = v % tt
        if c == 0 or 2 * c == tt:
            return 0
        if c <= mm:
            cl.append(c)
        else:
            cl.append(tt - c); sg *= -1
    if len(set(cl)) != mm:
        return 0
    perm = [mm - cl[j] for j in range(mm)]
    inv = sum(1 for i in range(mm) for j in range(i + 1, mm) if perm[i] > perm[j])
    return sg * (-1) ** mm if False else sg * (-1) ** inv


# ------------------------------------------------------- el lado IZQUIERDO ----------------------
_EPS = {}
def eps_nu(lam, nu, tt):
    """eps_nu = s_{lambda/nu}(mu_t).

    SIMPLIFICACION, y es la que hace el calculo instantaneo Y explica el rango {0,+-1}:
    en la orbita ENTERA de raices t-esimas,  sum_k h_k z^k = prod_j 1/(1 - zeta^j z) = 1/(1 - z^t),
    luego  h_k(mu_t) = [t | k]  y cero si k < 0.  Metiendolo en Jacobi-Trudi para el sesgado,

        eps = det( h_{lam_i - nu_j - i + j} )  =  det( [ t divide lam_i - nu_j - i + j ] ),

    un determinante de una matriz de CEROS Y UNOS.  Que su valor este en {0,+-1} deja de ser un
    hecho medido y pasa a ser una propiedad de ese determinante.  La version anterior expandia el
    Schur sesgado en t variables y sustituia: exacta pero de una fila por hora.
    """
    key = (tuple(lam), tuple(nu), tt)
    if key not in _EPS:
        n = len(lam)
        nn = list(nu) + [0] * (n - len(nu))
        M = matrix(ZZ, n, n, lambda i, j:
                   1 if (lam[i] - nn[j] - i + j) >= 0 and (lam[i] - nn[j] - i + j) % tt == 0 else 0)
        _EPS[key] = M.determinant()
    return _EPS[key]


def s_nu_reciproco(nu, rr):
    """s_nu(z_1,1/z_1,...,z_r,1/z_r) como dict sobre Z^r."""
    L = LaurentPolynomialRing(QQ, rr, 'z')
    zs = L.gens()
    alf = [g ** e for g in zs for e in (1, -1)]
    f = sch[Partition(list(nu))]
    pol = f.expand(2 * rr)
    q = L(pol(*alf))
    return {tuple(e) if hasattr(e, '__iter__') else (e,): c
            for e, c in zip(q.exponents(), q.coefficients()) if c != 0}


t, r = 4, 2
m, R, N = (t - 2) // 2, r + (t - 2) // 2, t + 2 * r
CASOS = [(10, 9, 7, 4, 3, 2, 1, 0), (12, 11, 10, 5, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
         (13, 9, 8, 7, 5, 4, 2, 0), (14, 13, 11, 4, 3, 2, 1, 0)]
DELTA = list(range(N - 1, -1, -1))

print("=" * 122)
print("EL CUADRADO CONMUTATIVO   --   t=%d, r=%d, N=%d" % (t, r, N))
print("=" * 122)
print("")
print("  %-28s | #nu | eps en {0,+-1} | C0 izq | C0 der | (*) coef a coef | señuelo t'" % "beta")
print("  " + "-" * 116)
res = []
for b in CASOS:
    lam = tuple(b[i] - DELTA[i] for i in range(N))
    lam = tuple(x for x in lam if x > 0)
    Phi = phi_bialternante(b, t, r)
    if Phi in (None, "NO-POL"):
        print("  %-28s | ---" % str(b)); continue

    # ---------- lado IZQUIERDO: primero la raiz -----------------------------------------------
    NUS = [nu for nu in Partitions(sum(lam)).list() if False]      # se enumeran por contencion
    NUS = []
    for k in range(sum(lam) + 1):
        for nu in Partitions(k, max_length=2 * r).list():
            # la comprobacion de LONGITUD va antes de indexar: al reves, nu mas largo que lambda
            # revienta con IndexError en vez de descartarse.
            if len(nu) <= len(lam) and all(nu[i] <= lam[i] for i in range(len(nu))):
                NUS.append(tuple(nu))
    eps = {}
    for nu in NUS:
        e = eps_nu(lam, nu, t)
        if e != 0:
            eps[nu] = e
    ok_eps = all(v in (1, -1) for v in eps.values())
    A_izq = defaultdict(lambda: 0)
    for nu, e in eps.items():
        S = s_nu_reciproco(nu, r)
        bb, resto = pelar_sp(S, r)
        for mu, c in bb.items():
            A_izq[mu] += e * c
    A_izq = {mu: c for mu, c in A_izq.items() if c != 0}

    # ---------- lado DERECHO: primero ramificar ------------------------------------------------
    Psi = phi_bialternante(b, 2, R)
    B, _ = pelar_branching(Psi, m, r)
    A_der = defaultdict(lambda: 0)
    for (eta, mu), bb in B.items():
        v = tau(eta, t, m)
        if v:
            A_der[mu] += bb * v
    A_der = {mu: c for mu, c in A_der.items() if c != 0}

    def recon(A):
        out = {}
        for mu, c in A.items():
            for k, v in sp_char(mu, r).items():
                out[k] = out.get(k, 0) + c * v
        return {k: v for k, v in out.items() if v != 0}
    P = {k: v for k, v in Phi.items()}
    c0i = (recon(A_izq) == P)
    c0d = (recon(A_der) == P)
    igual = ({k: int(v) for k, v in A_izq.items()} == {k: int(v) for k, v in A_der.items()})

    # ---------- C3 señuelo: los eps del t equivocado --------------------------------------------
    epsf = {}
    for nu in NUS:
        e = eps_nu(lam, nu, t + 2)
        if e != 0:
            epsf[nu] = e
    A_f = defaultdict(lambda: 0)
    for nu, e in epsf.items():
        S = s_nu_reciproco(nu, r)
        bb, _ = pelar_sp(S, r)
        for mu, c in bb.items():
            A_f[mu] += e * c
    A_f = {mu: c for mu, c in A_f.items() if c != 0}
    dif = (A_f != A_der)

    print("  %-28s | %3d | %-14s | %-6s | %-6s | %-15s | %s"
          % (str(b), len(eps), "si" if ok_eps else "*** NO ***",
             "ok" if c0i else "FALLA", "ok" if c0d else "FALLA",
             "IGUAL" if igual else "*** DISTINTO ***",
             "discrepa" if dif else "*** COINCIDE: no mide ***"))
    sys.stdout.flush()
    res.append({"beta": [int(x) for x in b], "lambda": [int(x) for x in lam],
                "n_eps": int(len(eps)), "eps_pm1": bool(ok_eps),
                "C0_izq": bool(c0i), "C0_der": bool(c0d), "igual": bool(igual),
                "senuelo_discrepa": bool(dif)})

print("")
n = len(res)
print("  n = %d formas" % n)
print("  C1  eps_nu en {0,+-1} (enunciado del companion) : %d/%d" % (sum(x["eps_pm1"] for x in res), n))
print("  C0  cada lado reconstruye Phi                   : izq %d/%d, der %d/%d"
      % (sum(x["C0_izq"] for x in res), n, sum(x["C0_der"] for x in res), n))
print("  C2  (*) coeficiente a coeficiente               : %d/%d" % (sum(x["igual"] for x in res), n))
print("  C3  señuelo con los eps del t equivocado        : discrepa en %d/%d"
      % (sum(x["senuelo_discrepa"] for x in res), n))
print("")
if n and all(x["igual"] and x["C0_izq"] and x["C0_der"] for x in res):
    print("  -> EL CUADRADO CONMUTA en toda la poblacion medida.")
else:
    print("  -> NO conmuta en alguna forma: mirar esa fila antes de decir nada.")
json.dump(res, open("commuting_square_DUMP.json", "w"), indent=1)
print("")
print("=" * 122)
print("DONE")
