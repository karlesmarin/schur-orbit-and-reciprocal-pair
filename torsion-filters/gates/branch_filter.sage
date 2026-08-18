# -*- coding: utf-8 -*-
# LA RUTA DEL BRANCHING FILTRADO.   Phi_{t,r} = sum_mu ( sum_eta B_{eta,mu} tau_t(eta) ) sp_mu.
# 15 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 12, P12.  El propone factorizar en dos pasos:
#
#   (1)  Psi_R := Phi_{2,R}  =  sum_{eta,mu}  B_{eta,mu} sp_eta(y) sp_mu(z)     [branching Sp_{2m} x Sp_{2r}]
#   (2)  Phi_{t,r} = Psi_R |_{y = (xi,...,xi^m)}  =  sum_mu ( sum_eta B_{eta,mu} tau_t(eta) ) sp_mu(z)
#   (3)  A_mu = sum_eta B_{eta,mu} tau_t(eta),   tau_t(eta) en {0,+-1}  por (T), ya medida en torsion_filter.
#
# Dado (1), el paso (2) es formalmente una tautologia -- PERO (1) no lo es, y la cadena entera
# tampoco: hay que ver que nuestro Phi_{2,R} SE DEJA escribir asi con B enteros, y que el resultado
# coincide monomio a monomio con el Phi_{t,r} que calculamos por otro camino.  Eso es lo que se mide.
#
# Y de paso sale la columna que el pide para P13:  que eta acompaña al mu maximal, y que PARED lo mata.
#
# CONTROLES
#   C0  FATAL.  Reconstruccion del branching: sum B_{eta,mu} sp_eta sp_mu tiene que devolver Psi_R
#       EXACTO.  Resto no nulo => la descomposicion no existe o el pelado no termino.
#   C1  FATAL.  B entero.  Si sale un B racional, (1) es falsa tal como el la escribe.
#   C2  FATAL.  La ruta de EL contra el Phi_{t,r} directo (bialternante en t), monomio a monomio.
#   C3  SEÑUELO.  La misma ruta con tau' = "todo eta sobrevive con +1" (o sea, ignorando el filtro).
#       Tiene que DISCREPAR de Phi_{t,r}.  Si coincidiera, el filtro no estaria haciendo nada y C2 no
#       mediria nada.
#   C4  SEÑUELO.  La ruta con tau'' = tau del t EQUIVOCADO (t+2 en vez de t).  Tiene que discrepar.
#   C5  no vacuidad: n impreso siempre, y las formas que se anulan aparte.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage branch_filter.sage

import itertools, json

# ------------------------------------------------------------------ Phi por el bialternante -----
def phi_bialternante(beta, t, nvar):
    """Phi_{t,nvar} = det(x_i^{beta_j})/det(x_i^{delta_j}),  x = (raices t-esimas, z_k, 1/z_k)."""
    N = t + 2 * nvar
    K = CyclotomicField(t) if t > 2 else QQ
    zeta = K.gen() if t > 2 else K(-1)
    L = LaurentPolynomialRing(K, nvar, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(t)] + [g ** e for g in zs for e in (1, -1)]
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
        return "NO-POLINOMIO"
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        e = tuple(e) if hasattr(e, '__iter__') else (e,)
        if c != 0:
            out[e] = c
    return out


# ------------------------------------------------------------------ caracteres sp ----------------
_SP = {}
def sp_char(mu, r):
    key = (tuple(mu), r)
    if key not in _SP:
        W = WeylCharacterRing("C%d" % r)
        el = W(W.space().from_vector(vector(list(mu))))
        d = {}
        for wt, mult in el.weight_multiplicities().items():
            k = tuple(int(v) for v in wt.to_vector())
            d[k] = d.get(k, 0) + mult
        _SP[key] = d
    return _SP[key]


def sp_producto(eta, mu, m, r):
    """sp_eta(y) * sp_mu(z) como dict sobre Z^{m+r}."""
    a, b = sp_char(eta, m), sp_char(mu, r)
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            k = e1 + e2
            out[k] = out.get(k, 0) + c1 * c2
    return out


def pelar_branching(P, m, r, tope=4000):
    """P = sum B_{eta,mu} sp_eta(y) sp_mu(z).  Devuelve (dict {(eta,mu): B}, resto)."""
    P = {e: c for e, c in P.items() if c != 0}
    out = {}
    for _ in range(tope):
        if not P:
            return out, {}
        dom = [e for e in P
               if list(e[:m]) == sorted(e[:m], reverse=True) and (m == 0 or min(e[:m]) >= 0)
               and list(e[m:]) == sorted(e[m:], reverse=True) and min(e[m:]) >= 0]
        if not dom:
            return out, P
        top = max(dom, key=lambda e: (sum(e), e))
        B = P[top]
        eta, mu = tuple(top[:m]), tuple(top[m:])
        out[(eta, mu)] = out.get((eta, mu), 0) + B
        for k, v in sp_producto(eta, mu, m, r).items():
            nv = P.get(k, 0) - B * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P


# ------------------------------------------------------------------ el filtro (T) ----------------
def clase(a, t, m):
    c = a % t
    if c == 0 or 2 * c == t:
        return None
    return (c, +1) if c <= m else (t - c, -1)

def tau(eta, t, m):
    """(valor, tipo de pared).  Por la regla (T)+(T'), ya validada contra dos rutas en torsion_filter."""
    a = [eta[j] + (m - (j + 1) + 1) for j in range(m)]
    cs = [clase(v, t, m) for v in a]
    if any(c is None for c in cs):
        return 0, "pared a_i=0 o t/2"
    cl = [c for c, _ in cs]
    if len(set(cl)) != m:
        return 0, "pared a_i=+-a_j"
    s = 1
    for _, sg in cs:
        s *= sg
    perm = [m - cl[j] for j in range(m)]
    inv = sum(1 for i in range(m) for j in range(i + 1, m) if perm[i] > perm[j])
    return s * (-1) ** inv, "-"


def ruta_de_el(B, t, m, r, filtro):
    """A_mu = sum_eta B_{eta,mu} filtro(eta).  Devuelve dict {mu: A_mu} sin ceros."""
    A = {}
    for (eta, mu), b in B.items():
        v = filtro(eta)
        if v:
            A[mu] = A.get(mu, 0) + b * v
    return {mu: a for mu, a in A.items() if a != 0}


def a_monomios(A, r):
    out = {}
    for mu, a in A.items():
        for k, v in sp_char(mu, r).items():
            out[k] = out.get(k, 0) + a * v
    return {k: v for k, v in out.items() if v != 0}


def dominantes_maximales(A):
    """los mu maximales en el orden de dominancia dentro del soporte."""
    S = list(A)
    return [mu for mu in S
            if not any(nu != mu and all(sum(nu[:k + 1]) >= sum(mu[:k + 1]) for k in range(len(mu)))
                       for nu in S)]


# ================================================================== CASOS =======================
CASOS = {
    4: [(18, 17, 11, 8, 7, 6, 1, 0), (10, 9, 7, 4, 3, 2, 1, 0), (14, 13, 11, 4, 3, 2, 1, 0),
        (12, 11, 10, 5, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
        (13, 9, 8, 7, 5, 4, 2, 0), (19, 17, 11, 8, 7, 6, 1, -1), (21, 17, 11, 8, 7, 6, 1, -3)],
    6: [(13, 11, 9, 7, 5, 4, 3, 2, 1, 0), (15, 13, 11, 8, 6, 5, 3, 2, 1, 0),
        (17, 13, 11, 9, 7, 5, 3, 2, 1, 0), (12, 11, 10, 9, 8, 5, 3, 2, 1, 0)],
}
r = 2
FILAS = []
for t in (4, 6):
    m = (t - 2) // 2
    R = r + m
    print("")
    print("=" * 118)
    print("t = %d   m = %d   R = %d   ->   branching  Sp_%d  >  Sp_%d x Sp_%d" % (t, m, R, 2 * R, 2 * m, 2 * r))
    print("=" * 118)
    print("")
    print("  beta                            | #B  | C0 resto | C1 enteros | C2 su ruta == Phi | C3 señuelo sin filtro | C4 señuelo t' ")
    print("  " + "-" * 148)
    for b in CASOS[t]:
        Psi = phi_bialternante(b, 2, R)          # el objeto de t=2, en R pares
        Phi = phi_bialternante(b, t, r)          # el objeto real, camino independiente
        if Psi in (None, "NO-POLINOMIO") or Phi in (None, "NO-POLINOMIO"):
            print("  %-31s | ---  denominador nulo / no polinomio" % str(b))
            continue
        Psi = {k: QQ(v) for k, v in Psi.items()}
        B, resto = pelar_branching(Psi, m, r)
        B = {k: v for k, v in B.items() if v != 0}
        enteros = all(QQ(v).denominator() == 1 for v in B.values())
        A = ruta_de_el(B, t, m, r, lambda e: tau(e, t, m)[0])
        suyo = a_monomios(A, r)
        real = {k: v for k, v in Phi.items()}
        ok2 = (suyo == real)
        # señuelo C3: sin filtro (todo eta cuenta con +1)
        s3 = a_monomios(ruta_de_el(B, t, m, r, lambda e: 1), r)
        # señuelo C4: el tau del t equivocado
        tp = t + 2
        mp = (tp - 2) // 2
        s4 = a_monomios(ruta_de_el(B, t, m, r, lambda e: tau(tuple(list(e) + [0] * (mp - m))[:mp], tp, mp)[0]), r)
        FILAS.append((t, b, B, A, real))
        print("  %-31s | %3d | %8s | %10s | %17s | %21s | %s"
              % (str(b), len(B), "0" if not resto else "*** %d ***" % len(resto),
                 "OK" if enteros else "*** NO ***",
                 "OK" if ok2 else "*** DISTINTO ***",
                 "discrepa (ok)" if s3 != real else "*** COINCIDE: no mide ***",
                 "discrepa (ok)" if s4 != real else "*** COINCIDE: no mide ***"))

# ================================================================== P13 =========================
print("")
print("=" * 118)
print("P13 -- LA COLUMNA QUE EL PIDE:  que eta acompaña al mu maximal, y que pared mata a los de arriba")
print("=" * 118)
print("")
print("  t | beta                            | mu_max     | A    | eta que sobreviven en mu_max        | eta MAS ALTO del branching / su pared")
print("  " + "-" * 148)
for (t, b, B, A, real) in FILAS:
    m = (t - 2) // 2
    if not A:
        print("  %d | %-31s | Phi == 0" % (t, str(b)))
        continue
    maxi = dominantes_maximales(A)
    for mu in maxi:
        vivos = [(eta, B[(eta, mu)], tau(eta, t, m)[0])
                 for (eta, nu) in B if nu == mu and tau(eta, t, m)[0] != 0]
        todos = sorted([eta for (eta, nu) in B if nu == mu], key=lambda e: (sum(e), e), reverse=True)
        alto = todos[0] if todos else None
        pared = tau(alto, t, m)[1] if alto else "-"
        print("  %d | %-31s | %-10s | %-4s | %-35s | %s  [%s]"
              % (t, str(b), str(mu), str(A[mu]), str(vivos)[:35], str(alto), pared))

print("")
print("  n = %d formas medidas en total (%d con t=4, %d con t=6)"
      % (len(FILAS), sum(1 for f in FILAS if f[0] == 4), sum(1 for f in FILAS if f[0] == 6)))
print("")
print("=" * 118)
print("DONE")
