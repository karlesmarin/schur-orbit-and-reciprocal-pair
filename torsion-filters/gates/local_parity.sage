# -*- coding: utf-8 -*-
# SU LEMA (10):  los r_i SIN reordenar, y la anatomia de la cancelacion.  15 de agosto de 2026.
#
# DOS COSAS, Y LA PRIMERA ES UNA COMPUERTA.
#
# (A) EL LEMA (10), QUE EL DERIVA Y NO CITA.  Afirma que bajo el doble entrelazado mu << Lambda no
#     hace falta el rearrangement global de Yacobi, porque
#
#         x_i = min(Lambda_i, mu_{i-1}),      y_i = max(Lambda_{i+1}, mu_i),
#         con  mu_0 = +infinito,  mu_R = 0,  Lambda_{R+1} = 0,
#         luego  r_i = min(Lambda_i, mu_{i-1}) - max(Lambda_{i+1}, mu_i).
#
#     Si vale, la supervivencia pasa de "ordenar 2R numeros y mirar paridades" a R TESTS LOCALES, y
#     eso es lo que hace atacable una prueba combinatoria.  Como es SUYO y no de Yacobi, se verifica
#     contra el rearrangement de verdad antes de usarlo.  Es exactamente lo que la regla pide.
#
# (B) LA ANATOMIA.  Sus tres escenarios para la suma que da +-1, sobre los c_Lambda = a_Lambda * eps:
#         A  cancelacion por magnitud   {7,-7,3,-3,1}  -> hay involucion al nivel Lambda
#         B  cancelacion tras atomizar  {5,-3,-1}      -> no hay pairing de Lambda, quiza de objetos
#         C  telescopaje                q1-q2, q2-q3.. -> colapsa a un borde
#     Y su observacion de conteo, que mata nuestra P18 tal como estaba: con 10 supervivientes, una
#     involucion con parejas y UN punto fijo necesitaria 10-1 = 9 par.  Cierto.
#
# CONTROLES
#   C0  FATAL.  (10) contra el rearrangement de Yacobi, r_i a r_i, en TODOS los pares (Lambda,mu)
#       con mu << Lambda.  Un solo desacuerdo y no se usa.
#   C1  se prueba tambien en pares que NO entrelazan, para ver si (10) falla ahi -- su enunciado lo
#       condiciona al entrelazado, y hay que comprobar que la condicion es necesaria y no adorno.
#   C2  la anatomia se imprime ENTERA para cada forma: los c_Lambda uno a uno, no un resumen.
#   C3  no vacuidad: n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage local_parity.sage

import itertools, json, sys
from collections import defaultdict

t, r = 4, 2
m, R, N = 1, 3, 8
INF = 10 ** 9


def yacobi_sort(mu, Lam):
    """los r_i de Yacobi: reordenacion decreciente de (mu, Lambda^+)."""
    n = len(Lam)
    seq = sorted(list(mu) + list(Lam) + [0], reverse=True)
    return [seq[2 * i] - seq[2 * i + 1] for i in range(n)]


def local_r(mu, Lam):
    """SU (10): r_i = min(Lam_i, mu_{i-1}) - max(Lam_{i+1}, mu_i), sin ordenar nada."""
    n = len(Lam)
    L = list(Lam) + [0]
    M = [INF] + list(mu) + [0] * (n + 2)
    return [min(L[i], M[i]) - max(L[i + 1], M[i + 1]) for i in range(n)]


def entrelaza(mu, Lam):
    n = len(Lam)
    Lp = list(Lam) + [0]
    return all(Lp[i] >= mu[i] >= Lp[i + 2] for i in range(n - 1))


# ================================================================== (A) la compuerta ============
print("=" * 122)
print("C0/C1  SU LEMA (10) CONTRA EL REARRANGEMENT DE YACOBI")
print("=" * 122)
print("")
TOPE = 9
ok = mal = n_ent = n_no = mal_no = 0
ejemplos = []
for L1 in range(TOPE + 1):
    for L2 in range(L1 + 1):
        for L3 in range(L2 + 1):
            Lam = (L1, L2, L3)
            for u1 in range(TOPE + 1):
                for u2 in range(u1 + 1):
                    mu = (u1, u2)
                    a, b = yacobi_sort(mu, Lam), local_r(mu, Lam)
                    if entrelaza(mu, Lam):
                        n_ent += 1
                        if a == b:
                            ok += 1
                        else:
                            mal += 1
                            if len(ejemplos) < 3:
                                ejemplos.append((Lam, mu, a, b))
                    else:
                        n_no += 1
                        mal_no += (a != b)
print("   pares (Lambda,mu) con mu << Lambda        : %d" % n_ent)
print("   (10) coincide con el rearrangement        : %d   discrepa: %d" % (ok, mal))
print("   pares SIN entrelazado                     : %d   de los que (10) discrepa: %d (%.1f %%)"
      % (n_no, mal_no, 100.0 * mal_no / max(1, n_no)))
print("")
if mal == 0:
    print("   -> SU LEMA (10) VALE en todo el rango, y la condicion de entrelazado NO es adorno:")
    print("      fuera de ella la formula local se separa del rearrangement en el %.0f %% de los casos."
          % (100.0 * mal_no / max(1, n_no)))
else:
    print("   -> (10) FALLA.  Ejemplos:")
    for Lam, mu, a, b in ejemplos:
        print("      Lam=%s mu=%s   Yacobi=%s   suyo=%s" % (str(Lam), str(mu), str(a), str(b)))
sys.stdout.flush()

# ================================================================== (B) la anatomia =============
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


def eps_de(rs):
    if any(x % 2 for x in rs):
        return 0
    return (-1) ** (sum(rs) // 2)


CASOS = [(10, 9, 7, 4, 3, 2, 1, 0), (12, 11, 10, 5, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
         (13, 9, 8, 7, 5, 4, 2, 0), (14, 13, 11, 4, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0),
         (18, 17, 11, 8, 7, 6, 1, 0)]
print("")
print("=" * 122)
print("C2  LA ANATOMIA DE LA CANCELACION:  los c_Lambda = a_Lambda * eps  en mu_max, uno a uno")
print("=" * 122)
RES = []
for b in CASOS:
    Psi = phi_bialternante(b, 2, R)
    if Psi in (None, "NO-POL"):
        continue
    aL, _ = pelar(Psi, R)
    aL = {k: int(v) for k, v in aL.items() if v != 0}
    # mu_max por la formula (13)
    MUS = set()
    for Lam in aL:
        for u1 in range(Lam[0] + 1):
            for u2 in range(u1 + 1):
                if entrelaza((u1, u2), Lam):
                    MUS.add((u1, u2))
    A = {}
    cs = {}
    for mu in MUS:
        tot, lista = 0, []
        for Lam, a in aL.items():
            if not entrelaza(mu, Lam):
                continue
            e = eps_de(local_r(mu, Lam))
            if e:
                tot += a * e
                lista.append((Lam, a * e))
        if tot:
            A[mu] = tot
            cs[mu] = lista
    S = list(A)
    maxi = [mu for mu in S if not any(nu != mu and all(sum(nu[:k + 1]) >= sum(mu[:k + 1])
                                                       for k in range(r)) for nu in S)]
    if len(maxi) != 1:
        print("   %-28s | *** mu_max no unico: %s ***" % (str(b), str(maxi))); continue
    mm = maxi[0]
    valores = sorted((c for _, c in cs[mm]), key=lambda x: (-abs(x), x))
    pos = sorted([v for v in valores if v > 0], reverse=True)
    neg = sorted([v for v in valores if v < 0])
    empareja = (sorted(abs(v) for v in pos) == sorted(abs(v) for v in neg)) if len(pos) != len(neg) else None
    print("")
    print("   %-28s  mu_max=%s   A=%d   %d supervivientes" % (str(b), str(mm), A[mm], len(valores)))
    print("      c_Lambda: %s" % str(valores)[:110])
    print("      positivos %s   negativos %s   suma %d" % (str(pos)[:44], str(neg)[:44], sum(valores)))
    RES.append({"beta": [int(x) for x in b], "mu_max": [int(x) for x in mm], "A": int(A[mm]),
                "c": [int(v) for v in valores]})
    sys.stdout.flush()

print("")
print("  n = %d formas" % len(RES))
print("  su argumento de conteo: una involucion con parejas y UN punto fijo exige |S| impar.")
for x in RES:
    par = "impar" if len(x["c"]) % 2 else "PAR -> imposible"
    print("     %-28s |S| = %2d  (%s)" % (str(tuple(x["beta"])), len(x["c"]), par))
json.dump(RES, open("local_parity_DUMP.json", "w"), indent=1)
print("")
print("=" * 122)
print("DONE")
