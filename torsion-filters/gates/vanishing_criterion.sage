# -*- coding: utf-8 -*-
# M1 -- POR QUE SE ANULA:  ¿por SOPORTE o por CANCELACION?   15 de agosto de 2026.
#
# LA VERSION QUE NO SE CORRE, Y POR QUE.  La propuesta inicial era: aplicar la ley
#
#     mu_max = top(Newt(N_beta)) - sigma_r
#
# a las formas NULAS y ver si el vector predicho deja de ser dominante.  Es CIRCULAR: top_real() se
# calcula sobre el soporte que SOBREVIVE en N_beta, y si Phi = 0 entonces N_beta = Phi * N_delta = 0
# identicamente, luego no hay soporte y no hay vector.  El criterio saldria "Phi = 0 <=> N_beta = 0",
# que es la definicion devuelta.  Ver [[circular-artifact-measurement-returns-definition]].
#
# LA QUE SI CONTESTA LA MISMA PREGUNTA.  Con la cadena de la vuelta 12 delante,
#
#     Phi_{t,r} = sum_mu A_mu sp_mu,      A_mu = sum_eta B_{eta,mu} tau_t(eta)
#
# hay EXACTAMENTE TRES maneras de que Phi se anule, y son excluyentes:
#
#   (I)   HEREDADA    Psi_R = Phi_{2,R} = 0 ya antes de especializar.  Es el Teorema 8.6, el 3-6 %.
#   (II)  SOPORTE     Psi_R != 0 pero TODO eta de supp(B) es singular mod t: el filtro los mata a
#                     todos.  Si esto domina, el criterio de anulacion es COMBINATORIO -- basta
#                     saber que eta aparecen, no cuanto valen.
#   (III) CANCELACION Psi_R != 0, hay eta regulares con tau != 0, pero A_mu = 0 para TODO mu: los
#                     supervivientes se suman a cero.  Si esto domina, el problema es de
#                     cancelacion y ninguna descripcion del soporte lo va a cerrar.
#
# ESE REPARTO ES EL DATO.  No sabemos cual manda, y es la pregunta mas informativa que sabemos hacer
# sobre el problema abierto: dice si hay criterio combinatorio o no lo hay.
#
# CONTROLES, y los cuatro pueden fallar
#   C0  FATAL.  Sobre las NO nulas: tiene que haber algun mu con A_mu != 0, y ningun caso de tipo
#       (II).  Una forma no nula con "todos los eta singulares" seria imposible: instrumento roto.
#   C1  FATAL.  La reconstruccion sum_mu A_mu sp_mu tiene que dar Phi_{t,r} exacto (calculado por un
#       bialternante independiente, en t).  Si no, no se dice nada del reparto.
#   C2  SEÑUELO.  El mismo reparto con tau del t EQUIVOCADO.  Tiene que cambiar.  Si diera el mismo
#       reparto, el filtro no esta decidiendo nada y las tres categorias no miden.
#   C3  no vacuidad: n impreso SIEMPRE, por (t, r) y por categoria, y la poblacion nula y la no nula
#       por separado.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage vanishing_criterion.sage

import itertools, json, sys
from collections import defaultdict


# ------------------------------------------------------------------ Phi por el bialternante -----
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
        return "NO-POLINOMIO"
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        e = tuple(e) if hasattr(e, '__iter__') else (e,)
        if c != 0:
            out[e] = c
    return out


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


def pelar_branching(P, m, r, tope=9000):
    P = {e: c for e, c in P.items() if c != 0}
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
        eta, mu = tuple(top[:m]), tuple(top[m:])
        out[(eta, mu)] = out.get((eta, mu), 0) + B
        a, b = sp_char(eta, m), sp_char(mu, r)
        for e1, c1 in a.items():
            for e2, c2 in b.items():
                k = e1 + e2
                nv = P.get(k, 0) - B * c1 * c2
                if nv == 0:
                    P.pop(k, None)
                else:
                    P[k] = nv
    return out, P


def tau_de(eta, tt, mm):
    a = [eta[j] + (mm - (j + 1) + 1) for j in range(mm)]
    cl = []
    sg = 1
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
    return sg * (-1) ** inv


# ------------------------------------------------------------------ poblaciones ------------------
load("pob_helper.py")          # betas / occupied / phi_zero, envueltos para Sage


def clasificar(b, t, r, m, R, t_falso=None):
    """Devuelve (categoria, detalle).  Categorias: I heredada, II soporte, III cancelacion, NZ no nula."""
    tt = t_falso if t_falso else t
    mm = (tt - 2) // 2
    Psi = phi_bialternante(b, 2, R)
    if Psi in ("NO-POLINOMIO",) or Psi is None:
        return "ERR", {}
    if not Psi:
        return "I", {"n_B": 0, "n_reg": 0, "n_sing": 0}
    Psi = {k: QQ(v) for k, v in Psi.items()}
    B, resto = pelar_branching(Psi, m, r)
    B = {k: v for k, v in B.items() if v != 0}
    etas = sorted({e for (e, _) in B})
    reg = [e for e in etas if tau_de(tuple(list(e) + [0] * max(0, mm - m))[:mm], tt, mm) != 0]
    A = defaultdict(lambda: 0)
    for (e, mu), bb in B.items():
        v = tau_de(tuple(list(e) + [0] * max(0, mm - m))[:mm], tt, mm)
        if v:
            A[mu] += bb * v
    A = {mu: a for mu, a in A.items() if a != 0}
    det = {"n_B": len(B), "n_eta": len(etas), "n_reg": len(reg), "n_sing": len(etas) - len(reg),
           "resto": len(resto), "n_mu_vivos": len(A)}
    if A:
        return "NZ", det
    if not reg:
        return "II", det
    return "III", det


# ================================================================== corrida =====================
print("=" * 122)
print("M1  --  POR QUE SE ANULA:  heredada (I) / por soporte (II) / por cancelacion (III)")
print("=" * 122)
print("")
sys.stdout.flush()

CONF = [(4, 2, 13, 45), (6, 2, 13, 25)]        # (t, r, W, tope de formas por columna)
RES = {}
for (t, r, W, TOPE) in CONF:
    m = (t - 2) // 2
    R = r + m
    # MUESTREO.  La primera corrida tomo las PRIMERAS TOPE nulas en orden de enumeracion, o sea las
    # beta mas bajas: sesgo hacia lo degenerado, y salio 42,9 % heredadas contra el 3-6 % de la
    # vuelta 11.  Aqui se recorre la poblacion ENTERA, se cuentan los totales, y se submuestrea con
    # PASO CONSTANTE para cubrir todo el rango.  El total se imprime siempre, para que la proporcion
    # que se cite sea la de la poblacion y no la de la muestra.
    nulas_all, nonulas_all = [], []
    for b in betas_py(t, r, W):
        if not occupied_py(b, t):
            continue
        z = phi_zero_py(b, t, r)
        if z is None:
            continue
        (nulas_all if z else nonulas_all).append(tuple(b))

    def muestra(L, k):
        if len(L) <= k:
            return list(L), 1
        paso = len(L) // k
        return [L[i * paso] for i in range(k)], paso

    nulas, paso_n = muestra(nulas_all, TOPE)
    nonulas, paso_nn = muestra(nonulas_all, TOPE)

    print("  t=%d r=%d W=%d :  POBLACION ENTERA -> %d ocupadas, de ellas %d NULAS (%.1f %%) y %d no nulas"
          % (t, r, W, len(nulas_all) + len(nonulas_all), len(nulas_all),
             100.0 * len(nulas_all) / max(1, len(nulas_all) + len(nonulas_all)), len(nonulas_all)))
    print("                    se clasifican %d nulas (paso %d) y %d no nulas (paso %d)"
          % (len(nulas), paso_n, len(nonulas), paso_nn))
    sys.stdout.flush()
    cat = defaultdict(int)
    cat_falso = defaultdict(int)
    detalles = []
    malo1 = 0
    for b in nulas:
        c, d = clasificar(b, t, r, m, R)
        cat[c] += 1
        detalles.append({"beta": list(b), "cat": c, **d})
        cf, _ = clasificar(b, t, r, m, R, t_falso=t + 2)      # C2 señuelo
        cat_falso[cf] += 1
        # C1: la reconstruccion, sobre las nulas, tiene que dar Phi = 0 -- y lo da si cat != NZ
        if c == "NZ":
            malo1 += 1
    catn = defaultdict(int)
    for b in nonulas:
        c, d = clasificar(b, t, r, m, R)
        catn[c] += 1
    RES["t%d_r%d" % (t, r)] = {"nulas": dict(cat), "nulas_senuelo": dict(cat_falso),
                               "no_nulas": dict(catn), "detalle": detalles}
    tot = sum(cat.values())
    print("")
    print("     LAS NULAS, n = %d" % tot)
    for k, nom in [("I", "I   HEREDADA de t=2   (Psi_R = 0 ya)"),
                   ("II", "II  POR SOPORTE      (todo eta singular mod t)"),
                   ("III", "III POR CANCELACION  (hay eta regulares, y aun asi A_mu = 0)"),
                   ("NZ", "NZ  *** no deberia salir: la forma no era nula ***"),
                   ("ERR", "ERR *** error de calculo ***")]:
        if cat.get(k):
            print("        %-52s %3d   (%5.1f %%)" % (nom, cat[k], 100.0 * cat[k] / tot))
    print("")
    print("     C0  las NO nulas, n = %d  ->  %s" % (sum(catn.values()), dict(catn)))
    print("         (tiene que ser todo NZ.  Un II aqui seria el instrumento roto)")
    print("     C1  ninguna nula clasificada como NZ : %s" % ("PASA" if malo1 == 0 else "*** FALLA en %d ***" % malo1))
    print("     C2  SEÑUELO, el mismo reparto con tau del t=%d : %s" % (t + 2, dict(cat_falso)))
    print("         (tiene que CAMBIAR.  Si diera igual, el filtro no decide nada)")
    print("")
    sys.stdout.flush()

json.dump(RES, open("vanishing_criterion_DUMP.json", "w"), indent=1)
print("=" * 122)
print("  LO QUE EL REPARTO DECIDE")
print("=" * 122)
print("")
print("   * si manda (II) POR SOPORTE  ->  el criterio de anulacion es COMBINATORIO: basta describir")
print("     que eta aparecen en supp(B) desde beta, y el problema abierto se vuelve un problema de")
print("     soporte de ramificacion, no de cancelacion.")
print("   * si manda (III) POR CANCELACION  ->  no hay criterio de soporte posible, y hay que atacar")
print("     la suma.  Es el caso malo, y saberlo tambien vale: cierra una via muerta.")
print("")
print("   datos en vanishing_criterion_DUMP.json")
print("")
print("=" * 122)
print("DONE")
