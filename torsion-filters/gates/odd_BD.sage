# -*- coding: utf-8 -*-
# EL IMPAR EN SU GRUPO:  B_{R'} -> B_{m'} x D_r,  la propuesta 3 de su reseña.   16 de agosto de 2026.
#
# DE DONDE SALE.  Su reseña (vuelta 14, punto 1) dice que el pipeline de §4 no cubre el impar y que la
# lectura correcta no es "el par con una pared menos" sino OTRO symmetric pair:
#
#       t par:    C_R  -> C_m  x C_r          (lo que llevamos hecho)
#       t impar:  B_R' -> B_m' x D_r          (esto, y sin orbit algebra)
#
# Y el reparto del alfabeto lo decide donde cae el punto fijo.  Para t impar,
#
#       mu_t = {1} u {zeta^{+-j}}_{j=1..m'},      m' = (t-1)/2,
#
# hay UN solo punto fijo de x -> 1/x, y va con la TORSION:  (1, zeta^{+-1..m'}) son 2m'+1 = t
# numeros, o sea el toro de SO_{2m'+1} = B_{m'}, y las r parejas libres z^{+-} son el toro de
# SO_{2r} = D_r.  Juntos, 2R'+1 = t+2r = N numeros: el toro de B_{R'}, R' = m'+r.
#
# ESO CORRIGE LO QUE HICIMOS EN odd_companion.sage, que puso la torsion en C_{m'} -- es decir, dejo
# el 1 fuera del bloque de torsion.  Aqui el 1 entra, y el filtro pasa a ser un caracter ORTOGONAL
# IMPAR evaluado en un elemento cuyos autovalores son TODAS las raices t-esimas de la unidad.
#
# LA CADENA, entera:
#
#   Phi_{t,r}  =  s_lambda(1, zeta, ..., zeta^{t-1}, z^{+-})
#              =  Phi_{1,R'} evaluado en  w = (zeta^1..zeta^{m'}, z_1..z_r)
#              =  sum_Lambda a^B_Lambda  o^{B_R'}_Lambda(w)                    [Littlewood impar]
#              =  sum_{Lambda,eta,mu} a^B_Lambda B_{Lambda;eta,mu} o^B_eta(y) o^D_mu(z)   [ramificacion]
#              =  sum_mu A^D_mu o^D_mu(z),      A^D_mu = sum_{Lambda,eta} a^B_Lambda B_{...} tau^B(eta)
#
# y la diferencia conceptual con el par es que aqui NO hay caracter virtual ni plegado: la
# restriccion GL_N -> SO_N es honesta y a^B_Lambda son multiplicidades de verdad.
#
# LO QUE SE MIDE
#   D1  el filtro impar  tau^B_t(eta) = o^{B_m'}_eta(1, zeta^{+-1}, ..., zeta^{+-m'}).
#       ¿cae en {0,+-1}?  ¿cuantos eta sobreviven?
#   D2  la regla (T^B) PREDICHA antes de correr:  con A_j = 2 eta_j + 2(m'-j) + 1  (los 2a_j de
#       Weyl, enteros impares),  tau != 0  <=>  los A_j son no nulos mod t y ocupan una vez cada
#       clase de {+-1,...,+-m'} mod t.  Dos especies de pared, no tres: la pared A = t/2 NO EXISTE
#       porque t es impar.  Se mide en LAS DOS direcciones por separado.
#   D3  (H) en el impar: ¿hay un unico mu maximal y vale A^D_mu = +-1?  Es la mitad que convertiria
#       la observacion en un "parity dichotomy theorem".  OJO: se compara sobre mu^+ = (mu_1,..,|mu_r|)
#       -- el peso de O(2r) --, porque en D_r cada mu con mu_r != 0 viene con su quiral mu* y los dos
#       son incomparables en dominancia; comparar sin plegar deja el conjunto de maximales VACIO.
#       (Ese fue el fallo de la primera corrida: 0 maximales en 8 de 9 formas.)
#   D4  la quiralidad: Phi es invariante bajo z_i <-> 1/z_i por separado (O(2r), no solo SO(2r)),
#       luego mu y mu* = (mu_1,...,-mu_r) tienen que salir con el MISMO coeficiente.  Si no, el
#       peldaño D esta mal montado.
#
# CONTROLES, y todos pueden fallar
#   C0  FATAL.  ruta 1 (bialternante en z, pelado en D_r) == ruta 2 (la cadena entera), coeficiente
#       a coeficiente, en todas las mu.  Si esto falla, el reparto B x D es el equivocado.
#   C1  el pelado en B_{R'} tiene que dar resto 0 y coeficientes ENTEROS.  n impreso siempre.
#   C2  SEÑUELO t': el mismo circuito con tau calculado en una raiz (t+2)-esima.  Tiene que
#       DISCREPAR; si empata, el filtro no esta decidiendo nada (es el empate de la P15 en el par).
#   C3  SEÑUELO sin filtro: tau == 1 para todo eta.  Tiene que DISCREPAR -- si no, la torsion no
#       hace nada y todo esto es Littlewood.
#   C4  SEÑUELO del reparto: la torsion en C_{m'} (lo de odd_companion, el 1 fuera).  Tiene que
#       DISCREPAR, y esa discrepancia ES la correccion suya medida.
#   C5  no vacuidad: n de Lambda, de eta y de mu impreso en cada forma.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage odd_BD.sage

import json
import sys
from collections import defaultdict

# ---------------------------------------------------------------- el objeto, por bialternante

def phi_bialt(beta, tt, nvar):
    """s_beta(zeta^0..zeta^{tt-1}, z_1^{+-},...,z_nvar^{+-}) como dict {exponente: coef}."""
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
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        k = tuple(int(v) for v in e) if hasattr(e, '__iter__') else (int(e),)
        cc = QQ(c) if c in QQ else None
        if cc is None:
            return "NO-RAC"
        if cc != 0:
            out[k] = cc
    return out


# ---------------------------------------------------------------- caracteres y pelado

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


def dominante(e, typ):
    if typ == "B":
        return list(e) == sorted(e, reverse=True) and min(e) >= 0
    # D: e_1 >= ... >= e_{r-1} >= |e_r|
    f = list(e)
    return all(f[i] >= f[i + 1] for i in range(len(f) - 2)) and f[-2] >= abs(f[-1]) if len(f) >= 2 \
        else f[0] >= 0


def pelar(P, typ, rk, tope=6000):
    """descompone el polinomio P en caracteres de tipo typ y rango rk.  Devuelve (coefs, resto)."""
    P = {e: QQ(c) for e, c in P.items() if c != 0}
    out = {}
    for _ in range(tope):
        if not P:
            return out, {}
        dom = [e for e in P if dominante(e, typ)]
        if not dom:
            return out, P
        mu = max(dom, key=lambda e: (sum(abs(v) for v in e), e))
        c = P[mu]
        out[mu] = out.get(mu, 0) + c
        for k, v in car(typ, rk, mu).items():
            nv = P.get(k, 0) - c * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return out, P


# ---------------------------------------------------------------- el filtro impar

def tau_B(mp, eta, tt, salto=1):
    """o^{B_mp}_eta evaluado en el elemento con autovalores 1, w^{+-1}, ..., w^{+-mp},
       w = raiz tt-esima primitiva.  salto=1 es el nuestro; salto!=1 es señuelo."""
    K = CyclotomicField(tt)
    z = K.gen()
    s = K(0)
    for wt, mult in car("B", mp, eta).items():
        s += mult * z ** (sum(salto * (i + 1) * wt[i] for i in range(mp)) % tt)
    return s


def tau_C(mp, eta, tt):
    """SEÑUELO C4: el filtro de odd_companion -- simplectico C_mp, el 1 fuera de la torsion."""
    K = CyclotomicField(tt)
    z = K.gen()
    s = K(0)
    for wt, mult in car("C", mp, eta).items():
        s += mult * z ** (sum((i + 1) * wt[i] for i in range(mp)) % tt)
    return s


def regla_TB(eta, mp, tt):
    """(T^B) PREDICHA:  A_j = 2 eta_j + 2(mp-j) + 1 no nulos mod tt y distintos salvo signo."""
    A = [2 * int(eta[j]) + 2 * (mp - j - 1) + 1 for j in range(mp)]
    cl = []
    for a in A:
        c = a % tt
        if c == 0:
            return False
        cl.append(min(c, tt - c))
    return len(set(cl)) == mp


# ---------------------------------------------------------------- la ramificacion

_BR = {}
def branch_BD(Rp, mp, rr, Lam):
    key = (Rp, mp, rr, tuple(int(v) for v in Lam))
    if key not in _BR:
        W = WeylCharacterRing("B%d" % Rp)
        X = WeylCharacterRing("B%dxD%d" % (mp, rr))
        br = branching_rule("B%d" % Rp, "B%dxD%d" % (mp, rr), "orthogonal_sum")
        el = W(W.space().from_vector(vector([Integer(v) for v in Lam])))
        res = el.branch(X, rule=br)
        d = {}
        for wt, c in res.monomial_coefficients().items():
            v = [int(u) for u in wt.to_vector()]
            d[(tuple(v[:mp]), tuple(v[mp:]))] = int(c)
        _BR[key] = d
    return _BR[key]


# ---------------------------------------------------------------- el peso superior, plegado a O(2r)

def top_Omu(A):
    """A es {mu de D_r: coef}.  Devuelve (mu^+ maximal unico o None, coef, n de maximales).
       Se pliega a mu^+ = (mu_1,...,|mu_r|) porque mu y su quiral mu* son incomparables."""
    Ap = {}
    for mu, c in A.items():
        k = tuple(list(mu[:-1]) + [abs(mu[-1])])
        Ap[k] = c            # A(mu) = A(mu*) ya verificado en D4
    S = list(Ap)
    maxi = [m for m in S if not any(
        n != m and all(sum(n[:k + 1]) >= sum(m[:k + 1]) for k in range(len(m))) for n in S)]
    if len(maxi) == 1:
        return maxi[0], Ap[maxi[0]], 1
    return None, None, len(maxi)


def betas_estrictas(N, tope):
    """todas las beta estrictamente decrecientes con min 0 y max <= tope.
       Phi es invariante por traslacion: el producto del alfabeto vale 1 para t impar."""
    out = []
    for c in Combinations(range(1, tope + 1), N - 1):
        out.append(tuple(sorted(c, reverse=True)) + (0,))
    return out


# ---------------------------------------------------------------- las formas

# C0 (la cadena entera, cara) sobre estas:
CASOS = [
    (3, 2, [(9, 7, 5, 3, 2, 1, 0), (10, 9, 6, 4, 2, 1, 0), (11, 8, 6, 4, 3, 1, 0),
            (8, 6, 5, 3, 2, 1, 0), (12, 10, 7, 5, 3, 2, 0)]),
    (5, 2, [(10, 9, 7, 6, 4, 3, 2, 1, 0), (11, 9, 8, 6, 5, 3, 2, 1, 0),
            (12, 10, 8, 7, 5, 4, 2, 1, 0)]),
    (3, 3, [(12, 10, 8, 6, 4, 3, 2, 1, 0)]),
]

# (H) impar (solo ruta 1, barata) sobre la poblacion entera:
POBLACION = [(3, 2, 9), (5, 2, 10)]

RES = []
print("=" * 124)
print("EL IMPAR EN SU GRUPO:  B_{R'} -> B_{m'} x D_r     (propuesta 3 de la reseña)")
print("=" * 124)

for (t, r, betas) in CASOS:
    mp = (t - 1) // 2
    Rp = mp + r
    N = t + 2 * r
    if r < 2:
        continue
    print("")
    print("-" * 124)
    print("  t=%d  (m'=%d)   r=%d   R'=%d   N=%d" % (t, mp, r, Rp, N))
    print("-" * 124)

    # ---- D1/D2: el filtro, una vez por (t, m')
    ETAS = [tuple(e) for k in range(0, 4 * t + 1) for e in Partitions(k, max_length=mp)]
    ETAS = [e + (0,) * (mp - len(e)) for e in ETAS]
    val = {}
    for e in ETAS:
        v = tau_B(mp, e, t)
        vv = QQ(v) if v in QQ else None
        val[e] = vv
    no_rac = [e for e in ETAS if val[e] is None]
    vivos = [e for e in ETAS if val[e] is not None and val[e] != 0]
    fuera = [e for e in vivos if abs(val[e]) != 1]
    pred_ok = sum(1 for e in ETAS if val[e] is not None and (val[e] != 0) == regla_TB(e, mp, t))
    print("  D1  eta probados %d | tau no racional: %d | tau != 0: %d | |tau| != 1 entre los vivos: %d"
          % (len(ETAS), len(no_rac), len(vivos), len(fuera)))
    print("  D2  (T^B) predicha acierta %d/%d   (cero-falsos %d, nocero-falsos %d)"
          % (pred_ok, len(ETAS),
             sum(1 for e in ETAS if val[e] == 0 and regla_TB(e, mp, t)),
             sum(1 for e in ETAS if val[e] is not None and val[e] != 0 and not regla_TB(e, mp, t))))
    sys.stdout.flush()

    print("")
    print("  beta                          |  #Lam | resto B | ent | #mu | C0  | quiral | mu_max              |  A")
    print("  " + "-" * 120)
    for b in betas:
        # ---- ruta 1: la verdad
        P1 = phi_bialt(b, t, r)
        if P1 in (None, "NO-POL", "NO-RAC"):
            print("  %-29s |  ruta 1 falla: %s" % (str(b), P1))
            continue
        A1, rest1 = pelar(P1, "D", r)
        A1 = {k: QQ(v) for k, v in A1.items() if v != 0}

        # ---- ruta 2: la cadena
        P2 = phi_bialt(b, 1, Rp)
        if P2 in (None, "NO-POL", "NO-RAC"):
            print("  %-29s |  ruta 2 falla: %s" % (str(b), P2))
            continue
        aB, restB = pelar(P2, "B", Rp)
        aB = {k: QQ(v) for k, v in aB.items() if v != 0}
        enteros = all(v.denominator() == 1 for v in aB.values())

        A2 = defaultdict(lambda: QQ(0))
        A2s = {"t2": defaultdict(lambda: QQ(0)), "sf": defaultdict(lambda: QQ(0)),
               "C": defaultdict(lambda: QQ(0))}
        for Lam, a in aB.items():
            for (eta, mu), c in branch_BD(Rp, mp, r, Lam).items():
                tv = tau_B(mp, eta, t)
                tv = QQ(tv) if tv in QQ else None
                if tv is None:
                    continue
                A2[mu] += a * c * tv
                # señuelos
                s2 = tau_B(mp, eta, t + 2)
                A2s["t2"][mu] += a * c * (QQ(s2) if s2 in QQ else 0)
                A2s["sf"][mu] += a * c
                sC = tau_C(mp, eta, t)
                A2s["C"][mu] += a * c * (QQ(sC) if sC in QQ else 0)
        A2 = {k: v for k, v in A2.items() if v != 0}

        c0 = (A1 == A2)
        sen = {k: ({kk: vv for kk, vv in d.items() if vv != 0} != A1) for k, d in A2s.items()}

        # D4: quiralidad
        quir = all(A1.get(mu, 0) == A1.get(mu[:-1] + (-mu[-1],), 0) for mu in A1 if mu[-1] != 0)

        # D3: (H) impar
        mm, amm, nmax = top_Omu(A1)
        print("  %-29s | %5d | %7s | %3s | %3d | %-3s | %-6s | %-19s | %s"
              % (str(b), len(aB), "0" if not restB else "NO 0", "si" if enteros else "NO",
                 len(A1), "ok" if c0 else "FALLA", "ok" if quir else "NO",
                 str(mm) if mm else "no unico (%d)" % nmax,
                 str(amm) if mm else "-"))
        sys.stdout.flush()
        RES.append({"t": int(t), "r": int(r), "beta": [int(x) for x in b],
                    "n_Lambda": int(len(aB)), "resto_B_cero": bool(not restB),
                    "aB_enteros": bool(enteros), "n_mu": int(len(A1)), "C0": bool(c0),
                    "quiral": bool(quir),
                    "mu_max": [int(x) for x in mm] if mm else None,
                    "n_maximales": int(nmax),
                    "A_mu_max": int(amm) if mm else None,
                    "senuelo_t2_discrepa": bool(sen["t2"]),
                    "senuelo_sin_filtro_discrepa": bool(sen["sf"]),
                    "senuelo_C_discrepa": bool(sen["C"]),
                    "n_eta_vivos": int(len(vivos)), "n_eta": int(len(ETAS)),
                    "TB_acierta": int(pred_ok)})

# ---------------------------------------------------------------- (H) impar, poblacion entera
POB = []
for (t, r, tope) in POBLACION:
    N = t + 2 * r
    BB = betas_estrictas(N, tope)
    n_nula = n_uni = n_pm1 = n_multi = 0
    peores = []
    for b in BB:
        P1 = phi_bialt(b, t, r)
        if P1 in (None, "NO-POL", "NO-RAC"):
            continue
        if not P1:
            n_nula += 1
            continue
        A1, rest1 = pelar(P1, "D", r)
        A1 = {k: QQ(v) for k, v in A1.items() if v != 0}
        if rest1 or not A1:
            n_nula += 1
            continue
        mm, amm, nmax = top_Omu(A1)
        if mm is None:
            n_multi += 1
            peores.append((b, nmax))
            continue
        n_uni += 1
        if abs(amm) == 1:
            n_pm1 += 1
        else:
            peores.append((b, str(mm) + " -> A=" + str(amm)))
    print("")
    print("  (H) IMPAR sobre la POBLACION ENTERA   t=%d, r=%d, beta con max <= %d, min 0 : %d formas"
          % (t, r, tope, len(BB)))
    print("      Phi = 0                          : %d" % n_nula)
    print("      Phi != 0 con mu^+ maximal UNICO  : %d" % n_uni)
    print("      Phi != 0 con varios maximales    : %d" % n_multi)
    print("      de los unicos, |A| = 1           : %d   (los que no: %s)"
          % (n_pm1, str(peores[:6]) if peores else "ninguno"))
    sys.stdout.flush()
    POB.append({"t": int(t), "r": int(r), "tope": int(tope), "n_formas": int(len(BB)),
                "n_nula": int(n_nula), "n_unico": int(n_uni), "n_multi": int(n_multi),
                "n_pm1": int(n_pm1),
                "excepciones": [[list(map(int, x[0])), str(x[1])] for x in peores[:20]]})

print("")
print("=" * 124)
print("  RESUMEN")
print("=" * 124)
if RES:
    print("  formas: %d | C0 ok: %d | resto B cero: %d | a^B enteros: %d | quiralidad ok: %d"
          % (len(RES), sum(1 for d in RES if d["C0"]), sum(1 for d in RES if d["resto_B_cero"]),
             sum(1 for d in RES if d["aB_enteros"]), sum(1 for d in RES if d["quiral"])))
    uni = [d for d in RES if d["mu_max"] is not None]
    print("  (H) impar: mu maximal unico en %d/%d | de esos, |A| = 1 en %d"
          % (len(uni), len(RES), sum(1 for d in uni if abs(d["A_mu_max"]) == 1)))
    print("  señuelos que DISCREPAN (tienen que discrepar todos): t+2 -> %d/%d | sin filtro -> %d/%d | C_m' -> %d/%d"
          % (sum(1 for d in RES if d["senuelo_t2_discrepa"]), len(RES),
             sum(1 for d in RES if d["senuelo_sin_filtro_discrepa"]), len(RES),
             sum(1 for d in RES if d["senuelo_C_discrepa"]), len(RES)))
print("")
print("  LECTURA, escrita ANTES de correr:")
print("   * si C0 pasa en todas -> el reparto B_{m'} x D_r es el correcto y su punto 1 esta medido.")
print("   * si ademas mu_max es unico con |A| = 1 -> el impar tiene su propia (H) y hay dicotomia.")
print("   * si C0 falla -> el 1 no va con la torsion, y hay que probar D_{m'} x B_r.")
json.dump({"C0": RES, "poblacion": POB}, open("odd_BD_DUMP.json", "w"), indent=1)
print("")
print("=" * 124)
print("DONE")
