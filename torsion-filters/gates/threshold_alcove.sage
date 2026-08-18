# -*- coding: utf-8 -*-
# .ES EL EMPATE DEL SENUELO UN UMBRAL DE NIVEL?   prob:threshold, 17 de agosto de 2026.
#
# EL PROBLEMA, TAL COMO EL PAPER LO PLANTEA.  El senuelo "el tau del t equivocado" empata en 3 de
# las 13 formas, y el paper conjetura que eso es el fenomeno de umbral conocido en fusion afin:
#
#     .coinciden dos filtros de orden distinto sobre Phi_{t,r} exactamente cuando el soporte de
#      B_{eta,mu} cae dentro del alcoba del menor de los dos?
#
# LO PRIMERO QUE HAY QUE MIRAR ES LA PREMISA, y es lo que este guion mide antes que nada.  El
# senuelo NO es "el mismo grupo a nivel mas alto": es tau al orden t+2 sobre C_{m+1}, con eta
# rellenado con un cero.  El alcoba de C_n al orden l es  {eta : eta_1 < l/2 - n}  --- en
# coordenadas de pesos fundamentales  sum m_i = eta_1 ---, luego
#
#     verdadero : l = t,   n = m   ->  eta_1 < t/2 - m         = 1
#     senuelo   : l = t+2, n = m+1 ->  eta_1 < (t+2)/2 - (m+1) = 1
#
# LOS DOS ESTAN EN EL NIVEL MINIMO.  Si eso sale confirmado, la pregunta no puede responderse como
# esta escrita: no hay un "menor de los dos" alcobas que separe nada, y el umbral de fusion afin no
# es el mecanismo.  Queda entonces la pregunta honesta, que es T2.
#
# LO QUE SE MIDE
#   T0  FATAL  reproducir el reparto de branch_filter: 3 empates y 10 discrepancias de 13.
#   T1  los dos alcobas, calculados de la desigualdad y no afirmados.
#   T2  .empatan porque los dos filtros COINCIDEN PUNTO A PUNTO sobre supp(B), o porque difieren
#       y aun asi las sumas con peso salen iguales?  Son dos mecanismos distintos y solo el
#       primero es un umbral.
#   T3  max eta_1 sobre supp(B), por forma.  Si los empates fueran un umbral, tendrian que ser
#       las formas de soporte mas bajo, y con un corte limpio.
#   T4  el conjunto donde los dos filtros coinciden punto a punto, barrido sobre una caja de eta:
#       .es una condicion de tipo alcoba (eta_1 acotado) o no lo es?
#
# CONTROLES
#   C0  T0 es fatal.  Si no reproduce 3/10, este guion no esta midiendo la misma poblacion.
#   C1  SENUELO DEL SENUELO: el mismo analisis con el orden t+4 sobre C_{m+2}.  Si el reparto de
#       empates fuese el mismo, lo medido seria "estar equivocado" y no el orden concreto.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage threshold_alcove.sage

import itertools, json

# --------------------------------------------------------------------------------------------
# Las seis funciones siguientes van COPIADAS VERBATIM de branch_filter.sage, no importadas: ese
# guion corre sus trece formas al cargarse y su salida esta archivada.  Copiar permite diff.
# --------------------------------------------------------------------------------------------
def phi_bialternante(beta, t, nvar):
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
    a, b = sp_char(eta, m), sp_char(mu, r)
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            k = e1 + e2
            out[k] = out.get(k, 0) + c1 * c2
    return out


def pelar_branching(P, m, r, tope=4000):
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


def clase(a, t, m):
    c = a % t
    if c == 0 or 2 * c == t:
        return None
    return (c, +1) if c <= m else (t - c, -1)


def tau(eta, t, m):
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


r = 2
CASOS = {
    4: [(18, 17, 11, 8, 7, 6, 1, 0), (10, 9, 7, 4, 3, 2, 1, 0), (14, 13, 11, 4, 3, 2, 1, 0),
        (12, 11, 10, 5, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
        (13, 9, 8, 7, 5, 4, 2, 0), (19, 17, 11, 8, 7, 6, 1, -1), (21, 17, 11, 8, 7, 6, 1, -3)],
    6: [(13, 11, 9, 7, 5, 4, 3, 2, 1, 0), (15, 13, 11, 8, 6, 5, 3, 2, 1, 0),
        (17, 13, 11, 9, 7, 5, 3, 2, 1, 0), (12, 11, 10, 9, 8, 5, 3, 2, 1, 0)],
}


def filtro_desplazado(eta, t, m, salto):
    """tau al orden t+salto sobre C_{m+salto/2}, con eta rellenado de ceros.  El senuelo."""
    tp = t + salto
    mp = (tp - 2) // 2
    return tau(tuple(list(eta) + [0] * (mp - len(eta)))[:mp], tp, mp)[0]


def alcoba_cota(l, n):
    """el alcoba de C_n al orden l es {eta_1 < l/2 - n}; devuelve la cota."""
    return Rational(l) / 2 - n


print("=" * 108)
print("EL EMPATE DEL SENUELO: .UMBRAL DE NIVEL O CANCELACION?   prob:threshold")
print("=" * 108)
print("")

print("  T1  los dos alcobas, de la desigualdad:")
for t in (4, 6):
    m = (t - 2) // 2
    print("     t=%d  m=%d : filtro verdadero  C_%d al orden %d  ->  eta_1 < %s"
          % (t, m, m, t, alcoba_cota(t, m)))
    for salto in (2, 4):
        tp, mp = t + salto, (t + salto - 2) // 2
        print("              senuelo  t+%d      C_%d al orden %d  ->  eta_1 < %s"
              % (salto, mp, tp, alcoba_cota(tp, mp)))
print("")

filas, agregado = [], {"empate": 0, "discrepa": 0, "pw_igual_en_empate": 0,
                       "pw_igual_en_discrepa": 0, "n": 0}
for t in sorted(CASOS):
    m = (t - 2) // 2
    R = r + m
    for beta in CASOS[t]:
        P = phi_bialternante(list(beta), 2, R)          # el objeto con todos los pares libres
        if P is None or P == "NO-POLINOMIO":
            continue
        B, resto = pelar_branching(dict(P), m, r)
        if resto:
            print("  !! peeling con resto en beta=%s -- fila descartada" % (beta,))
            continue
        B = {k: v for k, v in B.items() if v != 0}
        sop = sorted({eta for (eta, mu) in B})
        Phi = phi_bialternante(list(beta), t, r)
        A_ok = ruta_de_el(B, t, m, r, lambda e: tau(e, t, m)[0])
        mon_ok = a_monomios(A_ok, r)
        veredicto = {}
        for salto in (2, 4):
            A_x = ruta_de_el(B, t, m, r, lambda e: filtro_desplazado(e, t, m, salto))
            mon_x = a_monomios(A_x, r)
            empata = (mon_ok == mon_x)
            # T2: .coinciden PUNTO A PUNTO sobre el soporte?
            dif = [e for e in sop if tau(e, t, m)[0] != filtro_desplazado(e, t, m, salto)]
            veredicto[salto] = (empata, len(dif), len(sop))
        emp2, dif2, nsop = veredicto[2]
        maxeta1 = max([e[0] for e in sop]) if sop else 0
        agregado["n"] += 1
        agregado["empate" if emp2 else "discrepa"] += 1
        if emp2 and dif2 == 0:
            agregado["pw_igual_en_empate"] += 1
        if (not emp2) and dif2 == 0:
            agregado["pw_igual_en_discrepa"] += 1
        filas.append({"t": int(t), "beta": [int(b) for b in beta], "n_B": int(len(B)),
                      "n_sop": int(nsop), "max_eta1": int(maxeta1), "empata_t2": bool(emp2),
                      "difieren_pw_t2": int(dif2), "empata_t4": bool(veredicto[4][0]),
                      "difieren_pw_t4": int(veredicto[4][1])})

print("  T0/T2/T3  por forma:")
print("     t | beta                              | #B  |sop| max eta_1 | senuelo t+2   | dif p.a p. | senuelo t+4")
print("     " + "-" * 118)
for f in filas:
    print("     %d | %-33s | %3d | %2d |    %3d    | %-13s | %4d de %2d | %s"
          % (f["t"], str(tuple(f["beta"])), f["n_B"], f["n_sop"], f["max_eta1"],
             "EMPATA" if f["empata_t2"] else "discrepa (ok)", f["difieren_pw_t2"], f["n_sop"],
             "EMPATA" if f["empata_t4"] else "discrepa (ok)"))
print("")
print("  T0  FATAL  empates con el senuelo t+2 : %d de %d   (branch_filter dice 3 de 13)"
      % (agregado["empate"], agregado["n"]))
print("  T2  de los que EMPATAN, cuantos coinciden punto a punto sobre supp(B) : %d de %d"
      % (agregado["pw_igual_en_empate"], agregado["empate"]))
print("  T2  de los que DISCREPAN, cuantos coinciden punto a punto              : %d de %d"
      % (agregado["pw_igual_en_discrepa"], agregado["discrepa"]))
print("  C1  SENUELO DEL SENUELO, empates con t+4 : %d de %d"
      % (sum(1 for f in filas if f["empata_t4"]), agregado["n"]))
print("")

# ---------------------------------------------------------------- T4: donde coinciden los dos ----
REGION = {}
print("  T4  el conjunto {eta : tau_t(eta) == tau_{t+2}(eta,0)}, barrido sobre eta_1 <= 12:")
for t in (4, 6):
    m = (t - 2) // 2
    iguales, distintos, por_eta1 = [], [], {}
    for eta in itertools.product(range(13), repeat=m):
        if any(eta[i] < eta[i + 1] for i in range(m - 1)):
            continue
        a, b = tau(eta, t, m)[0], filtro_desplazado(eta, t, m, 2)
        (iguales if a == b else distintos).append(eta)
        por_eta1.setdefault(eta[0], [0, 0])[0 if a == b else 1] += 1
    cotas = sorted(por_eta1)
    solo_ig = [e for e in cotas if por_eta1[e][1] == 0]
    print("     t=%d : coinciden en %d de %d pesos.  Por eta_1, (coinciden, difieren):"
          % (t, len(iguales), len(iguales) + len(distintos)))
    print("        %s" % {e: tuple(por_eta1[e]) for e in cotas[:13]})
    print("        eta_1 con coincidencia TOTAL: %s  ->  %s"
          % (solo_ig, "es un corte inicial (aspecto de alcoba)"
             if solo_ig == list(range(len(solo_ig))) else "NO es un corte inicial"))
    REGION[int(t)] = {"por_eta1": {int(e): [int(por_eta1[e][0]), int(por_eta1[e][1])] for e in cotas},
                      "total_iguales": int(len(iguales)),
                      "total": int(len(iguales) + len(distintos)),
                      "eta1_totalmente_iguales": [int(e) for e in solo_ig],
                      "corte_inicial": bool(solo_ig == list(range(len(solo_ig))))}
print("")

json.dump({"filas": filas, "agregado": {k: int(v) for k, v in agregado.items()},
           "region_pw": REGION}, open("threshold_alcove_DUMP.json", "w"), indent=1)
print("  LECTURA, escrita ANTES de correr:")
print("   * si T1 da la MISMA cota para los dos, la pregunta de prob:threshold no puede")
print("     responderse como esta escrita: los dos filtros estan en el nivel minimo.")
print("   * si T2 dice que los que empatan coinciden PUNTO A PUNTO, el empate es del filtro y")
print("     hay que buscar la region de coincidencia, que es T4.")
print("   * si T2 dice que difieren punto a punto y aun asi empatan, el empate es una")
print("     CANCELACION y no hay umbral ninguno: prob:threshold hay que reescribirlo.")
print("")
print("=" * 108)
print("DONE")
