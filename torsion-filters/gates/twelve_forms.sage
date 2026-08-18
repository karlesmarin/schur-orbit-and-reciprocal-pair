# -*- coding: utf-8 -*-
# LAS DOCE.  El residuo del problema abierto, y es pequeño.   15 de agosto de 2026.
#
# DE DONDE SALE, Y POR QUE SON DOCE.  A W=13, t=4, r=2 hay 1716 formas.  De ellas:
#
#     275  NO ocupadas   -> Phi = 0 SIEMPRE, por el criterio clasico de t-core (type-A de julio:
#                           alguna clase residual mod t vacia => bloque no cuadrado => 0)
#    1441  ocupadas      -> de estas, 21 se anulan igual
#                             9  HEREDADAS  (Phi_{2,R} = 0 ya antes de especializar)
#                            12  POR CANCELACION  <-- EL RESIDUO.  Estas doce.
#
# O sea: de los 296 ceros, 275 (93 %) los ve el criterio clasico y 9 mas los ve t=2.  El problema
# abierto son DOCE formas de 1716.  Ese es el objeto de este guion.
#
# QUE SE BUSCA.  Un discriminante: que tienen las 12 que no tengan las 1420 no nulas ocupadas.
# Con n=12 y muchos rasgos, ajustar ruido es trivial ([[a-coincidence-at-one-parameter]]), asi que
# TODA regla candidata se evalua ACTO SEGUIDO sobre la poblacion no nula completa.  Una regla que no
# se testee contra las 1420 no se reporta.
#
# RASGOS QUE SE VUELCAN, por forma (todos calculables desde beta, sin conocer la respuesta):
# AVISO (15-ago, tras quotient_split.sage): la columna "tamaños cuota" y "cuota_pesos" de este
# guion se calculan desde un beta-set RECOMPUTADO desde lambda con n variable, y por eso NO son
# canonicas -- ver §5 de RESIDUO_12_FORMAS.md.  La columna t_core SI es valida (el core no depende
# de n, verificado).  El beta-set canonico es beta, de tamaño N = t+2r.
#
#   lambda = beta - delta, su longitud, si esta en el RANGO ESTABLE l(lambda) <= r+1 (§14 de julio)
#   si lambda es rectangulo (k^n) con k impar   <- el vanisher del rango estable de julio
#   el t-core y el t-quotient de lambda         <- el lenguaje del motor R
#   la paridad |lambda| mod 2 y |beta| mod 2
#   el numero de eta regulares en supp(B), y cuantos mu tienen A_mu != 0 antes de cancelar
#
# CONTROLES
#   C0  FATAL.  Las 12 tienen que salir con Phi = 0 recalculado aqui por bialternante independiente.
#   C1  FATAL.  Las 1420 de control tienen que salir Phi != 0.  Si alguna diera 0, la lista esta mal.
#   C2  CADA regla candidata se aplica a la poblacion NO NULA y se reporta su tasa de falsos
#       positivos.  Una regla con 0 falsos positivos sobre 1420 es un candidato; con >0 se descarta
#       y se dice cuantos.
#   C3  no vacuidad: n impreso siempre, y las 12 filas VISIBLES una a una -- son doce, caben.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage twelve_forms.sage

import itertools, json, sys
from collections import defaultdict

load("pob_helper.py")


def phi_de(beta, tt, nvar):
    N = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    L = LaurentPolynomialRing(K, nvar, 'z')
    zs = L.gens()
    x = [L(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(N - 1, -1, -1))
    def alt(ex):
        return matrix(L, N, N, lambda i, j: x[i] ** ex[j]).determinant()
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


# ------------------------------------------------------------------ rasgos ----------------------
def lam_de(beta, N):
    d = list(range(N - 1, -1, -1))
    return tuple(beta[i] - d[i] for i in range(N))

def core_quotient(lam, tt):
    """t-core y t-quotient por el beta-set: clases residuales de beta = lam + delta."""
    n = len(lam)
    beta = [lam[i] + (n - 1 - i) for i in range(n)]
    clases = defaultdict(list)
    for b in beta:
        clases[b % tt].append(b)
    cuota = {}
    for q in range(tt):
        vs = sorted(clases[q], reverse=True)
        k = len(vs)
        cuota[q] = tuple((vs[i] - q) // tt - (k - 1 - i) for i in range(k))
    # el core: reempaquetar cada clase densa desde abajo
    nuevo = []
    for q in range(tt):
        k = len(clases[q])
        nuevo += [q + tt * i for i in range(k)]
    nuevo.sort(reverse=True)
    core = tuple(nuevo[i] - (len(nuevo) - 1 - i) for i in range(len(nuevo)))
    return tuple(x for x in core if x), {q: tuple(x for x in cuota[q] if x) for q in cuota}

def rasgos(beta, tt, r):
    N = tt + 2 * r
    lam = lam_de(beta, N)
    partes = tuple(x for x in lam if x != 0)
    L = len(partes)
    n_est = r + 1
    core, cuota = core_quotient(lam, tt)
    return {
        "beta": [int(x) for x in beta],
        "lambda": [int(x) for x in partes],
        "l_lambda": int(L),
        "estable": bool(L <= n_est),
        "rect_impar": bool(L == n_est and len(set(partes)) == 1 and partes and partes[0] % 2 == 1),
        "abs_lambda": int(sum(partes)),
        "par_lambda": int(sum(partes) % 2),
        "par_beta": int(sum(beta) % 2),
        "t_core": [int(x) for x in core],
        "t_core_vacio": bool(len(core) == 0),
        "t_core_len": int(len(core)),
        "cuota_tam": [int(len(cuota[q])) for q in sorted(cuota)],
        "cuota_pesos": [int(sum(cuota[q])) for q in sorted(cuota)],
    }


# ================================================================== corrida =====================
CONF = [(4, 2, 13), (6, 2, 13)]
SALIDA = {}
for (t, r, W) in CONF:
    N = t + 2 * r
    print("=" * 126)
    print("LAS QUE SOBREVIVEN AL CRITERIO CLASICO   --   t=%d  r=%d  W=%d  N=%d" % (t, r, W, N))
    print("=" * 126)
    print("")
    todas = [tuple(b) for b in betas_py(t, r, W)]
    ocup = [b for b in todas if occupied_py(b, t)]
    nulas = [b for b in ocup if phi_zero_py(b, t, r)]
    nonulas = [b for b in ocup if not phi_zero_py(b, t, r)]
    print("   poblacion: %d formas | %d ocupadas | %d nulas ocupadas | %d no nulas"
          % (len(todas), len(ocup), len(nulas), len(nonulas)))
    print("   (las %d NO ocupadas son cero por el criterio clasico de t-core y no entran aqui)"
          % (len(todas) - len(ocup)))
    print("")
    sys.stdout.flush()

    # C0 / C1 -- recalcular por bialternante independiente
    malo0 = sum(1 for b in nulas if phi_de(b, t, r))
    ctrl = nonulas[::max(1, len(nonulas) // 40)][:40]
    malo1 = sum(1 for b in ctrl if not phi_de(b, t, r))
    print("   C0  las %d nulas recalculadas por bialternante dan 0 : %s"
          % (len(nulas), "PASA" if malo0 == 0 else "*** FALLA en %d ***" % malo0))
    print("   C1  %d no nulas de control dan != 0                  : %s"
          % (len(ctrl), "PASA" if malo1 == 0 else "*** FALLA en %d ***" % malo1))
    print("")

    R_nulas = [rasgos(b, t, r) for b in nulas]
    R_no = [rasgos(b, t, r) for b in nonulas]
    SALIDA["t%d_r%d" % (t, r)] = {"nulas": R_nulas, "n_no_nulas": len(R_no),
                                  "no_nulas_muestra": R_no[::max(1, len(R_no) // 60)][:60]}

    print("   LAS %d, UNA A UNA:" % len(nulas))
    print("   %-30s | %-20s | l | est | rect | |l| | par | t-core        | tamaños cuota" % ("beta", "lambda"))
    print("   " + "-" * 122)
    for x in R_nulas:
        print("   %-30s | %-20s | %d | %-3s | %-4s | %3d | %d   | %-13s | %s"
              % (str(tuple(x["beta"])), str(tuple(x["lambda"])), x["l_lambda"],
                 "si" if x["estable"] else "no", "si" if x["rect_impar"] else "no",
                 x["abs_lambda"], x["par_lambda"], str(tuple(x["t_core"])), x["cuota_tam"]))
    print("")
    sys.stdout.flush()

    # ---------------------------------------------------------- C2: reglas candidatas -----------
    print("   C2  CADA regla que separa las nulas se evalua sobre LAS %d NO NULAS." % len(R_no))
    print("       Una regla solo es candidata si su tasa de falsos positivos es 0 sobre esas %d." % len(R_no))
    print("")
    REGLAS = [
        ("t_core vacio",            lambda x: x["t_core_vacio"]),
        ("l(lambda) > r+1 (inestable)", lambda x: not x["estable"]),
        ("rectangulo impar (§14)",   lambda x: x["rect_impar"]),
        ("|lambda| par",             lambda x: x["par_lambda"] == 0),
        ("|lambda| impar",           lambda x: x["par_lambda"] == 1),
        ("Sigma beta par",           lambda x: x["par_beta"] == 0),
        ("cuota equilibrada",        lambda x: len(set(x["cuota_tam"])) == 1),
        ("|t_core| <= 1",            lambda x: x["t_core_len"] <= 1),
    ]
    print("   %-32s | cubre nulas | falsos + sobre no nulas | veredicto" % "regla")
    print("   " + "-" * 110)
    for nom, f in REGLAS:
        cn = sum(1 for x in R_nulas if f(x))
        fp = sum(1 for x in R_no if f(x))
        if cn == len(R_nulas) and fp == 0:
            v = "*** CANDIDATA: cubre todas y 0 falsos ***"
        elif cn == len(R_nulas):
            v = "cubre todas pero %d falsos -> NO separa" % fp
        elif fp == 0:
            v = "0 falsos pero solo cubre %d/%d -> parcial" % (cn, len(R_nulas))
        else:
            v = "no"
        print("   %-32s | %5d/%-5d | %22d | %s" % (nom, cn, len(R_nulas), fp, v))
    print("")
    sys.stdout.flush()

json.dump(SALIDA, open("twelve_forms_DUMP.json", "w"), indent=1)
print("=" * 126)
print("  volcado completo en twelve_forms_DUMP.json  (las nulas enteras + muestra de las no nulas)")
print("=" * 126)
print("DONE")
