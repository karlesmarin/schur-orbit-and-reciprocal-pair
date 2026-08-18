# -*- coding: utf-8 -*-
# LA CAPA DE SIGNOS DEL FILTRO: .que modulo ve?   16 de agosto de 2026.
#
# DE DONDE SALE.  Cruzando dos cosas de hoy que parecian sueltas:
#
#   (a)  galois_sign.sage encontro que en tipo B la formula cerrada vale SALVO UNA CONSTANTE
#            tau^B_t = eps_t . delta(a),        eps_t = +1, +1, -1, +1   en t = 3, 5, 7, 9
#        y esa constante no la enuncia el paper (prop:oddfilter solo da |tau^B| = 1).
#   (b)  el mismo guion partio  gamma_t(k) = sgn(sigma_k) . prod eps_j(k),  con  prod eps_j = (k/t)
#        por el lema de Gauss.  El factor  sgn(sigma_k)  quedo sin identificar.
#
# LOS DOS SON EL MISMO OBJETO: el det de una permutacion con signo del sistema medio {1..n}.
# Y hay un atajo que los vuelve ARITMETICA PURA, sin representaciones:
#
#        tau^B_t(0) = 1   (la representacion trivial vale 1 en cualquier punto)
#   =>   eps_t = delta( A_rho )   con   A_rho = (t-2, t-4, ..., 3, 1),
#
# o sea el det de la permutacion con signo que las clases plegadas de los IMPARES < t definen.
# Tres lineas de Python, para cualquier t.
#
# POR QUE IMPORTA.  Si eps_t tiene forma cerrada, prop:oddfilter sube de "|tau^B_t| = 1" a un SIGNO
# EXPLICITO, que es exactamente lo que le falta.  Y si el modulo que gobierna la capa de signos es
# t mod 8 mientras el soporte ve t mod 2, entonces la paridad del paper aparece una CUARTA vez y con
# otro modulo -- sec:closing dice ahora "una paridad vista tres veces".
#
# LO QUE SE MIDE
#   S1  eps_t para t impar hasta 61, contra candidatos de forma cerrada.
#   S2  sgn(sigma_k) como caracter de (Z/t)^x: .trivial, Jacobi, u otro?  Y su dependencia de t mod 8.
#   S3  gamma_t = sgn(sigma_k).(k/t), tabulado igual.
#   S4  el lema de Gauss como control: prod eps_j(k) == (k/t) para todo t, toda unidad.
#
# CONTROLES
#   C0  S4 es fatal: si el lema de Gauss falla, el plegado esta mal implementado.
#   C1  se imprime la tabla cruda por t, no solo el veredicto.
#   C2  el ajuste de forma cerrada se busca sobre t <= 31 y se PREDICE sobre 33..61, para que el
#       ajuste no pueda validarse a si mismo.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python sign_layer.py

import json
import sys
from math import gcd


def jacobi(a, n):
    """simbolo de Jacobi (a/n), n impar positivo."""
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def plegar(v, t):
    v %= t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def sgn_perm(perm):
    n = len(perm)
    s = 1
    visto = [False] * n
    for i in range(n):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = perm[j]
            L += 1
        if L % 2 == 0:
            s = -s
    return s


def delta_de(a, t):
    """det de la permutacion con signo que las clases plegadas de a definen; 0 si no es permutacion."""
    n = len(a)
    cl, ep = [], []
    for v in a:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, n + 1)):
        return 0, None, None
    s = sgn_perm([c - 1 for c in cl])
    pe = 1
    for e in ep:
        pe *= e
    return s * pe, s, pe


TS = list(range(3, 62, 2))
RES = []
print("=" * 112)
print("LA CAPA DE SIGNOS DEL FILTRO IMPAR")
print("=" * 112)
print("")
print("  t   n  t%8   eps_t   sgn(sigma_.)   gamma_.        (2/t)  (-1)^{n(n+1)/2}")
print("  " + "-" * 100)

gauss_ok = gauss_n = 0
for t in TS:
    n = (t - 1) // 2
    A_rho = [t - 2 * (j + 1) for j in range(n)]          # (t-2, t-4, ..., 1)
    eps_t, s_rho, pe_rho = delta_de(A_rho, t)

    unidades = [k for k in range(1, t) if gcd(k, t) == 1]
    sgn_s, prod_e, gam = {}, {}, {}
    for k in unidades:
        cl, ep = [], []
        for j in range(1, n + 1):
            c_, e_ = plegar(k * j, t)
            cl.append(c_)
            ep.append(e_)
        assert sorted(cl) == list(range(1, n + 1)), (t, k)
        s = sgn_perm([c - 1 for c in cl])
        pe = 1
        for e in ep:
            pe *= e
        sgn_s[k], prod_e[k] = s, pe
        gam[k] = s * pe
        gauss_n += 1
        if pe == jacobi(k, t):
            gauss_ok += 1

    def nombre(d):
        if all(v == 1 for v in d.values()):
            return "trivial"
        if all(d[k] == jacobi(k, t) for k in unidades):
            return "Jacobi"
        return "OTRO"

    ns, ng = nombre(sgn_s), nombre(gam)
    cand2 = jacobi(2, t)
    candn = 1 if (n * (n + 1) // 2) % 2 == 0 else -1
    print("  %2d %3d   %d   %+3d     %-9s      %-9s      %+3d     %+3d"
          % (t, n, t % 8, eps_t, ns, ng, cand2, candn))
    RES.append({"t": t, "n": n, "tmod8": t % 8, "eps_t": eps_t,
                "sgn_sigma": ns, "gamma": ng, "jacobi2": cand2, "cand_n": candn,
                "sgn_rho": s_rho, "prod_eps_rho": pe_rho})

print("")
print("  C0  LEMA DE GAUSS (fatal):  prod eps_j(k) == (k/t)  :  %d de %d" % (gauss_ok, gauss_n))

# ---------------------------------------------------------------- S1: forma cerrada de eps_t
print("")
print("=" * 112)
print("S1  .tiene eps_t forma cerrada?   Se AJUSTA sobre t <= 31 y se PREDICE sobre 33..61.")
print("=" * 112)

AJUSTE = [r for r in RES if r["t"] <= 31]
PRED = [r for r in RES if r["t"] > 31]

CANDS = {
    "(2/t)  = (-1)^{(t^2-1)/8}": lambda r: r["jacobi2"],
    "-(2/t)": lambda r: -r["jacobi2"],
    "(-1)^{n(n+1)/2}": lambda r: r["cand_n"],
    "(-1)^{floor(n/2)}": lambda r: 1 if (r["n"] // 2) % 2 == 0 else -1,
    "(-1)^{floor((t+1)/4)}": lambda r: 1 if ((r["t"] + 1) // 4) % 2 == 0 else -1,
    "por t mod 8: {1:+,3:+,5:+,7:-}": lambda r: -1 if r["tmod8"] == 7 else 1,
    "trivial (+1)": lambda r: 1,
}
buenos = []
for nombre_c, f in CANDS.items():
    ok_a = sum(1 for r in AJUSTE if f(r) == r["eps_t"])
    ok_p = sum(1 for r in PRED if f(r) == r["eps_t"])
    marca = "  <== SOBREVIVE" if (ok_a == len(AJUSTE) and ok_p == len(PRED)) else ""
    if ok_a == len(AJUSTE) and ok_p == len(PRED):
        buenos.append(nombre_c)
    print("   %-32s  ajuste %2d/%2d   prediccion %2d/%2d%s"
          % (nombre_c, ok_a, len(AJUSTE), ok_p, len(PRED), marca))

# ---------------------------------------------------------------- S2/S3: el caracter, por t mod 8
print("")
print("=" * 112)
print("S2/S3  el caracter sgn(sigma) y gamma, tabulados por t mod 8")
print("=" * 112)
tab = {}
for r in RES:
    tab.setdefault(r["tmod8"], []).append((r["t"], r["sgn_sigma"], r["gamma"]))
for k in sorted(tab):
    ss = set(x[1] for x in tab[k])
    gg = set(x[2] for x in tab[k])
    print("   t = %d (mod 8) :  t = %-28s  sgn(sigma) = %-20s  gamma = %s"
          % (k, str([x[0] for x in tab[k]])[:28], str(sorted(ss)), str(sorted(gg))))

# ---------------------------------------------------------------- S5: LA COMPUERTA DE LITERATURA
# Lo de arriba no es nuestro.  Se comprueba contra los enunciados publicados, que es la unica forma
# de que la compuerta valga algo:
#
#   Gauss (version de Jacobi):        prod eps_j(k) = (k/t)                    <- C0, ya arriba
#   Pan 2006, Teorema 1:              sgn(sigma_k)  = (k/t)^{(t+1)/2}
#     [H. Pan, "A remark on Zolotarev's theorem", arXiv:math/0601026;
#      enunciado citado en Z.-W. Sun, Finite Fields Appl. 59 (2019) 246-283, p. 2]
#   y por composicion:                gamma_t(k)    = (k/t)^{(t+3)/2}
#
# Y nuestro eps_t NO es un objeto nuevo: el vector congelado es
#     A_rho = (t-2, t-4, ..., 1)  =  (-2) . (1, 2, ..., n)   mod t,
# luego  eps_t = gamma_t(-2),  y Pan + Gauss lo dan cerrado.
print("")
print("=" * 112)
print("S5  COMPUERTA DE LITERATURA: .lo dicen ya Gauss y Pan?")
print("=" * 112)
pan_ok = pan_n = gam_ok = eps_g_ok = eps_c_ok = eps_n = 0
for t in TS:
    n = (t - 1) // 2
    unidades = [k for k in range(1, t) if gcd(k, t) == 1]
    for k in unidades:
        cl, ep = [], []
        for j in range(1, n + 1):
            c_, e_ = plegar(k * j, t)
            cl.append(c_)
            ep.append(e_)
        s = sgn_perm([c - 1 for c in cl])
        pe = 1
        for e in ep:
            pe *= e
        J = jacobi(k, t)
        pan_n += 1
        if s == (J ** ((t + 1) // 2) if J == -1 else 1) or s == pow(J, (t + 1) // 2):
            pan_ok += 1
        if s * pe == pow(J, (t + 3) // 2):
            gam_ok += 1
    # eps_t contra gamma_t(-2)  y contra la forma cerrada
    A_rho = [t - 2 * (j + 1) for j in range(n)]
    e_med = delta_de(A_rho, t)[0]
    kk = (-2) % t
    cl, ep = [], []
    for j in range(1, n + 1):
        c_, e_ = plegar(kk * j, t)
        cl.append(c_)
        ep.append(e_)
    g_m2 = sgn_perm([c - 1 for c in cl])
    for e in ep:
        g_m2 *= e
    eps_n += 1
    if e_med == g_m2:
        eps_g_ok += 1
    if e_med == pow(jacobi(kk, t), (t + 3) // 2):
        eps_c_ok += 1

print("   Pan 2006 Thm 1   sgn(sigma_k) == (k/t)^{(t+1)/2}   : %d de %d" % (pan_ok, pan_n))
print("   corolario        gamma_t(k)   == (k/t)^{(t+3)/2}   : %d de %d" % (gam_ok, pan_n))
print("   A_rho == (-2).(1..n)  =>  eps_t == gamma_t(-2)     : %d de %d" % (eps_g_ok, eps_n))
print("   y por tanto      eps_t == (-2/t)^{(t+3)/2}         : %d de %d" % (eps_c_ok, eps_n))
print("")
print("   VEREDICTO: si estas cuatro salen completas, la capa de signos es Gauss + Pan y no hay")
print("   nada nuestro en el mecanismo.  Lo nuestro es solo que el vector congelado de tipo B es")
print("   el sistema medio multiplicado por la unidad -2.")

json.dump({"por_t": RES, "gauss": [gauss_ok, gauss_n], "eps_cerrada": buenos,
           "pan": [pan_ok, pan_n], "gamma_cerrada": [gam_ok, pan_n],
           "eps_es_gamma_de_menos2": [eps_g_ok, eps_n], "eps_cerrada_jacobi": [eps_c_ok, eps_n]},
          open("sign_layer_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("DONE")
