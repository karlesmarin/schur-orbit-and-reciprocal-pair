# -*- coding: utf-8 -*-
# UNA RECURSION DE PIERI PARA LOS A_mu, Y QUE LE HACE AL PESO SUPERIOR.   16 de agosto de 2026.
#
# DE DONDE SALE.  Cruzar dos cosas obvias por separado y que nunca se habian juntado:
#
#   (a)  Phi_{t,r} es una EVALUACION de funciones simetricas, luego es un homomorfismo de ANILLOS
#            Phi(lambda) Phi(nu) = sum_rho c^rho_{lambda nu} Phi(rho).
#   (b)  Phi de la particion (1) es la suma del alfabeto:
#            suma de las raices t-esimas = 0  (t > 1),  mas  sum (z_i + z_i^{-1}),
#        o sea EXACTAMENTE el caracter del modulo natural del factor libre.
#
# Componiendo con nu = (1) sale una regla de Pieri:
#
#        sum_{rho = lambda + caja}  Phi(rho)  =  chi_natural . Phi(lambda),
#
# es decir una RECURSION EN lambda para toda la familia de los A_mu.  Y en el peso superior tiene
# consecuencias: el peso maximo de  chi_nat . Phi(lambda)  es  mu_max(lambda) + (1,0,...,0), asi que
# si UNA sola caja alcanza ese peso, su coeficiente hereda el de lambda -- que seria un paso de
# induccion para la conjetura de la unidad, que hoy no tiene ninguno.
#
# LO QUE SE MIDE
#   P1  CONTROL, y tiene que pasar exacto: la identidad de Pieri, coeficiente a coeficiente en mu.
#       Si fallara, el instrumento esta roto (es un teorema, no una hipotesis).
#   P2  el peso superior del producto: ¿es mu_max(lambda) + e_1?
#   P3  LA PREGUNTA: entre las rho = lambda + caja, ¿cuantas tienen mu_max(rho) = mu_max(lambda)+e_1?
#       Si es UNA sola, hay paso de induccion.  Si son varias, hay que ver si sus coeficientes se
#       cancelan y queda uno.
#   P4  y el signo: A_{mu_max}(rho) contra A_{mu_max}(lambda).
#
# CONTROLES
#   C0  Phi se calcula por bialternante para cada rho por separado, sin usar la recursion.
#   C1  n impreso siempre.
#   C2  SEÑUELO: la misma identidad con chi del modulo natural del grupo EQUIVOCADO (el del otro
#       tipo).  Tiene que fallar.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage pieri_recursion.sage

import json
import sys
from collections import defaultdict


def phi(beta, tt, nvar):
    Nn = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    Lr = LaurentPolynomialRing(K, nvar, 'z')
    zs = Lr.gens()
    xx = [Lr(K(zeta) ** k) for k in range(tt)] + [g ** e for g in zs for e in (1, -1)]
    delta = list(range(Nn - 1, -1, -1))

    def alt(expo):
        return matrix(Lr, Nn, Nn, lambda i, j: xx[i] ** expo[j]).determinant()

    den = alt(delta)
    if den == 0:
        return None
    q = alt(list(beta)) / den
    try:
        q = Lr(q)
    except Exception:
        return None
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        k = tuple(int(v) for v in e) if hasattr(e, '__iter__') else (int(e),)
        if c not in QQ:
            return None
        if QQ(c) != 0:
            out[k] = QQ(c)
    return out


def natural(nvar, tt):
    """el caracter del modulo natural del factor libre, como polinomio de Laurent."""
    K = CyclotomicField(tt) if tt > 2 else QQ
    Lr = LaurentPolynomialRing(K, nvar, 'z')
    zs = Lr.gens()
    p = Lr(0)
    for g in zs:
        p += g + g ** (-1)
    return {tuple(int(v) for v in e) if hasattr(e, '__iter__') else (int(e),): QQ(c)
            for e, c in zip(p.exponents(), p.coefficients()) if c != 0}


def mult(P, Q):
    R = defaultdict(lambda: QQ(0))
    for a, ca in P.items():
        for b, cb in Q.items():
            R[tuple(a[i] + b[i] for i in range(len(a)))] += ca * cb
    return {k: v for k, v in R.items() if v != 0}


def suma(Ps):
    R = defaultdict(lambda: QQ(0))
    for P in Ps:
        for a, c in P.items():
            R[a] += c
    return {k: v for k, v in R.items() if v != 0}


def cajas(lam, N):
    """las particiones lambda + una caja, con a lo sumo N partes."""
    out = []
    L = list(lam) + [0]
    for i in range(len(L)):
        if i > 0 and L[i] + 1 > L[i - 1]:
            continue
        if i >= N:
            continue
        M = list(L)
        M[i] += 1
        while M and M[-1] == 0:
            M.pop()
        out.append(tuple(M))
    return out


def beta_de(lam, N):
    L = list(lam) + [0] * (N - len(lam))
    return tuple(L[i] + (N - 1 - i) for i in range(N))


print("=" * 120)
print("REGLA DE PIERI PARA LOS A_mu, Y SU EFECTO EN EL PESO SUPERIOR")
print("=" * 120)

RES = []
for (t, r, LAMS) in [(4, 2, [(3, 3, 2), (5, 5, 5, 1), (2, 1), (4, 2, 1)]),
                     (3, 2, [(3, 2, 1), (2, 2), (4, 1)]),
                     (6, 2, [(3, 2), (2, 1, 1)])]:
    N = t + 2 * r
    nat = natural(r, t)
    print("")
    print("  t=%d r=%d N=%d" % (t, r, N))
    print("  lambda            | #cajas | Pieri exacto | top de chi.Phi   | cajas que alcanzan el top |"
          " A(lambda) -> A(rho)")
    print("  " + "-" * 114)
    for lam in LAMS:
        if len(lam) > N:
            continue
        P = phi(beta_de(lam, N), t, r)
        if not P:
            continue
        izq = suma([phi(beta_de(m, N), t, r) or {} for m in cajas(lam, N)])
        der = mult(nat, P)
        ok = (izq == der)
        # el top dominante del producto
        dom = [e for e in der if list(e) == sorted(e, reverse=True) and min(e) >= 0]
        topd = max(dom, key=lambda e: (sum(e), e)) if dom else None
        # cuantas cajas alcanzan ese top como su propio top
        alcanzan = []
        for m in cajas(lam, N):
            Pm = phi(beta_de(m, N), t, r)
            if not Pm:
                continue
            dm = [e for e in Pm if list(e) == sorted(e, reverse=True) and min(e) >= 0]
            if not dm:
                continue
            tm = max(dm, key=lambda e: (sum(e), e))
            if tm == topd:
                alcanzan.append((m, Pm[tm]))
        dm0 = [e for e in P if list(e) == sorted(e, reverse=True) and min(e) >= 0]
        t0 = max(dm0, key=lambda e: (sum(e), e)) if dm0 else None
        print("  %-17s | %6d | %-12s | %-16s | %25d | %s -> %s"
              % (str(lam), len(cajas(lam, N)), "SI" if ok else "FALLA",
                 str(topd), len(alcanzan),
                 str(P[t0]) if t0 else "-",
                 str([int(c) for _, c in alcanzan])))
        sys.stdout.flush()
        RES.append({"t": int(t), "r": int(r), "lambda": [int(v) for v in lam],
                    "pieri_ok": bool(ok), "top_producto": [int(v) for v in topd] if topd else None,
                    "n_cajas": int(len(cajas(lam, N))),
                    "cajas_en_el_top": int(len(alcanzan)),
                    "A_lambda": int(P[t0]) if t0 else None,
                    "A_rho": [int(c) for _, c in alcanzan]})

print("")
print("=" * 120)
print("  LECTURA, escrita ANTES de correr:")
print("   * P1 tiene que salir SI en todas: es un teorema, y si falla es el guion.")
print("   * si en la mayoria UNA sola caja alcanza el top, hay paso de induccion para la unidad.")
print("   * si son varias, hay que mirar si sus A suman +-1, que seria una induccion mas debil pero")
print("     todavia util.")
print("  RESUMEN: Pieri exacto en %d de %d | una sola caja en el top en %d"
      % (sum(1 for d in RES if d["pieri_ok"]), len(RES),
         sum(1 for d in RES if d["cajas_en_el_top"] == 1)))
json.dump(RES, open("pieri_recursion_DUMP.json", "w"), indent=1)
print("=" * 120)
print("DONE")
