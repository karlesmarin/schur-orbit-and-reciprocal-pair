# -*- coding: utf-8 -*-
# LA RECURSION DE PIERI SOBRE POBLACION, Y DONDE SE ROMPE.   16 de agosto de 2026.
#
# DE DONDE SALE.  pieri_recursion.sage midio, sobre nueve formas, que la identidad
#
#        sum_{rho = lambda + caja}  Phi(rho)  =  chi_natural . Phi(lambda)
#
# vale siempre (es un teorema: Phi es una evaluacion, luego un homomorfismo de anillos), y que en
# 7 de 9 UNA SOLA caja alcanza el peso superior del producto, heredando el coeficiente.  Eso seria
# un paso de induccion para la conjetura de la unidad.  Aqui se mide sobre poblacion y, sobre todo,
# se CARACTERIZA la excepcion: en las dos formas donde fallaba, ninguna caja alcanzaba el top, o
# sea varias se cancelaban arriba.
#
# LO QUE SE MIDE, sobre todas las lambda de una caja
#   Q1  el reparto de "cuantas cajas alcanzan el top del producto": 0, 1, 2, ...
#   Q2  cuando es 1: ¿hereda el coeficiente?  |A(rho)| == |A(lambda)|?
#   Q3  cuando es 0 o >1: ¿que tienen esas lambda en comun?  Se registra |lambda|, si Phi(lambda)
#       se anula, y cuantas de las rho se anulan -- que es el sospechoso natural.
#
# CONTROLES
#   C0  la identidad de Pieri se comprueba en TODAS, no solo en una muestra.  Es el control fatal.
#   C1  n impreso siempre y el reparto completo, no un promedio.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage pieri_population.sage

import json
import sys
from collections import Counter, defaultdict


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
    K = CyclotomicField(tt) if tt > 2 else QQ
    Lr = LaurentPolynomialRing(K, nvar, 'z')
    p = Lr(0)
    for g in Lr.gens():
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


def top(P, r):
    dom = [e for e in P if list(e) == sorted(e, reverse=True) and min(e) >= 0]
    return max(dom, key=lambda e: (sum(e), e)) if dom else None


def cajas(lam, N):
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
print("LA RECURSION DE PIERI SOBRE POBLACION")
print("=" * 120)

RES = []
for (t, r, tope) in [(4, 2, 9), (3, 2, 9), (6, 2, 8)]:
    N = t + 2 * r
    nat = natural(r, t)
    _c = {}

    def PH(lam):
        if lam not in _c:
            _c[lam] = phi(beta_de(lam, N), t, r)
        return _c[lam]

    LAMS = []
    for k in range(0, tope + 1):
        for e in Partitions(k, max_length=N):
            LAMS.append(tuple(e))
    pieri_ok = 0
    rep = Counter()
    hereda = 0
    con1 = 0
    excep = []
    n = 0
    for lam in LAMS:
        P = PH(lam)
        if not P:
            continue
        cs = cajas(lam, N)
        Ps = [PH(m) or {} for m in cs]
        if suma(Ps) != mult(nat, P):
            continue_ok = False
        else:
            continue_ok = True
        n += 1
        if continue_ok:
            pieri_ok += 1
        td = top(mult(nat, P), r)
        alc = [(m, Pm) for m, Pm in zip(cs, Ps) if Pm and top(Pm, r) == td]
        rep[len(alc)] += 1
        if len(alc) == 1:
            con1 += 1
            t0 = top(P, r)
            if t0 is not None and abs(alc[0][1][td]) == abs(P[t0]):
                hereda += 1
        elif len(excep) < 6:
            nulas = sum(1 for Pm in Ps if not Pm)
            excep.append((lam, len(alc), len(cs), nulas))
    print("")
    print("  t=%d r=%d  lambda hasta tamaño %d :  %d formas" % (t, r, tope, n))
    print("     C0  la identidad de Pieri : %d de %d" % (pieri_ok, n))
    print("     Q1  cajas que alcanzan el top del producto : %s" % dict(sorted(rep.items())))
    print("     Q2  con UNA sola caja: %d, y hereda |A| en %d de ellas" % (con1, hereda))
    if excep:
        print("     Q3  excepciones (lambda, #alcanzan, #cajas, #rho nulas): %s" % str(excep))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "n": int(n), "pieri_ok": int(pieri_ok),
                "reparto": {str(k): int(v) for k, v in rep.items()},
                "con_una": int(con1), "hereda": int(hereda),
                "excepciones": [[list(map(int, e[0])), int(e[1]), int(e[2]), int(e[3])] for e in excep]})

print("")
print("=" * 120)
print("  LECTURA, escrita ANTES de correr: si el reparto se concentra en 1 y hereda casi siempre,")
print("  hay paso de induccion y la conjetura de la unidad deja de estar sin ruta.")
json.dump(RES, open("pieri_population_DUMP.json", "w"), indent=1)
print("=" * 120)
print("DONE")
