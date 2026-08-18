# -*- coding: utf-8 -*-
# UNA REGLA DE SELECCION QUE SALE DE CRUZAR DOS FORMULAS.   16 de agosto de 2026.
#
# DE DONDE SALE.  De cruzar formulas: estas dos nunca se habian compuesto:
#
#   (a)  el lado ORBITA:   eps_{lambda,nu} = s_{lambda/nu}(mu_t) != 0  exige que lambda/nu se pueda
#        teselar por t-ribbons, luego  t | |lambda| - |nu|.
#   (b)  el lado LIBRE:    b_{nu,mu} != 0  (la restriccion GL_2r -> Sp_2r de Littlewood) exige que
#        nu/mu se rellene con columnas pares, luego  |nu| - |mu|  es PAR.
#
# Componiendo por el cuadrado conmutativo,  A_mu != 0  =>  |mu| = |lambda| - t k - 2 j.  Y ahi entra
# la paridad de t, otra vez y por tercera via:
#
#     t PAR   ->  t k  es par  ->  |lambda| = |mu|  (mod 2)   REGLA DE SELECCION
#     t IMPAR ->  t k puede ser impar  ->  NO hay regla
#
# Si eso vale, Phi_{t,r} con t par vive entero en UNA clase de paridad de pesos, y eso es un criterio
# de anulacion barato: cualquier mu de la paridad equivocada tiene A_mu = 0 sin mirar nada mas.
#
# LO QUE SE MIDE
#   S1  para cada beta: el conjunto de |mu| mod 2 sobre los mu con A_mu != 0, y si es un solo valor.
#   S2  si lo es, ¿coincide con |lambda| mod 2?
#   S3  en el impar: ¿aparecen LAS DOS paridades?  Si no aparecieran, la regla seria mas fuerte de
#       lo que el argumento da y habria que entender por que.
#
# CONTROLES
#   C0  A_mu se calcula por bialternante + pelado, sin usar ni eps ni b: la regla se deduce de dos
#       formulas y se comprueba con una tercera ruta.
#   C1  n impreso siempre, y el reparto completo, no solo el veredicto.
#   C2  el caso par y el impar en el mismo barrido, para que la comparacion sea a igualdad de todo
#       lo demas.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage selection_rule.sage

import json
import sys
from collections import Counter


def phi_bialt(beta, tt, nvar):
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
        return "NO-POL"
    out = {}
    for e, c in zip(q.exponents(), q.coefficients()):
        k = tuple(int(v) for v in e) if hasattr(e, '__iter__') else (int(e),)
        if c not in QQ:
            return "NO-RAC"
        if QQ(c) != 0:
            out[k] = QQ(c)
    return out


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
    f = list(e)
    if typ in ("B", "C"):
        return f == sorted(f, reverse=True) and min(f) >= 0
    if len(f) < 2:
        return f[0] >= 0
    return all(f[i] >= f[i + 1] for i in range(len(f) - 2)) and f[-2] >= abs(f[-1])


def pelar(P, typ, rk, tope=8000):
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


def betas(N, tope):
    return [tuple(sorted(c, reverse=True)) + (0,) for c in Combinations(range(1, tope + 1), N - 1)]


print("=" * 118)
print("REGLA DE SELECCION POR PARIDAD DEL TAMAÑO   --   cruce del teorema de los ribbons con Littlewood")
print("=" * 118)
print("   t | r | tipo | formas | una sola paridad de |mu| | y coincide con |lambda| | las DOS paridades")
print("   " + "-" * 112)

RES = []
for (t, r, tope, tipo) in [(2, 2, 9, "C"), (4, 2, 11, "C"), (6, 2, 11, "C"),
                           (3, 2, 9, "D"), (5, 2, 10, "D")]:
    N = t + 2 * r
    delta = list(range(N - 1, -1, -1))
    n = una = coincide = dos = 0
    ejemplos = []
    for b in betas(N, tope):
        P = phi_bialt(b, t, r)
        if P in (None, "NO-POL", "NO-RAC") or not P:
            continue
        A, rest = pelar(P, tipo, r)
        A = {k: QQ(v) for k, v in A.items() if v != 0}
        if rest or not A:
            continue
        n += 1
        lam = [b[i] - delta[i] for i in range(N)]
        lsz = sum(v for v in lam if v > 0)
        par = set(sum(abs(v) for v in mu) % 2 for mu in A)
        if len(par) == 1:
            una += 1
            if list(par)[0] == lsz % 2:
                coincide += 1
            elif len(ejemplos) < 3:
                ejemplos.append((b, lsz % 2, list(par)[0]))
        else:
            dos += 1
    print("   %2d | %d | %-4s | %6d | %24d | %23d | %d"
          % (t, r, "par" if tipo == "C" else "impar", n, una, coincide, dos))
    if ejemplos:
        print("        una sola paridad pero NO la de |lambda|: %s" % str(ejemplos))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "tipo": tipo, "formas": int(n),
                "una_paridad": int(una), "coincide_con_lambda": int(coincide),
                "dos_paridades": int(dos)})

print("")
print("=" * 118)
print("  LECTURA, escrita ANTES de correr:")
print("   * si en el PAR sale una sola paridad SIEMPRE y coincide con |lambda|, la regla vale y es")
print("     un criterio de anulacion gratis: todo mu de la otra paridad tiene A_mu = 0.")
print("   * si en el IMPAR salen las dos, la regla es genuinamente de paridad de t -- tercera via")
print("     independiente por la que la paridad decide.")
print("   * si en el IMPAR tambien saliera una sola, el argumento da menos de lo que pasa y hay que")
print("     entender por que: seria un enunciado mas fuerte, no mas debil.")
json.dump(RES, open("selection_rule_DUMP.json", "w"), indent=1)
print("=" * 118)
print("DONE")
