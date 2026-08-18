# -*- coding: utf-8 -*-
# EL LOCUS DE ANULACION ES UN IDEAL, Y ESO LE PROHIBE COSAS.   16 de agosto de 2026.
#
# DE DONDE SALE.  Encadenando dos cosas que ya estan:
#
#   (a)  Phi es una EVALUACION, luego un homomorfismo de anillos.  De hecho es la composicion de
#        tres homomorfismos: Littlewood (restriccion), branching (restriccion), y -- por el teorema
#        de fusion -- la proyeccion de fusion minima.  Los tres son mapas de anillos.
#   (b)  luego  ker Phi  es un IDEAL, y la regla de Pieri  s_lambda . s_(1) = sum_{cajas} s_rho  da
#
#            Phi(lambda) = 0   ==>   sum_{rho = lambda + caja} Phi(rho) = 0.
#
# CONSECUENCIA QUE SE PUEDE FALSAR: una forma NULA no puede tener EXACTAMENTE UN vecino no nulo.
# Si lo tuviera, la suma seria ese unico Phi(rho) != 0, contradiccion.  Y mas fino: los Phi(rho) de
# los vecinos no nulos tienen que cancelarse entre ellos, o sea el locus de anulacion impone
# relaciones LINEALES sobre sus vecinos, no solo sobre si mismo.
#
# Eso es una restriccion sobre el locus que el criterio del Paper I no lleva escrita, y sale gratis.
#
# LO QUE SE MIDE
#   Z1  sobre TODAS las lambda nulas de la caja: el reparto de "cuantos vecinos lambda+caja son NO
#       nulos".  El valor 1 tiene que estar AUSENTE.
#   Z2  el control fatal: sum_{cajas} Phi(rho) = 0 en todas las nulas, coeficiente a coeficiente.
#   Z3  y al reves, sobre las NO nulas: el reparto de vecinos nulos, para tener con que comparar.
#   Z4  cuantas nulas tienen TODOS los vecinos nulos (el locus seria "cerrado hacia arriba" ahi).
#
# CONTROLES
#   C0  Z2 es un teorema: si falla, el guion esta mal, no el enunciado.
#   C1  n impreso siempre y repartos completos.
#   C2  SEÑUELO: el mismo conteo sobre lambda NO nulas -- ahi el 1 SI puede aparecer, y si no
#       apareciera nunca el enunciado no estaria diciendo nada sobre las nulas.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage zeros_propagate.sage

import json
import sys
from collections import Counter, defaultdict


def phi(beta, tt, nvar):
    Nn = tt + 2 * nvar
    K = CyclotomicField(tt) if tt > 2 else QQ
    zeta = K.gen() if tt > 2 else K(-1)
    Lr = LaurentPolynomialRing(K, nvar, 'z')
    xx = [Lr(K(zeta) ** k) for k in range(tt)] + [g ** e for g in Lr.gens() for e in (1, -1)]
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


print("=" * 116)
print("EL LOCUS DE ANULACION COMO IDEAL:  una forma nula no puede tener UN SOLO vecino no nulo")
print("=" * 116)

RES = []
for (t, r, tope) in [(4, 2, 12), (6, 2, 11), (3, 2, 12), (5, 2, 11)]:
    N = t + 2 * r
    _c = {}

    def PH(lam):
        if lam not in _c:
            _c[lam] = phi(beta_de(lam, N), t, r)
        return _c[lam]

    nulas, nonulas = [], []
    for k in range(0, tope + 1):
        for e in Partitions(k, max_length=N):
            lam = tuple(e)
            P = PH(lam)
            if P is None:
                continue
            (nulas if not P else nonulas).append(lam)

    rep_nula = Counter()
    rep_nonula = Counter()
    suma_ok = 0
    todos_nulos = 0
    for lam in nulas:
        cs = cajas(lam, N)
        Ps = [PH(m) for m in cs]
        nz = sum(1 for P in Ps if P)
        rep_nula[nz] += 1
        if nz == 0:
            todos_nulos += 1
        # control fatal: la suma tiene que ser cero
        S = defaultdict(lambda: QQ(0))
        for P in Ps:
            if P:
                for a, c in P.items():
                    S[a] += c
        if not {k: v for k, v in S.items() if v != 0}:
            suma_ok += 1
    for lam in nonulas[:400]:
        cs = cajas(lam, N)
        nz = sum(1 for m in cs if PH(m))
        rep_nonula[nz] += 1

    print("")
    print("  t=%d r=%d  lambda hasta tamaño %d :  %d nulas, %d no nulas"
          % (t, r, tope, len(nulas), len(nonulas)))
    print("     Z2  CONTROL FATAL, sum_{cajas} Phi(rho) = 0 : %d de %d" % (suma_ok, len(nulas)))
    print("     Z1  vecinos NO nulos de una forma nula : %s" % dict(sorted(rep_nula.items())))
    print("         -> ¿aparece el 1? %s" % ("SI, y eso refutaria el enunciado"
                                             if rep_nula.get(1) else "NO"))
    print("     Z4  nulas con TODOS los vecinos nulos : %d" % todos_nulos)
    print("     C2  vecinos NO nulos de una forma NO nula (muestra) : %s"
          % dict(sorted(rep_nonula.items())))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "tope": int(tope),
                "n_nulas": int(len(nulas)), "n_nonulas": int(len(nonulas)),
                "suma_cero": int(suma_ok),
                "reparto_nulas": {str(k): int(v) for k, v in rep_nula.items()},
                "reparto_nonulas": {str(k): int(v) for k, v in rep_nonula.items()},
                "todos_vecinos_nulos": int(todos_nulos)})

print("")
print("=" * 116)
print("  LECTURA, escrita ANTES de correr: si el 1 no aparece nunca entre las nulas y si aparece")
print("  entre las no nulas, el locus de anulacion tiene una restriccion combinatoria gratis que el")
print("  criterio del Paper I no lleva escrita.")
json.dump(RES, open("zeros_propagate_DUMP.json", "w"), indent=1)
print("=" * 116)
print("DONE")
