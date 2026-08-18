# -*- coding: utf-8 -*-
# ¿VALE |A_mu| <= 1 PARA TODO mu, O SOLO PARA EL MAXIMO?   15 de agosto de 2026.
#
# DE DONDE SALE, Y POR QUE AHORA.  Yacobi (arXiv:0907.3247), verificado hoy en su abstract:
#
#   "each multiplicity space that arises in the restriction of an irreducible representation of
#    Sp_{2n} to Sp_{2n-2} is canonically an irreducible module for the n-fold product of SL_2 ...
#    This induces a canonical decomposition of the multiplicity spaces into one dimensional spaces"
#
# Nuestro A_mu es el caracter del espacio de multiplicidad M_mu evaluado en el elemento de torsion.
# Si M_mu es un producto tensorial de SL_2-modulos, ese caracter es un PRODUCTO de valores que, por
# nuestro propio lema con m=1, son 0 o +-1.  Luego la ruta de Yacobi PREDICE
#
#       |A_mu| <= 1   PARA TODO mu,   no solo para el maximo.
#
# Es una prediccion falsable y barata, y es la primera vez que la estructura del branching y la del
# filtro se tocan en vez de estar simplemente compuestas.  Si aparece |A_mu| >= 2, la ruta no cierra
# y hay que saberlo ANTES de escribir una linea sobre ella.
#
# CONTROLES
#   C0  FATAL.  sum_mu A_mu sp_mu tiene que reproducir Phi_{t,r} monomio a monomio.  Sin eso el
#       histograma de |A_mu| no es de nada.
#   C1  se imprime el histograma ENTERO de |A_mu|, no un maximo: un maximo esconde la forma.
#   C2  t=4 (m=1, el caso de Yacobi) Y t=6 (m=2, donde su teorema NO aplica directamente).  Si el
#       acotamiento fuera un accidente de la poblacion y no de m=1, saldria igual en los dos, y eso
#       tambien hay que verlo.
#   C3  no vacuidad: n impreso siempre, y el numero de mu por forma.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage all_A_mu.sage

import itertools, json, sys
from collections import defaultdict


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


def pelar_branching(P, m, r, tope=9000):
    P = {e: QQ(c) for e, c in P.items() if c != 0}
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


def tau(eta, tt, mm):
    a = [eta[j] + (mm - (j + 1) + 1) for j in range(mm)]
    cl, sg = [], 1
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


CASOS = {
    4: [(18, 17, 11, 8, 7, 6, 1, 0), (10, 9, 7, 4, 3, 2, 1, 0), (14, 13, 11, 4, 3, 2, 1, 0),
        (12, 11, 10, 5, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
        (13, 9, 8, 7, 5, 4, 2, 0)],
    6: [(13, 11, 9, 7, 5, 4, 3, 2, 1, 0), (12, 11, 10, 9, 8, 4, 3, 2, 1, 0),
        (14, 13, 10, 6, 5, 4, 3, 2, 1, 0), (12, 10, 9, 8, 5, 4, 3, 2, 1, 0)],
}
r = 2
print("=" * 120)
print("EL HISTOGRAMA COMPLETO DE |A_mu|,  no solo el del peso maximo")
print("=" * 120)
RES = {}
for t in (4, 6):
    m = (t - 2) // 2
    R = r + m
    hist = defaultdict(int)
    malo0 = n = 0
    print("")
    print("  t=%d  m=%d  R=%d   (Yacobi aplica directamente solo a m=1)" % (t, m, R))
    print("  %-30s | #mu | |A_mu| observados                  | C0" % "beta")
    print("  " + "-" * 108)
    for b in CASOS[t]:
        Psi = phi_bialternante(b, 2, R)
        Phi = phi_bialternante(b, t, r)
        if Psi in (None, "NO-POL") or Phi in (None, "NO-POL"):
            continue
        n += 1
        B, resto = pelar_branching({k: QQ(v) for k, v in Psi.items()}, m, r)
        A = defaultdict(lambda: 0)
        for (eta, mu), bb in B.items():
            v = tau(eta, t, m)
            if v:
                A[mu] += bb * v
        A = {mu: a for mu, a in A.items() if a != 0}
        # C0: reconstruir
        rec = {}
        for mu, a in A.items():
            for k, v in sp_char(mu, r).items():
                rec[k] = rec.get(k, 0) + a * v
        rec = {k: v for k, v in rec.items() if v != 0}
        ok = (rec == {k: v for k, v in Phi.items()})
        malo0 += (not ok)
        vals = sorted({abs(int(a)) for a in A.values()})
        for a in A.values():
            hist[abs(int(a))] += 1
        print("  %-30s | %3d | %-34s | %s"
              % (str(b), len(A), str(vals)[:34], "ok" if ok else "*** FALLA ***"))
        sys.stdout.flush()
    print("")
    print("  histograma de |A_mu| sobre las %d formas: %s"
          % (n, dict(sorted(hist.items()))))
    mayores = sum(v for k, v in hist.items() if k > 1)
    print("  pesos con |A_mu| >= 2: %d de %d" % (mayores, sum(hist.values())))
    if mayores == 0:
        print("  -> |A_mu| <= 1 EN TODO PESO.  Es lo que la estructura de producto de Yacobi predice.")
    else:
        print("  -> HAY pesos con |A_mu| >= 2: el acotamiento es EXCLUSIVO del peso maximo, y la ruta")
        print("     del producto de SL_2 no da la conjetura tal cual.")
    print("  C0 reconstruccion: %s" % ("PASA en %d/%d" % (n - malo0, n) if malo0 == 0 else "*** FALLA en %d ***" % malo0))
    RES["t%d" % t] = {"hist": {str(k): int(v) for k, v in hist.items()}, "n": int(n),
                      "mayores": int(mayores)}
    sys.stdout.flush()

json.dump(RES, open("all_A_mu_DUMP.json", "w"), indent=1)
print("")
print("=" * 120)
print("DONE")
