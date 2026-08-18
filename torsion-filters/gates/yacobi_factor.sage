# -*- coding: utf-8 -*-
# ¿ES NUESTRO Sp_2 UNO DE LOS FACTORES SL_2 DE YACOBI, O NO?   15 de agosto de 2026.
#
# POR QUE.  all_A_mu.sage midio que a t=4 (m=1) hay NUEVE pesos con |A_mu| = 2, lo que refuta la
# ruta "M_mu es un producto tensorial de SL_2-modulos => el caracter en el elemento de torsion es un
# producto de 0,+-1 => |A_mu| <= 1".  Pero antes de devolverle ese negativo hay que separar dos
# cosas que se parecen y no son lo mismo:
#
#   (a) la idea es correcta y el enunciado |A_mu|<=1 simplemente no se sigue;  o
#   (b) nuestro Sp_2 NO es uno de los factores SL_2 de Yacobi -- por ejemplo es una diagonal, o el
#       espacio de multiplicidad relevante es el de OTRA ramificacion -- y entonces el caracter
#       nunca fue un producto, y el negativo mide mi aplicacion y no su propuesta.
#
# Mandar (b) como si fuera (a) seria devolverle un negativo falso.
#
# LO QUE DISTINGUE (a) DE (b), Y ES MEDIBLE.  Si nuestro Sp_2 fuera un factor tensorial, entonces
# para cada mu el espacio de multiplicidad M_mu = sum_eta B_{eta,mu} [sp_eta] seria, como
# Sp_2-modulo VIRTUAL, un multiplo entero de un unico irreducible por un entero d >= 1:
#
#       M_mu  =  d * [V_k]   para un solo k        <=>   B_{.,mu} tiene UN SOLO eta con B != 0
#
# porque un producto tensorial V_{d_1} (x) ... (x) V_{d_R} restringido a UN factor es
# (dim del resto) copias de ese factor.  Asi que el test es directo y no necesita a Yacobi:
#
#       ¿cuantos eta distintos aparecen en el bloque de cada mu?
#
# Si es siempre UNO, el Sp_2 se comporta como factor y (a) es la lectura.  Si son varios, el
# espacio de multiplicidad NO es un solo irreducible del factor, el caracter no es un producto, y
# (b) es la lectura -- el negativo seria sobre mi aplicacion.
#
# CONTROLES
#   C0  FATAL.  La reconstruccion sum_mu A_mu sp_mu = Phi_{t,r}.
#   C1  se imprime el histograma del NUMERO de eta por mu, no un promedio.
#   C2  se separan los mu maximales del resto: la conjetura es sobre la punta, y puede que la punta
#       se comporte distinto del cuerpo -- que es justo lo que all_A_mu sugiere.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage yacobi_factor.sage

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


def pelar(P, m, r, tope=9000):
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
        out[(tuple(top[:m]), tuple(top[m:]))] = out.get((tuple(top[:m]), tuple(top[m:])), 0) + B
        a, b = sp_char(tuple(top[:m]), m), sp_char(tuple(top[m:]), r)
        for e1, c1 in a.items():
            for e2, c2 in b.items():
                k = e1 + e2
                nv = P.get(k, 0) - B * c1 * c2
                if nv == 0:
                    P.pop(k, None)
                else:
                    P[k] = nv
    return out, P


def tau4(eta):
    k = eta[0]
    return 0 if k % 2 else (-1) ** (k // 2)


t, r, m, R = 4, 2, 1, 3
CASOS = [(18, 17, 11, 8, 7, 6, 1, 0), (10, 9, 7, 4, 3, 2, 1, 0), (14, 13, 11, 4, 3, 2, 1, 0),
         (12, 11, 10, 5, 3, 2, 1, 0), (14, 13, 11, 8, 3, 2, 1, 0), (12, 11, 10, 9, 7, 2, 1, 0),
         (13, 9, 8, 7, 5, 4, 2, 0)]

print("=" * 118)
print("¿CUANTOS eta POR mu?  --  t=4, m=1: si fuera UNO, el Sp_2 actua como factor tensorial")
print("=" * 118)
print("")
print("  %-30s | #mu | eta por mu (histograma) | en el mu MAXIMAL | |A_max| | C0" % "beta")
print("  " + "-" * 112)
hist_todos, hist_max = defaultdict(int), defaultdict(int)
malo = n = 0
for b in CASOS:
    Psi = phi_bialternante(b, 2, R)
    Phi = phi_bialternante(b, t, r)
    if Psi in (None, "NO-POL") or Phi in (None, "NO-POL"):
        continue
    n += 1
    B, resto = pelar(Psi, m, r)
    B = {k: v for k, v in B.items() if v != 0}
    porm = defaultdict(list)
    for (eta, mu), bb in B.items():
        porm[mu].append((eta, bb))
    A = {}
    for mu, L in porm.items():
        a = sum(bb * tau4(eta) for eta, bb in L)
        if a:
            A[mu] = a
    rec = {}
    for mu, a in A.items():
        for k, v in sp_char(mu, r).items():
            rec[k] = rec.get(k, 0) + a * v
    rec = {k: v for k, v in rec.items() if v != 0}
    ok = (rec == {k: v for k, v in Phi.items()})
    malo += (not ok)
    S = list(A)
    maxi = [mu for mu in S if not any(nu != mu and all(sum(nu[:k + 1]) >= sum(mu[:k + 1])
                                                       for k in range(r)) for nu in S)]
    h = defaultdict(int)
    for mu in A:
        h[len(porm[mu])] += 1
        hist_todos[len(porm[mu])] += 1
    nmax = [len(porm[mu]) for mu in maxi]
    for v in nmax:
        hist_max[v] += 1
    print("  %-30s | %3d | %-23s | %-16s | %-7s | %s"
          % (str(b), len(A), str(dict(sorted(h.items())))[:23], str(nmax),
             str([int(A[mu]) for mu in maxi]), "ok" if ok else "*** FALLA ***"))
    sys.stdout.flush()

print("")
print("  histograma de #eta por mu, TODOS los mu   : %s" % dict(sorted(hist_todos.items())))
print("  histograma de #eta por mu, solo MAXIMALES : %s" % dict(sorted(hist_max.items())))
print("  C0 reconstruccion: %s" % ("PASA %d/%d" % (n - malo, n) if malo == 0 else "*** FALLA %d ***" % malo))
print("")
print("  LECTURA, escrita antes de correr:")
print("   * si #eta = 1 siempre  -> el Sp_2 se comporta como factor tensorial, el caracter SI es un")
print("     producto, y el |A_mu|=2 refuta la ruta de verdad.")
print("   * si #eta > 1          -> el espacio de multiplicidad no es un solo irreducible del factor,")
print("     el caracter nunca fue un producto, y el negativo es sobre MI aplicacion, no sobre su idea.")
print("")
print("=" * 118)
print("DONE")
