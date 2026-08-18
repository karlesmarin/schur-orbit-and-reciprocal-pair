# -*- coding: utf-8 -*-
# LOS ETA QUE SOBREVIVEN EN EL PESO SUPERIOR, LADO IMPAR.   16 de agosto de 2026.
#
# DE DONDE SALE.  Vuelta 22, P24 y P26.  El corrige nuestra frase "la cancelacion entra en el paso
# del determinante y solo ahi": es demasiado fuerte.  Lo que entra ahi es branching ORDINARIO contra
# TWINED; pero la suma final  A_mu = sum_eta B_{eta,mu} tau^B(eta)  con tau^B en {0,+-1} PUEDE
# cancelar igual.  Y propone la medida que lo decide:
#
#     en mu_max^+ :   #{ eta : B_{eta,mu_max} != 0  y  tau^B_t(eta) != 0 },  con multiplicidades.
#
#   * si sale siempre 1 x 1  -> el Unit Theorem impar cae por extremal branching + NPP, sin
#     cancelacion de tableaux, y seria MUCHO mas facil que el par;
#   * si salen varios y todos menos uno cancelan -> el impar tambien necesita atomos;
#   * si salen varios con patron rigido -> hay que buscar un teorema PRV/extremal para el par
#     simetrico SO_{2R'+1} -> SO_{2m'+1} x SO_{2r}.
#
# COMPARACION QUE HAY QUE TENER DELANTE.  En el PAR, en la forma mayor, el mismo coeficiente
# necesitaba 25 terminos de hasta 798 con sumas parciales de 455.  Si el impar sale con uno o dos,
# la dicotomia se ve tambien en la aritmetica.
#
# LO QUE SE MIDE, por forma
#   A1  el numero de PARES (Lambda, eta) que contribuyen a mu_max^+ con tau != 0.
#   A2  el numero de eta DISTINTOS que sobreviven, y sus multiplicidades  sum_Lambda a_Lambda B.
#   A3  los terminos individuales  a_Lambda * B * tau  y su mayor valor absoluto, y la mayor suma
#       parcial acumulada -- que es lo que mide si hay cancelacion de verdad.
#   A4  el reparto de signos.
#
# CONTROLES
#   C0  FATAL.  la suma de los terminos tiene que dar el A^D_{mu_max} conocido por la ruta
#       bialternante independiente.
#   C1  se imprime la forma entera, no un agregado: un promedio escondería el caso interesante.
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage odd_atoms.sage

import json
import sys
from collections import defaultdict


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


_BR = {}
def branch_BD(Rp, mp, rr, Lam):
    key = (Rp, mp, rr, tuple(int(v) for v in Lam))
    if key not in _BR:
        W = WeylCharacterRing("B%d" % Rp)
        X = WeylCharacterRing("B%dxD%d" % (mp, rr))
        br = branching_rule("B%d" % Rp, "B%dxD%d" % (mp, rr), "orthogonal_sum")
        el = W(W.space().from_vector(vector([Integer(v) for v in Lam])))
        d = {}
        for wt, c in el.branch(X, rule=br).monomial_coefficients().items():
            v = [int(u) for u in wt.to_vector()]
            d[(tuple(v[:mp]), tuple(v[mp:]))] = int(c)
        _BR[key] = d
    return _BR[key]


def tau_B(mp, eta, tt):
    K = CyclotomicField(tt)
    z = K.gen()
    s = K(0)
    for wt, mult in car("B", mp, eta).items():
        s += mult * z ** (sum((i + 1) * wt[i] for i in range(mp)) % tt)
    return QQ(s) if s in QQ else None


def betas(N, tope):
    return [tuple(sorted(c, reverse=True)) + (0,) for c in Combinations(range(1, tope + 1), N - 1)]


CASOS = [(3, 2, 9), (5, 2, 10)]

print("=" * 126)
print("LOS ETA QUE SOBREVIVEN EN mu_max^+   --   lado IMPAR   (la medida que pide su P26)")
print("=" * 126)
RES = []
for (t, r, tope) in CASOS:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    N = t + 2 * r
    print("")
    print("  t=%d  m'=%d  r=%d  R'=%d  N=%d" % (t, mp, r, Rp, N))
    print("  beta                       | mu_max^+        |  A | #(Lam,eta) | #eta distintos | terminos"
          " | mayor |t| | mayor parcial | C0")
    print("  " + "-" * 122)
    hist_eta = defaultdict(lambda: 0)
    hist_term = defaultdict(lambda: 0)
    n_forma = 0
    for b in betas(N, tope):
        P1 = phi_bialt(b, t, r)
        if P1 in (None, "NO-POL", "NO-RAC") or not P1:
            continue
        A1, rest1 = pelar(P1, "D", r)
        A1 = {k: QQ(v) for k, v in A1.items() if v != 0}
        if rest1 or not A1:
            continue
        Ap = {}
        for mu, c in A1.items():
            Ap[tuple(list(mu[:-1]) + [abs(mu[-1])])] = c
        S = list(Ap)
        maxi = [m for m in S if not any(n2 != m and all(sum(n2[:k + 1]) >= sum(m[:k + 1])
                                                        for k in range(len(m))) for n2 in S)]
        if len(maxi) != 1:
            continue
        mm = maxi[0]
        P2 = phi_bialt(b, 1, Rp)
        aB, restB = pelar(P2, "B", Rp)
        if restB:
            continue
        pares = 0
        etas = defaultdict(lambda: 0)
        term = []
        for Lam, a in aB.items():
            for (eta, mu), c in branch_BD(Rp, mp, r, Lam).items():
                if tuple(list(mu[:-1]) + [abs(mu[-1])]) != mm:
                    continue
                tv = tau_B(mp, eta, t)
                if tv is None or tv == 0:
                    continue
                pares += 1
                etas[eta] += int(a * c)
                term.append(int(a * c * tv))
        # las contribuciones de mu y mu* se cuentan las dos cuando mu_r != 0
        fac = 2 if mm[-1] != 0 else 1
        suma = sum(term)
        c0 = (QQ(suma) == fac * Ap[mm])
        par = 0
        mx = 0
        for v in term:
            par += v
            mx = max(mx, abs(par))
        n_forma += 1
        hist_eta[len([e for e, v in etas.items() if v != 0])] += 1
        hist_term[len(term)] += 1
        if n_forma <= 12:
            print("  %-26s | %-15s | %2s | %10d | %14d | %8d | %9d | %13d | %s"
                  % (str(b), str(mm), str(Ap[mm]), pares,
                     len([e for e, v in etas.items() if v != 0]), len(term),
                     max([abs(v) for v in term]) if term else 0, mx, "ok" if c0 else "FALLA"))
            sys.stdout.flush()
        RES.append({"t": int(t), "r": int(r), "beta": [int(v) for v in b],
                    "mu_max": [int(v) for v in mm], "A": int(Ap[mm]),
                    "n_pares": int(pares),
                    "n_eta": int(len([e for e, v in etas.items() if v != 0])),
                    "n_terminos": int(len(term)),
                    "mayor_termino": int(max([abs(v) for v in term])) if term else 0,
                    "mayor_parcial": int(mx), "C0": bool(c0)})
    print("  ...")
    print("  FORMAS: %d | histograma de #eta distintos: %s | histograma de #terminos: %s"
          % (n_forma, dict(sorted(hist_eta.items())), dict(sorted(hist_term.items()))))
    print("  C0 pasa en: %d de %d" % (sum(1 for d in RES if d["t"] == t and d["C0"]),
                                      sum(1 for d in RES if d["t"] == t)))
    sys.stdout.flush()

print("")
print("=" * 126)
print("  LECTURA, escrita ANTES de correr:")
print("   * histograma concentrado en 1 -> el Unit Theorem impar cae por extremal branching + NPP.")
print("   * varios eta pero terminos pequeños -> cancelacion suave, y aun asi mucho mas facil que el par.")
print("   * terminos grandes como en el par -> la dicotomia NO alcanza a la aritmetica final, y hay")
print("     que decirlo asi en el paper.")
json.dump(RES, open("odd_atoms_DUMP.json", "w"), indent=1)
print("")
print("=" * 126)
print("DONE")
