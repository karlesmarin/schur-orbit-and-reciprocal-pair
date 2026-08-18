# -*- coding: utf-8 -*-
# LA ESTRUCTURA EXTREMAL DEL LADO IMPAR: ¿por que sale +-1?   16 de agosto de 2026.
#
# DE DONDE SALE.  Su vuelta 26 propone demostrar la MITAD IMPAR de la primitividad extremal, donde
# hay branching genuino, elemento principal y 123/123 casos medidos.  Pero antes de intentar una
# prueba hay que saber QUE hay que probar, y lo que sabemos no basta:
#
#   * 93 de 123 formas tienen UN SOLO eta en mu_max -- pero 30 tienen dos o tres;
#   * y "un solo eta" no es "un solo termino": varios Lambda pueden aportar el mismo eta.
#
# Asi que se abre la descomposicion ENTERA en mu_max, pareja a pareja:
#
#       A^D_{mu_max} = sum_{Lambda, eta}  a^B_Lambda . B^odd_{Lambda; eta, mu_max} . tau^B(eta)
#
# y se busca que distingue al sumando que sobrevive.  Las hipotesis candidatas, escritas ANTES:
#
#   H1  hay un Lambda MAXIMO (el de la componente extremal de Littlewood, con a = 1) y su pareja
#       extremal aporta +-1; los demas se cancelan entre si.
#   H2  no hay un Lambda distinguido, pero el sumando extremal tiene B = 1 y los otros vienen en
#       parejas de signo opuesto.
#   H3  ninguna de las dos, y la cancelacion es colectiva tambien aqui -- en cuyo caso la mitad
#       impar NO es mas facil y hay que decirlo.
#
# LO QUE SE MIDE, forma a forma
#   E1  la tabla (Lambda, a_Lambda, eta, B, tau, termino) en mu_max, ordenada por Lambda.
#   E2  si hay un Lambda maximo en dominancia entre los que contribuyen, y cual es su a_Lambda.
#   E3  el termino de ESE Lambda, contra la suma de todos los demas.
#   E4  a^B_Lambda de la componente extremal de Littlewood: ¿vale 1?
#
# CONTROLES
#   C0  FATAL: la suma reproduce A^D_{mu_max} calculado por la ruta bialternante.
#   C1  se imprime la tabla entera, no un resumen.
#   C2  n impreso siempre.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage odd_extremal.sage

import json
import sys
from collections import defaultdict


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


def pelar(P, typ, rk, tope=9000):
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
def branch(Rp, mp, rr, Lam):
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


def tauB(mp, eta, t):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car("B", mp, eta).items():
        s += mult * z ** (sum((i + 1) * wt[i] for i in range(mp)) % t)
    return QQ(s) if s in QQ else None


def domina(a, b):
    return all(sum(a[:k + 1]) >= sum(b[:k + 1]) for k in range(len(b)))


CASOS = [(3, 2, [(12, 10, 7, 5, 3, 2, 0), (11, 8, 6, 4, 3, 1, 0), (10, 9, 6, 4, 2, 1, 0),
                 (9, 8, 7, 5, 3, 2, 0), (12, 11, 7, 5, 3, 2, 0)]),
         (5, 2, [(12, 10, 8, 7, 5, 4, 2, 1, 0), (11, 10, 9, 6, 5, 3, 2, 1, 0)])]

print("=" * 122)
print("LA DESCOMPOSICION EXTREMAL DEL LADO IMPAR")
print("=" * 122)

RES = []
for (t, r, betas) in CASOS:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    N = t + 2 * r
    delta = list(range(N - 1, -1, -1))
    for b in betas:
        P1 = phi(b, t, r)
        if not P1:
            continue
        A1, rest = pelar(P1, "D", r)
        A1 = {k: QQ(v) for k, v in A1.items() if v != 0}
        if rest or not A1:
            continue
        Ap = {}
        for mu, c in A1.items():
            Ap[tuple(list(mu[:-1]) + [abs(mu[-1])])] = c
        S = list(Ap)
        maxi = [m for m in S if not any(n != m and domina(n, m) for n in S)]
        if len(maxi) != 1:
            continue
        mm = maxi[0]
        P2 = phi(b, 1, Rp)
        aB, rest2 = pelar(P2, "B", Rp)
        if rest2:
            continue
        lam = tuple(v for v in (b[i] - delta[i] for i in range(N)) if v > 0)
        filas = []
        for Lam, a in aB.items():
            for (eta, mu), c in branch(Rp, mp, r, Lam).items():
                if tuple(list(mu[:-1]) + [abs(mu[-1])]) != mm:
                    continue
                tv = tauB(mp, eta, t)
                if tv is None or tv == 0:
                    continue
                filas.append((Lam, int(a), eta, int(c), int(tv), int(a) * int(c) * int(tv)))
        fac = 2 if mm[-1] != 0 else 1
        suma = sum(f[5] for f in filas)
        # E2: hay un Lambda maximo entre los que contribuyen?
        Lams = sorted(set(f[0] for f in filas), reverse=True)
        maxL = [L for L in Lams if not any(M != L and domina(M, L) for M in Lams)]
        print("")
        print("  t=%d  beta=%s  lambda=%s" % (t, str(b), str(lam)))
        print("     mu_max^+ = %s   A = %s   (suma = %s, factor quiral %d)  %s"
              % (str(mm), str(Ap[mm]), str(suma), fac,
                 "C0 ok" if QQ(suma) == fac * Ap[mm] else "C0 FALLA"))
        print("     Lambda              a_L   eta        B   tau   termino")
        for f in sorted(filas, key=lambda x: (tuple(-v for v in x[0]), x[2])):
            marca = "  <-- Lambda maximo" if (f[0] in maxL and len(maxL) == 1) else ""
            print("     %-19s %4d   %-10s %3d %5d %8d%s"
                  % (str(f[0]), f[1], str(f[2]), f[3], f[4], f[5], marca))
        if len(maxL) == 1:
            tmax = sum(f[5] for f in filas if f[0] == maxL[0])
            tresto = suma - tmax
            print("     -> Lambda maximo unico %s : su aporte %d, el resto %d"
                  % (str(maxL[0]), tmax, tresto))
        else:
            print("     -> NO hay Lambda maximo unico: %s" % str(maxL))
        sys.stdout.flush()
        RES.append({"t": int(t), "beta": [int(v) for v in b], "mu_max": [int(v) for v in mm],
                    "A": int(Ap[mm]), "suma": int(suma), "n_filas": int(len(filas)),
                    "n_Lambda_max": int(len(maxL)),
                    "aporte_Lmax": int(sum(f[5] for f in filas if f[0] in maxL)) if len(maxL) == 1 else None,
                    "filas": [[list(map(int, f[0])), f[1], list(map(int, f[2])), f[3], f[4], f[5]]
                              for f in filas]})

print("")
print("=" * 122)
print("  LECTURA, escrita ANTES de correr:")
print("   * H1: si hay un Lambda maximo UNICO con a_L = 1 cuyo aporte es +-1 y el resto suma 0,")
print("     la mitad impar sale por un argumento extremal y se puede escribir.")
print("   * H2: si no, pero los demas vienen en parejas de signo opuesto, hay involucion local.")
print("   * H3: si la cancelacion es colectiva tambien aqui, el impar NO es mas facil y hay que")
print("     decirlo en el paper en vez de prometer una prueba.")
json.dump(RES, open("odd_extremal_DUMP.json", "w"), indent=1)
print("=" * 122)
print("DONE")
