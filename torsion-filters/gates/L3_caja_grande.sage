# -*- coding: utf-8 -*-
# (L3) DONDE PUEDE FALLAR.   16 de agosto de 2026.
#
# POR QUE SE REHACE.  (L3) dice que el unico Lambda con c(Lambda, mu_max) != 0 tiene a^B_Lambda = 1.
# Medida en beta_i <= 9/10/9 salio 132/132 -- pero ALLI TODOS los a^B_Lambda valen 1 (476 de 476),
# asi que el enunciado no puede ser falso en esa caja y la medida no dice nada.  En beta_i <= 13
# solo el 47% de los a^B valen 1, con coeficientes hasta 23.  Aqui se mide (L3) donde el fallo
# es POSIBLE.
#
# EL ATAJO QUE LO HACE VIABLE, y su control.  Ya no hace falta el branching por Lambda:
#   nu(Lambda,.)  es combinatoria pura           (prop:transversal)
#   c(Lambda,.) = nu / Delta_t  por sustitucion hacia atras   (verificado 245/245 en gkrs_L1)
# Eso quita la parte cara.  Pero el atajo esta verificado a otra escala, asi que:
#   C1  CONTROL FATAL: en las primeras formas se compara la c combinatoria con la c por BRANCHING.
#       Si discrepan, el atajo no vale a esta escala y el resto de la corrida no se mira.
#
# LO QUE SE MIDE
#   V0  DISPONIBILIDAD: .cuantas formas tienen ALGUN Lambda con a^B > 1?  Si son pocas, la caja
#       sigue sin testar (L3) y hay que subirla mas.  Va primero.
#   L3  a^B = 1 para el Lambda que aporta.
#   L2  ese Lambda es unico.
#   L3b reparto de a^B del que aporta, y del maximo disponible en su forma.
#   X   .hay correlacion?  se tabula (a^B del que aporta, max a^B de la forma).
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage L3_caja_grande.sage

import json
import sys
import itertools
from collections import Counter, defaultdict

CASOS = [(3, 2, 12), (5, 2, 12)]
N_CONTROL = 12          # formas en las que se cruza el atajo contra el branching


def plegar(v, t):
    v = int(v) % t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def sgn_perm(perm):
    n, s, visto = len(perm), 1, [False] * len(perm)
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
    return int(s)


def jacobi_i(a, n):
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


def eps_t(t, mp):
    e = jacobi_i((-2) % t, t) ** ((t + 3) // 2)
    return int(e * (1 if (mp * (mp - 1) // 2) % 2 == 0 else -1))


def delta_dec(A, t, mp):
    """delta en el convenio de lem:T: clases contra el orden DECRECIENTE."""
    cl, ep = [], []
    for v in A:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, mp + 1)):
        return 0
    s = sgn_perm([mp - c for c in cl])
    for e in ep:
        s *= e
    return int(s)


def enderezar_D(x):
    r = len(x)
    a = [abs(int(v)) for v in x]
    if len(set(a)) != r:
        return None
    idx = sorted(range(r), key=lambda i: -a[i])
    s = sgn_perm(list(idx))
    y = [int(x[i]) for i in idx]
    neg = sum(1 for v in y if v < 0)
    cero = any(v == 0 for v in y)
    y = [abs(v) for v in y]
    if neg % 2 == 1 and not cero:
        y[-1] = -y[-1]
    return (tuple(y), s)


def desplazar(x, paso, r):
    out = defaultdict(lambda: 0)
    for eps in itertools.product((1, -1), repeat=r):
        sg = 1
        for e in eps:
            sg *= e
        e2 = enderezar_D(tuple(int(x[j]) + paso * eps[j] for j in range(r)))
        if e2 is None:
            continue
        out[e2[0]] += sg * e2[1]
    return {k: v for k, v in out.items() if v != 0}


def cabeza(d):
    return max(d, key=lambda k: (sum(k), k))


def nu_de(Lam, t, r):
    """nu por transversales.  Indices x = V restringido al bloque libre."""
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    V = [2 * int(Lam[i]) + 2 * (Rp - i) - 1 for i in range(Rp)]
    d = {}
    for i, x in enumerate(V):
        d.setdefault(plegar(x, t)[0], []).append(i)
    if any(not d.get(j) for j in range(1, mp + 1)):
        return {}
    E = eps_t(t, mp)
    out = {}
    for pick in itertools.product(*[d[j] for j in range(1, mp + 1)]):
        S = frozenset(pick)
        if len(S) != mp:
            continue
        Sc = [i for i in range(Rp) if i not in S]
        A = sorted([V[i] for i in S], reverse=True)
        dv = delta_dec(A, t, mp)
        if not dv:
            continue
        for qui in (1, -1):
            libre = sorted([V[i] for i in Sc], reverse=True)
            libre[-1] *= qui
            orden = sorted(S, key=lambda i: -V[i]) + sorted(Sc, key=lambda i: -V[i])
            sg = sgn_perm([orden.index(i) for i in range(Rp)])
            if qui == -1:
                sg = -sg
            out[tuple(libre)] = out.get(tuple(libre), 0) + sg * E * dv
    return {k: v for k, v in out.items() if v != 0}


def dividir(nu, t, r, tope=40000):
    P = {k: int(v) for k, v in nu.items() if v != 0}
    c = {}
    for _ in range(tope):
        if not P:
            return c, {}
        y = cabeza(P)
        cand = None
        for eps in itertools.product((1, -1), repeat=r):
            e = enderezar_D(tuple(int(y[j]) - t * eps[j] for j in range(r)))
            if e is None:
                continue
            D = desplazar(e[0], t, r)
            if D and cabeza(D) == y:
                cand = (e[0], D)
                break
        if cand is None:
            return c, P
        x, D = cand
        if P[y] % D[y] != 0:
            return c, P
        cv = P[y] // D[y]
        c[x] = c.get(x, 0) + cv
        for k, v in D.items():
            nv = P.get(k, 0) - cv * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return c, P


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


def pelar(P, typ, rk, tope=40000):
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


def tauB_car(mp, eta, t):
    K = CyclotomicField(t)
    z = K.gen()
    s = K(0)
    for wt, mult in car("B", mp, eta).items():
        s += mult * z ** (sum((i + 1) * int(wt[i]) for i in range(mp)) % t)
    return int(QQ(s)) if s in QQ else None


def domina(a, b):
    return all(sum(a[:k + 1]) >= sum(b[:k + 1]) for k in range(len(b)))


print("=" * 112)
print("(L3) EN LA CAJA DONDE PUEDE FALLAR")
print("=" * 112)
sys.stdout.flush()

RES = []
for (t, r, tope) in CASOS:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    N = t + 2 * r
    rhoD2 = [2 * (r - j) - 2 for j in range(r)]
    n = l3 = l2 = v0 = c1_ok = c1_n = 0
    ap = Counter()
    cruce = Counter()
    fallos = []
    for b in [tuple(sorted(c, reverse=True)) + (0,) for c in Combinations(range(1, tope + 1), N - 1)]:
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
        maxi = [m for m in S if not any(nn != m and domina(nn, m) for nn in S)]
        if len(maxi) != 1:
            continue
        mm = maxi[0]
        P2 = phi(b, 1, Rp)
        aB, rest2 = pelar(P2, "B", Rp)
        if rest2 or not aB:
            continue
        aB = {k: int(v) for k, v in aB.items() if v != 0}
        n += 1
        maxa = max(abs(v) for v in aB.values())
        if maxa > 1:
            v0 += 1

        # c(Lambda, mu_max) por el atajo combinatorio
        cL = {}
        for Lam in aB:
            nu = nu_de(Lam, t, r)
            if not nu:
                continue
            cc, resto = dividir(nu, t, r)
            if resto:
                continue
            val = 0
            for x, cv in cc.items():
                mu = tuple((x[j] - rhoD2[j]) // 2 for j in range(r))
                if tuple(list(mu[:-1]) + [abs(mu[-1])]) == mm:
                    val += cv
            if val:
                cL[Lam] = val

        # C1: cruce contra el branching en las primeras formas
        if c1_n < N_CONTROL:
            c1_n += 1
            cB = {}
            for Lam in aB:
                s = 0
                for (eta, mu), cc2 in branch(Rp, mp, r, Lam).items():
                    mup = tuple(list(mu[:-1]) + [abs(mu[-1])])
                    if mup != mm:
                        continue
                    tv = tauB_car(mp, eta, t)
                    if tv:
                        s += int(cc2) * int(tv)
                if s:
                    cB[Lam] = s
            if set(cB) == set(cL):
                c1_ok += 1
            else:
                fallos.append({"C1": "el atajo discrepa del branching",
                               "beta": [int(x) for x in b]})

        l2 += 1 if len(cL) == 1 else 0
        if len(cL) == 1:
            Lst = list(cL)[0]
            a_star = abs(int(aB[Lst]))
            ap[a_star] += 1
            l3 += 1 if a_star == 1 else 0
            cruce[(a_star, maxa)] += 1
            if a_star != 1 and len(fallos) < 6:
                fallos.append({"L3 FALLA": True, "beta": [int(x) for x in b],
                               "Lambda": [int(x) for x in Lst], "a": a_star,
                               "max_a_de_la_forma": int(maxa)})
    print("")
    print("  t=%d r=%d  beta_i <= %d   (%d formas)" % (t, r, tope, n))
    print("     V0  formas con ALGUN a^B > 1 (la caja testa)  : %3d de %3d = %.0f%%"
          % (v0, n, 100.0 * v0 / n if n else 0))
    print("     C1  el atajo == el branching                  : %3d de %3d   <== FATAL"
          % (c1_ok, c1_n))
    print("     L2  el Lambda que aporta es unico             : %3d de %3d" % (l2, n))
    print("     L3  y tiene a^B = 1                           : %3d de %3d" % (l3, l2))
    print("         reparto de a^B del que aporta             : %s" % dict(sorted(ap.items())))
    print("         cruce (a^B del que aporta, max a^B forma) : %s"
          % dict(sorted(cruce.items(), key=lambda kv: kv[0])[:12]))
    if fallos:
        print("     !! %s" % json.dumps(fallos[:3]))
    sys.stdout.flush()
    RES.append({"t": int(t), "r": int(r), "tope": int(tope), "n": int(n), "V0": int(v0),
                "C1": [int(c1_ok), int(c1_n)], "L2": int(l2), "L3": int(l3),
                "reparto_aporta": {str(k): int(v) for k, v in ap.items()},
                "fallos": fallos})

json.dump(RES, open("L3_caja_grande_DUMP.json", "w"), indent=1)
print("")
print("=" * 112)
print("  LECTURA, escrita ANTES de correr:")
print("   * si V0 es alto y L3 sale completa, (L3) por fin esta testada y sobrevive.")
print("   * si V0 es bajo, la caja SIGUE sin testarla y hay que subirla otra vez.")
print("   * si L3 falla, el testigo dice que Lambda con a^B > 1 alcanza el extremo, y eso mata")
print("     (L3) tal como esta enunciada.")
print("=" * 112)
print("DONE")
