# -*- coding: utf-8 -*-
# INTERROGATORIO DE LA QUIRALIDAD.  .Que es simetrico y que antisimetrico?
#
# La duda: (L2) se enuncia en mu^+ (PLEGADO) y prop:transversal indexa nu por pesos de D_r SIN
# plegar, contando las dos quiralidades como dos puntos.  Si c(Lambda,mu) = c(Lambda,mu*), el 2 de
# (ii) desaparece al plegar y hay que decir en que soporte se cuenta.
#
# Prediccion a falsar, en este orden:
#   Q1  nu(Lambda, mu*) == -nu(Lambda, mu)      (ANTIsimetrica)
#   Q2  Delta_t es antisimetrico bajo quiralidad
#   Q3  luego  c(Lambda, mu*) == +c(Lambda, mu)  (SIMETRICA)
#
# Si Q3 falla, el paper esta comparando cosas distintas en (L2) y hay que rehacer el enunciado.
# Si Q1 falla, mi construccion de nu tiene el signo de la quiralidad mal puesto.

import itertools
from collections import defaultdict


def plegar(v, t):
    v %= t
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
    return s


def jacobi(a, n):
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
    return jacobi((-2) % t, t) ** ((t + 3) // 2) * (1 if (mp * (mp - 1) // 2) % 2 == 0 else -1)


def delta_dec(A, t, mp):
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
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    V = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
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


def dividir(nu, t, r, tope=20000):
    P = dict(nu)
    c = {}
    for _ in range(tope):
        P = {k: v for k, v in P.items() if v != 0}
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


print("=" * 96)
print("INTERROGATORIO DE LA QUIRALIDAD")
print("=" * 96)
q1 = q1n = q3 = q3n = 0
solo_una = 0
for (t, r, cota) in [(3, 2, 6), (5, 2, 4), (7, 2, 3), (3, 3, 4)]:
    Rp = (t - 1) // 2 + r
    for Lam in itertools.product(range(cota + 1), repeat=Rp):
        if any(Lam[i] < Lam[i + 1] for i in range(Rp - 1)):
            continue
        nu = nu_de(list(Lam), t, r)
        if not nu:
            continue
        # Q1: nu(mu*) == -nu(mu)
        for x, v in nu.items():
            xs = tuple(list(x[:-1]) + [-x[-1]])
            q1n += 1
            if x[-1] == 0:
                solo_una += 1
                continue
            q1 += 1 if nu.get(xs, 0) == -v else 0
        # Q3: c(mu*) == +c(mu)
        c, resto = dividir(nu, t, r)
        if resto:
            continue
        for x, v in c.items():
            xs = tuple(list(x[:-1]) + [-x[-1]])
            if x[-1] == 0:
                continue          # x = x*, no dice nada: NO cuenta en el denominador
            q3n += 1
            q3 += 1 if c.get(xs, 0) == v else 0
print("")
print("  Q1  nu(mu*) == -nu(mu)   (ANTIsimetrica) : %d de %d   (con ultima coord 0: %d)"
      % (q1, q1n, solo_una))
print("  Q3  c(mu*)  == +c(mu)    (SIMETRICA)     : %d de %d" % (q3, q3n))
print("")
print("  LECTURA: si Q1 y Q3 salen, nu es antisimetrica y c simetrica -- y entonces el '2' de")
print("  prop:transversal(ii) cuenta el soporte SIN PLEGAR, mientras (L2) vive en el plegado.")
print("  Hay que decirlo en el enunciado, o se estan comparando dos soportes distintos.")
print("=" * 96)
print("DONE")
