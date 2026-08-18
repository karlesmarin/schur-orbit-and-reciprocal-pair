# -*- coding: utf-8 -*-
# .Por que falla E7 en el par?  El argumento del impar no usa la paridad: v es estrictamente
# decreciente, S_min toma el valor MINIMO de v en cada clase, luego minimiza sum(v|_S) y maximiza el
# complemento.  Si en el par falla, o el argumento tiene un agujero o mi S_min no es lo que creo.

import itertools


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


def delta_C(a, t, m):
    cl, ep = [], []
    for v in a:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, m + 1)):
        return 0
    s = sgn_perm([m - c for c in cl])
    for e in ep:
        s *= e
    return int(s)


for (t, r, cota) in [(6, 2, 3)]:
    m = (t - 2) // 2
    R = m + r
    print("t=%d r=%d m=%d R=%d" % (t, r, m, R))
    malos = 0
    for Lam in itertools.product(range(cota + 1), repeat=R):
        if any(Lam[i] < Lam[i + 1] for i in range(R - 1)):
            continue
        v = [Lam[i] + R - i for i in range(R)]
        d = {}
        for i, x in enumerate(v):
            d.setdefault(plegar(x, t)[0], []).append(i)
        if any(not d.get(j) for j in range(1, m + 1)):
            continue
        nu = {}
        for S in itertools.combinations(range(R), m):
            Sc = [i for i in range(R) if i not in S]
            A = sorted([v[i] for i in S], reverse=True)
            dv = delta_C(A, t, m)
            if not dv:
                continue
            libre = tuple(sorted([v[i] for i in Sc], reverse=True))
            orden = sorted(S, key=lambda i: -v[i]) + sorted(Sc, key=lambda i: -v[i])
            sg = sgn_perm([orden.index(i) for i in range(R)])
            nu[libre] = nu.get(libre, 0) + sg * dv
        nu = {k: val for k, val in nu.items() if val}
        if not nu:
            continue
        top = max(nu, key=lambda k: (sum(k), k))
        S_min = frozenset(max(d[j]) for j in range(1, m + 1))
        Sc = [i for i in range(R) if i not in S_min]
        pred = tuple(sorted([v[i] for i in Sc], reverse=True))
        if pred != top:
            malos += 1
            if malos <= 3:
                sumas = sorted(((sum(k), k) for k in nu), reverse=True)
                print("  FALLO Lambda=%s  v=%s" % (list(Lam), v))
                print("        clases: %s" % {k: val for k, val in sorted(d.items())})
                print("        S_min=%s -> pred=%s (suma %d)" % (sorted(S_min), pred, sum(pred)))
                print("        top real =%s (suma %d)" % (top, sum(top)))
                print("        alturas de supp nu: %s" % sumas[:4])
                print("        .es pred un elemento de supp nu?  %s" % (pred in nu))
    print("  fallos:", malos)
