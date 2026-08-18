# LA PREDICCION CORREGIDA: no UN caracter simplectico, sino una COMBINACION ENTERA NO NEGATIVA
# (que es lo que da el branching GL -> O seguido del plegado).  Aritmetica entera, sin sympy.
from itertools import permutations, product
from collections import Counter


def pmul(a, b):
    o = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = ea + eb
            o[e] = o.get(e, 0) + ca * cb
    return {e: c for e, c in o.items() if c}


def padd(a, b):
    o = dict(a)
    for e, c in b.items():
        o[e] = o.get(e, 0) + c
    return {e: c for e, c in o.items() if c}


def sgn(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def det(M):
    n = len(M)
    acc = {}
    for p in permutations(range(n)):
        t = {0: sgn(p)}
        for i in range(n):
            t = pmul(t, M[i][p[i]])
            if not t:
                break
        if t:
            acc = padd(acc, t)
    return acc


def A(expos):                       # det(u_j^{b_i}) con u = (1, -1, z, 1/z)
    return det([[{0: 1}, {0: (-1) ** (b % 2)}, {b: 1}, {-b: 1}] for b in expos])


def sp_car(mu):                     # caracter de Sp(2): z^mu + z^{mu-2} + ... + z^{-mu}
    return {e: 1 for e in range(mu, -mu - 1, -2)}


def cond(lam, t=2):
    N = len(lam)
    beta = [lam[i] + (N - 1 - i) for i in range(N)]
    cl = {}
    for b in beta:
        cl.setdefault(b % t, []).append(b)
    if len(cl) < t:
        return None
    E = [k for k in cl if len(cl[k]) >= 2]
    if not E:
        return None
    S = sorted({v for k in E for v in cl[k]})
    C = S[0] + S[-1]
    sols = [k for k in range(t) if (2 * k - C) % t == 0]
    return (set(C - v for v in S) == set(S)) and len(sols) == 2 and all(k in E for k in sols)


N = 4
Ad = A([N - 1 - i for i in range(N)])
cnt = Counter()
ej = []
for lam in product(range(7, -1, -1), repeat=N):
    if any(lam[i] < lam[i + 1] for i in range(N - 1)):
        continue
    Ab = A([lam[i] + (N - 1 - i) for i in range(N)])
    if not Ab:
        cnt[("CERO", cond(list(lam)))] += 1
        continue
    q = {}
    rem = dict(Ab)
    dmax = max(Ad)
    dc = Ad[dmax]
    while rem:
        m = max(rem)
        c = rem[m]
        if c % dc:
            q = None
            break
        k = m - dmax
        q[k] = q.get(k, 0) + c // dc
        f = c // dc
        for e, cc in Ad.items():
            rem[e + k] = rem.get(e + k, 0) - f * cc
        rem = {e: v for e, v in rem.items() if v}
    if q is None:
        cnt[("no divide", None)] += 1
        continue
    coef = {}
    r2 = dict(q)
    malo = False
    while r2:
        m = max(r2)
        if m < 0:
            malo = True
            break
        c = r2[m]
        coef[m] = coef.get(m, 0) + c
        for e, cc in sp_car(m).items():
            r2[e] = r2.get(e, 0) - c * cc
        r2 = {e: v for e, v in r2.items() if v}
    neg = malo or any(v < 0 for v in coef.values())
    cnt[("con coeficientes NEGATIVOS" if neg else "suma NO negativa", cond(list(lam)))] += 1
    if len(ej) < 7:
        ej.append((lam, {k: v for k, v in sorted(coef.items(), reverse=True)}, neg))

print("  ejemplos (lambda -> multiplicidades en la base de caracteres de Sp(2)):")
for (l, c, neg) in ej:
    print("     %-16s %-38s %s" % (str(l), str(c), "*** NEGATIVOS ***" if neg else ""))
print("")
print("  resumen (tipo, criterio (i) y (ii)):")
for k in sorted(cnt, key=str):
    print("     %-30s criterio=%-6s : %d" % (k[0], str(k[1]), cnt[k]))
neg = sum(v for k, v in cnt.items() if "NEGATIV" in k[0])
ceros = sum(v for k, v in cnt.items() if k[0] == "CERO")
mal = sum(v for k, v in cnt.items() if k[0] == "CERO" and k[1] is False)
mal += sum(v for k, v in cnt.items() if k[0] != "CERO" and k[1] is True)
print("")
print("  restricciones con algun coeficiente NEGATIVO : %d" % neg)
print("  -> %s" % ("LA PREDICCION CORREGIDA AGUANTA: combinacion entera NO NEGATIVA de caracteres de Sp(2)"
                   if neg == 0 else "tambien falla"))
print("  cero <-> (i) y (ii)  (fuera de la rama (a)) : %d desacuerdos sobre %d ceros" % (mal, ceros))
