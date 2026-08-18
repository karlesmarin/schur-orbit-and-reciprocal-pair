# -*- coding: utf-8 -*-
# Authors: Carles Marin, Claude (AI assistant).
"""El t-cociente SOLO, ¿determina la terna (d1,d2,d3)?

El resumen dice "the arguments read off the t-quotient".  La Prop. 3.4 los da con los TAMANOS del
cociente y los DOS RESIDUOS.  Si dos lambda en P_N comparten el t-cociente (como tupla indexada) y
difieren en la terna, la frase del resumen esta incompleta.

Convenio del paper: beta_j = lambda_j + N - j, N = t+2.  Clase i = {beta_j == i mod t}.
Componente del cociente para una clase con valores b_1 > ... > b_n, b_k = t*a_k + i:
    lambda^(i) = (a_1 - (n-1), a_2 - (n-2), ..., a_n)      (la prueba de Prop. 3.4 con n=2 da
                                                            (a_1 - 1, a_2), que es este caso)
"""
import itertools


def data(lam, t):
    N = t + 2
    L = list(lam) + [0] * (N - len(lam))
    beta = [L[j] + N - 1 - j for j in range(N)]
    if len(set(beta)) != N:
        return None
    cls = {}
    for b in beta:
        cls.setdefault(b % t, []).append(b)
    if len(cls) < t:
        return None, beta, None          # degenerada
    quot = {}
    for i in sorted(cls):
        vals = sorted(cls[i], reverse=True)
        n = len(vals)
        a = [(v - i) // t for v in vals]
        comp = tuple(a[k] - (n - 1 - k) for k in range(n))
        quot[i] = tuple(x for x in comp if x > 0) or ()
    exc = sorted(i for i in cls if len(cls[i]) >= 2)
    if len(exc) == 2:
        rA, rB = exc
        A = sorted(cls[rA], reverse=True)
        B = sorted(cls[rB], reverse=True)
    else:
        rA = rB = exc[0]
        p, q, r = sorted(cls[exc[0]], reverse=True)
        A, B = [p, q], [q, r]
    d1 = A[0] - A[1]
    d2 = B[0] - B[1]
    d3 = abs(A[0] + A[1] - B[0] - B[1])
    return tuple(quot[i] for i in sorted(quot)), beta, tuple(sorted((d1, d2, d3)))


def parts(maxsize, maxlen):
    out = []
    def rec(rem, mx, cur):
        if len(cur) <= maxlen:
            out.append(tuple(cur))
        if not rem or len(cur) == maxlen:
            return
        for k in range(min(rem, mx), 0, -1):
            rec(rem - k, k, cur + [k])
    for s in range(0, maxsize + 1):
        rec(s, s, [])
    return set(out)


for t in (3, 4, 5, 6):
    N = t + 2
    seen = {}
    clash = []
    for lam in parts(22, N):
        d = data(lam, t)
        if d is None or d[0] is None or d[2] is None:
            continue
        quot, beta, tri = d
        if quot in seen and seen[quot][1] != tri:
            clash.append((seen[quot], (lam, tri, beta)))
        elif quot not in seen:
            seen[quot] = (lam, tri, beta)
    print("t=%d : %d cocientes distintos ; choques (mismo cociente, terna distinta): %d"
          % (t, len(seen), len(clash)))
    for a, b in clash[:2]:
        print("    lambda=%s  terna=%s  beta=%s" % (a[0], a[1], a[2]))
        print("    lambda=%s  terna=%s  beta=%s" % (b[0], b[1], b[2]))
        print("    -> mismo t-cociente, ternas DISTINTAS")
