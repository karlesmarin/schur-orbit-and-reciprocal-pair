# -*- coding: utf-8 -*-
# Authors: Carles Marin, Claude (AI assistant).
"""Verificacion INDEPENDIENTE del Teorema 3.1, escrita solo desde su enunciado.

Escrita el 2026-08-15 durante el repaso con lente de arbitro: NO usa ningun guion ni ninguna salida
del repositorio, solo el enunciado impreso del teorema y la definicion de bialternante.

Lado izquierdo : el bialternante s_lambda(1, zeta, ..., zeta^{t-1}, z, 1/z), a 60 digitos.
Lado derecho   : eps_lambda * sinh(d1 th/2) sinh(d2 th/2) sinh(d3 th/2) / (sinh^2(t th/2) sinh th),
                 con d1,d2,d3 y eps leidos del beta set como dice el teorema, y eps por eq:sign.

No se usa nada del repositorio: ni sus guiones, ni sus salidas.  Controles al final.
"""
import itertools
import mpmath as mp

mp.mp.dps = 60


def inv_count(w):
    return sum(1 for a in range(len(w)) for b in range(a + 1, len(w)) if w[a] > w[b])


def bialternant(beta, alphabet):
    N = len(beta)
    num = mp.matrix(N, N)
    den = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            num[i, j] = alphabet[i] ** beta[j]
            den[i, j] = alphabet[i] ** (N - 1 - j)
    return mp.det(num) / mp.det(den)


def closed_form(beta, t, z):
    """Devuelve (valor, etiqueta) segun el Teorema 3.1 leido literalmente."""
    N = len(beta)
    cls = {}
    for j, b in enumerate(beta):
        cls.setdefault(b % t, []).append((b, j))
    if len(cls) < t:
        return mp.mpf(0), "degenerado"
    exc = sorted(i for i in cls if len(cls[i]) >= 2)
    if len(exc) == 2:
        rA, rB = exc                                     # r_A <= r_B
        A = sorted(cls[rA], reverse=True)                # [(a1,jA1),(a2,.)]
        B = sorted(cls[rB], reverse=True)
    else:
        i0 = exc[0]
        p, q, r = sorted(cls[i0], reverse=True)          # p>q>r
        rA = rB = i0
        A, B = [p, q], [q, r]
    a1, jA1 = A[0]
    a2, _ = A[1]
    b1, jB1 = B[0]
    b2, _ = B[1]
    d1 = a1 - a2
    d2 = b1 - b2
    dt3 = a1 + a2 - b1 - b2
    d3 = abs(dt3)
    if dt3 == 0:
        return mp.mpf(0), "concentrico"
    S = [j for j in range(N) if j not in (jA1, jB1)]
    bS = [beta[j] % t for j in S]
    sgn = lambda x: (1 if x > 0 else (-1 if x < 0 else 0))
    eps = ((-1) ** (t + (N + 1) * N // 2)
           * (-1) ** ((jA1 + 1) + (jB1 + 1) + inv_count(bS))
           * sgn(a1 - b1) * sgn(dt3))
    th = mp.log(z)
    val = (eps * mp.sinh(d1 * th / 2) * mp.sinh(d2 * th / 2) * mp.sinh(d3 * th / 2)
           / (mp.sinh(t * th / 2) ** 2 * mp.sinh(th)))
    return val, "formula"


def partitions(maxsize, maxlen):
    out = []
    def rec(rem, mx, cur):
        out.append(tuple(cur))
        if not rem or len(cur) == maxlen:
            return
        for k in range(min(rem, mx), 0, -1):
            rec(rem - k, k, cur + [k])
    for s in range(maxsize + 1):
        rec(s, s, [])
    return sorted(set(out))


ZS = [mp.mpf(3) / 2, mp.mpf(7) / 5, mp.mpf(11) / 4]
print("  t  |lam|<=   formas   ceros   comprobadas   FALLOS")
grand_fail = 0
for t, MS in ((2, 12), (3, 12), (4, 11), (5, 10), (6, 9)):
    N = t + 2
    nf = nz = nc = bad = 0
    for lam in partitions(MS, N):
        L = list(lam) + [0] * (N - len(lam))
        beta = [L[j] + N - 1 - j for j in range(N)]
        nf += 1
        ok = True
        for z in ZS:
            alphabet = [mp.exp(2j * mp.pi * k / t) for k in range(t)] + [mp.mpf(z), 1 / mp.mpf(z)]
            lhs = bialternant(beta, beta and beta or beta) if False else bialternant(beta, alphabet)
            rhs, tag = closed_form(beta, t, z)
            if abs(lhs - rhs) > mp.mpf(10) ** (-25) * max(1, abs(lhs)):
                ok = False
        if not ok:
            bad += 1
            if bad <= 2:
                print("      FALLO lambda=%s beta=%s" % (lam, beta))
        _, tag = closed_form(beta, t, ZS[0])
        if tag != "formula":
            nz += 1
        nc += 1
    grand_fail += bad
    print("  %2d %6d %8d %7d %13d %8d" % (t, MS, nf, nz, nc, bad))
print()
print("TOTAL FALLOS: %d" % grand_fail)

# --- control: la formula TIENE que fallar si le cambio el signo a eps en una forma no nula
print()
print("CONTROL (tiene que fallar): mismo test con eps forzado a +1 siempre")
bad2 = 0
tested = 0
for t, MS in ((3, 9), (4, 9)):
    N = t + 2
    for lam in partitions(MS, N):
        L = list(lam) + [0] * (N - len(lam))
        beta = [L[j] + N - 1 - j for j in range(N)]
        rhs, tag = closed_form(beta, t, ZS[0])
        if tag != "formula":
            continue
        tested += 1
        z = ZS[0]
        alphabet = [mp.exp(2j * mp.pi * k / t) for k in range(t)] + [mp.mpf(z), 1 / mp.mpf(z)]
        lhs = bialternant(beta, alphabet)
        forced = abs(rhs)                       # eps := +1
        if abs(lhs - forced) > mp.mpf(10) ** (-25) * max(1, abs(lhs)):
            bad2 += 1
print("  el control falla en %d de %d formas no nulas (si fuera 0, el test no mide el signo)" % (bad2, tested))
