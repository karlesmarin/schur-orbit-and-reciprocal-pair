# -*- coding: utf-8 -*-
# VERIFICACION DEL SEGUNDO INFORME.   8 de septiembre de 2026.
#
# Cinco afirmaciones nuevas, todas medibles.  Se miden ANTES de escribir nada, y con senuelos.
#
#  V1  IDENTIDAD DE RESIDUOS.  Con P_x(X) = sum_i (X^{x_i} + X^{-x_i})  y
#      R_x(X) = sum_{alpha in Phi(C_m)} X^{(x,alpha)}  en Z[X]/(X^q - 1),
#            R_x = ( P_x(X)^2 + P_x(X^2) ) / 2  -  m.
#      Si vale, y si ademas R_ell = R_rho en los supervivientes, (B) queda DEMOSTRADO y deja de ser
#      una medida.
#  V2  p SE CAE:  chi_lambda(a_P^k) = chi_lambda(a_P^d)  con d = gcd(k, 2m+2).
#  V3  EL SIGNO POR RESIDUOS ORDENADOS:  sgn chi = (-1)^{F_q(ell mod q) - F_q(rho mod q)}  con
#      F_q(r) = #{i : 2r_i >= q} + #{i<j : r_i < r_j} + #{i<j : r_i + r_j >= q}.
#  V4  RECUENTO PARA TODO q, las dos ramas de la formula del informe.
#  V5  LOS DOS CONTRAEJEMPLOS:  (a) m=6,k=7: el grupo de la factorizacion es Sp6 x SO7 mientras el
#      centralizador es Sp6 x Sp6;  (b) Sp4, lambda=(1,0), q=7: Z_7(ell)/Z_7(rho) = 1 pero el
#      caracter vale 2cos(2pi/7)+2cos(4pi/7), luego la formula del valor NO vale en toda la recta.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python referee2.py > referee2_OUT.txt 2>&1

import sys
from fractions import Fraction
from itertools import product as iproduct
from math import comb, cos, factorial, gcd, pi

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 96


def dominants(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for rest in rec(k - 1, v):
                yield (v,) + rest
    return rec(m, top)


def pares_C_todas(v):
    """TODAS las raices de C_m, positivas y negativas."""
    m = len(v)
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            for s in (v[i] - v[j], v[j] - v[i], v[i] + v[j], -v[i] - v[j]):
                out.append(s)
    for i in range(m):
        out.append(2 * v[i]); out.append(-2 * v[i])
    return out


def pares_C(v):
    m = len(v)
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            out.append(v[i] - v[j]); out.append(v[i] + v[j])
    for i in range(m):
        out.append(2 * v[i])
    return out


def h_closed(q, d, N):
    base = [0] * (N + 1)
    k = 0
    while q * k <= N:
        c = 1
        for i in range(1, d):
            c = c * (k + i) // i
        base[q * k] = c
        k += 1
    return [base[n] - (2 * base[n - 1] if n >= 1 else 0) + (base[n - 2] if n >= 2 else 0)
            if d % 2 == 0 else base[n] - (base[n - 2] if n >= 2 else 0) for n in range(N + 1)]


def det_int(M):
    n = len(M)
    A = [r[:] for r in M]
    sg, prev = 1, 1
    for c in range(n):
        if A[c][c] == 0:
            p = next((r for r in range(c + 1, n) if A[r][c] != 0), None)
            if p is None:
                return 0
            A[c], A[p] = A[p], A[c]; sg = -sg
        for r in range(c + 1, n):
            for k in range(c + 1, n):
                A[r][k] = (A[r][k] * A[c][c] - A[r][c] * A[c][k]) // prev
            A[r][c] = 0
        prev = A[c][c]
    return sg * A[n - 1][n - 1]


def chi(lam, m, h):
    H = lambda n: (h[n] if 0 <= n < len(h) else 0)
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    return det_int(M) // 2


def vec_R(x, q):
    """el multiconjunto de residuos de (x,alpha) sobre TODAS las raices, como vector mod q."""
    v = [0] * q
    for a in pares_C_todas(x):
        v[a % q] += 1
    return v


def vec_R_identidad(x, q):
    """el mismo, por la identidad del informe: (P^2 + P(X^2))/2 - m."""
    m = len(x)
    P = [0] * q
    for xi in x:
        P[xi % q] += 1
        P[(-xi) % q] += 1
    P2 = [0] * q
    for i in range(q):
        if P[i]:
            for j in range(q):
                if P[j]:
                    P2[(i + j) % q] += P[i] * P[j]
    Pdos = [0] * q
    for i in range(q):
        Pdos[(2 * i) % q] += P[i]
    out = [(P2[i] + Pdos[i]) for i in range(q)]
    if any(o % 2 for o in out):
        return None
    out = [o // 2 for o in out]
    out[0] -= m
    return out


# --------------------------------------------------------------------------- V1
print(SEP)
print("V1 -- la identidad  R_x = (P_x^2 + P_x(X^2))/2 - m,  y  R_ell = R_rho en supervivientes")
mal_id = mal_eq = tot_id = tot_sup = 0
for m in range(2, 8):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    top = 8 if m <= 4 else 6
    for k in range(1, t):
        d = gcd(k, t); q = t // d
        if q < 2:
            continue
        r_q = sum(1 for b in pares_C(rho) if b % q == 0)
        Rrho = vec_R(rho, q)
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            tot_id += 1
            if vec_R_identidad(ell, q) != vec_R(ell, q):
                mal_id += 1
            if sum(1 for a in pares_C(ell) if a % q == 0) != r_q:
                continue
            tot_sup += 1
            if vec_R(ell, q) != Rrho:
                mal_eq += 1
print(f"   la identidad: {tot_id} casos, {mal_id} fallos")
print(f"   R_ell = R_rho en supervivientes: {tot_sup} casos, {mal_eq} fallos")

# --------------------------------------------------------------------------- V2
print(SEP)
print("V2 -- ¿se cae el numerador p?   chi(a_P^k) = chi(a_P^d)")
mal_p = tot_p = 0
for m in range(2, 8):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    top = 8 if m <= 4 else 6
    for k in range(1, t):
        d = gcd(k, t); q = t // d
        if q < 2 or k == d:
            continue
        hk = h_closed(q, d, 3 * (top + m) + 6)
        for lam in dominants(m, top):
            tot_p += 1
            if chi(lam, m, hk) != chi(lam, m, hk):   # mismo h: d y k dan el MISMO alfabeto
                mal_p += 1
# comparacion de verdad: el alfabeto de a_P^k depende solo de (d,q), luego h es el mismo objeto.
print("   el alfabeto de a_P^k solo depende de (d,q) por la Prop. del alfabeto: h_closed(q,d)")
print("   asi que chi(a_P^k) y chi(a_P^d) se calculan con la MISMA serie -> son iguales por")
print("   construccion.  Se comprueba contra el calculo directo desde los autovalores:")
import cmath
malv = totv = 0
for m in (2, 3, 4, 5):
    t = 2 * m + 2
    for k in range(1, t):
        d = gcd(k, t); q = t // d
        if q < 2:
            continue
        A = sorted(round(cmath.phase(cmath.exp(2j * cmath.pi * ((j * k) % t) / t)), 9)
                   for j in range(1, m + 1))
        B = sorted(round(cmath.phase(cmath.exp(2j * cmath.pi * ((j * d) % t) / t)), 9)
                   for j in range(1, m + 1))
        Aset = sorted(round(x, 6) for x in
                      [cmath.exp(2j * cmath.pi * ((j * k) % t) / t).real for j in range(1, m + 1)])
        Bset = sorted(round(x, 6) for x in
                      [cmath.exp(2j * cmath.pi * ((j * d) % t) / t).real for j in range(1, m + 1)])
        totv += 1
        if Aset != Bset:
            malv += 1
print(f"   espectros de a_P^k y a_P^d iguales como multiconjunto: {totv} pares, {malv} distintos")

# --------------------------------------------------------------------------- V3
print(SEP)
print("V3 -- el signo por residuos ordenados:  F_q")


def F_q(r, q):
    m = len(r)
    s = sum(1 for x in r if 2 * x >= q)
    for i in range(m):
        for j in range(i + 1, m):
            if r[i] < r[j]:
                s += 1
            if r[i] + r[j] >= q:
                s += 1
    return s


mal_s = tot_s = 0
fall = []
for m in range(2, 8):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    top = 8 if m <= 4 else 6
    for k in range(1, t):
        d = gcd(k, t); q = t // d
        if q < 2:
            continue
        h = h_closed(q, d, 3 * (top + m) + 6)
        r_q = sum(1 for b in pares_C(rho) if b % q == 0)
        Frho = F_q([x % q for x in rho], q)
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            if sum(1 for a in pares_C(ell) if a % q == 0) != r_q:
                continue
            v = chi(lam, m, h)
            if v == 0:
                continue
            tot_s += 1
            pred = 1 if (F_q([x % q for x in ell], q) - Frho) % 2 == 0 else -1
            if pred != (1 if v > 0 else -1):
                mal_s += 1
                if len(fall) < 4:
                    fall.append((m, k, q, lam, v))
print(f"   {tot_s} supervivientes, {mal_s} fallos")
for f in fall:
    print("    fallo:", f)

# --------------------------------------------------------------------------- V4
print(SEP)
print("V4 -- recuento para TODO q, las dos ramas del informe")


def S_formula(m, q):
    f = 1 if q % 2 else 2
    a = (q - f) // 2
    u, v = divmod(m, q)
    if a == 0:
        return None
    if v <= a:
        num = Fraction(factorial(m) * 2 ** (2 * u * a + v) * comb(a, v))
        den = (factorial(2 * u + 1) ** v) * (factorial(2 * u) ** (a - v)) * (factorial(u) ** f)
        return num / den
    num = Fraction(factorial(m) * 2 ** (a * (2 * u + 1)) * comb(a + f, v - a))
    den = (factorial(2 * u + 1) ** a) * (factorial(u) ** f) * ((u + 1) ** (v - a))
    return num / den


print(f"   {'m':>3} {'q':>4} {'S medido':>10} {'formula':>12} {'=?':>4}")
mal4 = tot4 = 0
for m in range(1, 5):
    for q in range(2, 13):
        if q ** m * (m * m) > 4_000_000:
            continue
        rho = [m - i for i in range(m)]
        r_q = sum(1 for b in pares_C(rho) if b % q == 0)
        S = 0
        for x in iproduct(range(q), repeat=m):
            ell = [sum(x[j] for j in range(i, m)) for i in range(m)]
            if sum(1 for a in pares_C(ell) if a % q == 0) == r_q:
                S += 1
        F = S_formula(m, q)
        if F is None:
            continue
        tot4 += 1
        okk = (F.denominator == 1 and int(F) == S)
        if not okk:
            mal4 += 1
        print(f"   {m:>3} {q:>4} {S:>10} {str(F):>12} {'SI' if okk else 'no':>4}")
print(f"   {tot4} pares medidos, {mal4} discrepancias")

# --------------------------------------------------------------------------- V5
print(SEP)
print("V5 -- los dos contraejemplos del informe")
# (a) m=6, k=7
m, k = 6, 7
t = 2 * m + 2
d = gcd(k, t); q = t // d
rho = [m - i for i in range(m)]
print(f"   (a) m={m}, k={k}: d={d}, q={q}")
ocup = {}
for x in rho:
    c = min(x % q, (-x) % q)
    ocup[c] = ocup.get(c, 0) + 1
print(f"       ocupaciones de rho por clase plegada: {ocup}")
esp = [((j * k) % t) for j in range(1, m + 1)]
print(f"       exponentes del espectro de a_P^k mod t: {sorted(esp)}  (t={t})")
# (b) Sp4, lambda=(1,0), q=7
m = 2
rho = [2, 1]
lam = (1, 0)
ell = [lam[i] + rho[i] for i in range(m)]
Z = [a for a in pares_C(ell) if a % 7 == 0]
Zr = [a for a in pares_C(rho) if a % 7 == 0]
val = 2 * cos(2 * pi / 7) + 2 * cos(4 * pi / 7)
print(f"   (b) Sp4, lambda={lam}, q=7: ell={ell}; (ell,alpha)={sorted(pares_C(ell))}")
print(f"       divisibles por 7: numerador {Z}, denominador {Zr} -> cociente = 1")
print(f"       caracter estandar en g_(1/7) = {val:.10f}   -> la formula del valor NO vale ahi")
print(SEP)
