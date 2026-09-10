# -*- coding: utf-8 -*-
# VERIFICACION DEL INFORME (2): la factorizacion en dimensiones y el recuento cerrado.
# 8 de septiembre de 2026.
#
# (1) FACTORIZACION.  Para lambda superviviente, el informe afirma
#        |chi| = kappa(lambda) * D_0 * D_{q/2} * prod_{clases no fijas} D_c
#     con  D_c = dim V^{GL_n}_{nu(c)}  en una clase no fija de n filas, construida asi: se orientan
#     las filas con signos para que eps_i ell_i == c (mod q), se ordenan y_i = (eps_i ell_i - c)/q
#     de mayor a menor, y nu_i = y_i - y_n - (n-i);
#        D_0 = dim V^{Sp_2n}_{nu(0)}  con  nu_i = ell_i/q - (n-i+1);
#        D_{q/2} = dim V^{SO_{2n+1}}_{nu(q/2)}  con  nu_i = ell_i/q - (n-i+1/2);
#     y kappa = 1 si d impar, 1 si d par y la plaza vacia esta en la clase 0, 2 en otro caso.
#
# (2) RECUENTO cerrado para q | 2m+2:
#        S_{C_m}(q) = m! 2^{ad} / ((d!)^a (s!)^f) * (1 si d impar; s(a+f) si d par),
#     con f = 1 (q impar) o 2 (q par), a = (q-f)/2, s = floor(d/2).
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python dims_recuento.py > dims_recuento_OUT.txt 2>&1

import sys
from fractions import Fraction
from itertools import product as iproduct
from math import factorial, gcd

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
            if d % 2 == 0 else base[n] - (base[n - 2] if n >= 2 else 0)
            for n in range(N + 1)]


def det_int(M):
    n = len(M)
    A = [row[:] for row in M]
    sign, prev = 1, 1
    for c in range(n):
        if A[c][c] == 0:
            piv = next((r for r in range(c + 1, n) if A[r][c] != 0), None)
            if piv is None:
                return 0
            A[c], A[piv] = A[piv], A[c]
            sign = -sign
        for r in range(c + 1, n):
            for kk in range(c + 1, n):
                A[r][kk] = (A[r][kk] * A[c][c] - A[r][c] * A[c][kk]) // prev
            A[r][c] = 0
        prev = A[c][c]
    return sign * A[n - 1][n - 1]


def chi(lam, m, h):
    H = lambda n: (h[n] if 0 <= n < len(h) else 0)
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    return det_int(M) // 2


def pares_C(v):
    m = len(v)
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            out.append(v[i] - v[j]); out.append(v[i] + v[j])
    for i in range(m):
        out.append(2 * v[i])
    return out


def cls(x, q):
    return min(x % q, (-x) % q)


def dim_gl(nu):
    """dim de la irrep de GL_n de peso maximo nu (particion)."""
    n = len(nu)
    l = [nu[i] + (n - i - 1) for i in range(n)]
    num = 1
    for i in range(n):
        for j in range(i + 1, n):
            num *= (l[i] - l[j])
    den = 1
    for j in range(1, n):
        den *= factorial(j)
    return num // den


def dim_sp(nu):
    """dim de la irrep de Sp(2n) de peso maximo nu."""
    n = len(nu)
    l = [Fraction(nu[i] + (n - i)) for i in range(n)]
    r = [Fraction(n - i) for i in range(n)]
    num = den = Fraction(1)
    for i in range(n):
        num *= l[i]; den *= r[i]
        for j in range(i + 1, n):
            num *= (l[i] - l[j]) * (l[i] + l[j])
            den *= (r[i] - r[j]) * (r[i] + r[j])
    v = num / den
    return int(v)


def dim_so_odd(nu):
    """dim de la irrep de SO(2n+1) de peso maximo nu (puede ser semientero)."""
    n = len(nu)
    l = [Fraction(nu[i]) + Fraction(2 * (n - i) - 1, 2) for i in range(n)]
    r = [Fraction(2 * (n - i) - 1, 2) for i in range(n)]
    num = den = Fraction(1)
    for i in range(n):
        num *= l[i]; den *= r[i]
        for j in range(i + 1, n):
            num *= (l[i] - l[j]) * (l[i] + l[j])
            den *= (r[i] - r[j]) * (r[i] + r[j])
    v = num / den
    return int(v)


print(SEP)
print("(1) FACTORIZACION EN DIMENSIONES DE GRUPOS CLASICOS")
print(SEP)
print(f"{'m':>3} {'k':>3} {'d':>3} {'q':>4} {'vivos':>7} {'factoriz.':>11} {'kappa=1 mal':>12}")
T = {"ok": 0, "bad": 0}
fallos = []
for m in range(2, 8):
    t = 2 * m + 2
    rho = [m - i for i in range(m)]
    Prho = pares_C(rho)
    top = 8 if m <= 4 else 6
    for k in range(1, t):
        d = gcd(k, t); q = t // d
        if q < 2:
            continue
        h = h_closed(q, d, 3 * (top + m) + 6)
        r_q = sum(1 for b in Prho if b % q == 0)
        fijas = {0} | ({q // 2} if q % 2 == 0 else set())
        vivos = ok = bad = kmal = 0
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            if sum(1 for a in pares_C(ell) if a % q == 0) != r_q:
                continue
            vivos += 1
            v = abs(chi(lam, m, h))
            porclase = {}
            for e in ell:
                porclase.setdefault(cls(e, q), []).append(e)
            D = 1
            for c, filas in sorted(porclase.items()):
                n = len(filas)
                if c == 0:
                    nu = sorted([f // q - (n - i + 1) for i, f in
                                 enumerate(sorted(filas, reverse=True), start=1)], reverse=True)
                    D *= dim_sp(nu)
                elif q % 2 == 0 and c == q // 2:
                    ff = sorted(filas, reverse=True)
                    nu = [Fraction(ff[i], q) - Fraction(2 * (n - i) - 1, 2) for i in range(n)]
                    D *= dim_so_odd(nu)
                else:
                    ys = sorted([(e - c) // q if (e - c) % q == 0 else -((e + c) // q)
                                 for e in filas], reverse=True)
                    nu = [ys[i] - ys[n - 1] - (n - 1 - i) for i in range(n)]
                    D *= dim_gl(nu)
            # kappa
            faltan = [c for c in range(q // 2 + 1)
                      if len(porclase.get(c, [])) < ((d // 2) if c in fijas else d)]
            if d % 2:
                kap = 1
            else:
                kap = 1 if (faltan and faltan[0] == 0) else 2
            ok, bad = (ok + 1, bad) if kap * D == v else (ok, bad + 1)
            if kap * D != v and 1 * D == v:
                kmal += 1
            if kap * D != v and len(fallos) < 6:
                fallos.append((m, k, d, q, lam, v, kap * D))
        if vivos == 0:
            continue
        T["ok"] += ok; T["bad"] += bad
        print(f"{m:>3} {k:>3} {d:>3} {q:>4} {vivos:>7} {ok:>5}/{vivos:<5} {kmal:>12}"
              + ("" if bad == 0 else "   <-- FALLA"))
print(f"   TOTAL factorizacion: {T['ok']} aciertos, {T['bad']} fallos")
if fallos:
    for f in fallos:
        print("    fallo:", f)

print(SEP)
print("(2) RECUENTO CERRADO  S_{C_m}(q)  para q | 2m+2")
print(f"   {'m':>3} {'q':>4} {'d':>3} {'S medido':>10} {'formula':>10} {'=?':>4}")
mal2 = 0
for m in range(1, 6):
    t = 2 * m + 2
    for q in range(2, t + 1):
        if t % q:
            continue
        d = t // q
        f = 1 if q % 2 else 2
        a = (q - f) // 2
        s = d // 2
        F = Fraction(factorial(m) * (2 ** (a * d)),
                     (factorial(d) ** a) * (factorial(s) ** f))
        F = F * (1 if d % 2 else s * (a + f))
        F = int(F) if F.denominator == 1 else F
        S = 0
        rho = [m - i for i in range(m)]
        r_q = sum(1 for b in pares_C(rho) if b % q == 0)
        for x in iproduct(range(q), repeat=m):
            ell = [sum(x[j] for j in range(i, m)) for i in range(m)]
            if sum(1 for aa in pares_C(ell) if aa % q == 0) == r_q:
                S += 1
        marca = "SI" if S == F else "no"
        if S != F:
            mal2 += 1
        print(f"   {m:>3} {q:>4} {d:>3} {S:>10} {F:>10} {marca:>4}")
print(SEP)
print(f"VEREDICTO: factorizacion {'OK' if T['bad'] == 0 else 'FALLA'};"
      f"  recuento cerrado {'OK' if mal2 == 0 else '%d discrepancias' % mal2}")
