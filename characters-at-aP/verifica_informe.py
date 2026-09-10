# -*- coding: utf-8 -*-
"""Las cuatro afirmaciones matematicas del informe, recomputadas desde cero.
Caracteres por Jacobi-Trudi simplectico con el h del alfabeto del propio articulo."""
import sys
from math import gcd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def h_closed(q, d, N):
    """coeficientes de H(z)=(1-z)^2/(1-z^q)^d  (d par)  o  (1-z^2)/(1-z^q)^d  (d impar)."""
    base = [0] * (N + 1)
    k = 0
    while q * k <= N:
        c = 1
        for i in range(1, d):
            c = c * (k + i) // i
        base[q * k] = c
        k += 1
    if d % 2 == 0:
        return [base[n] - (2 * base[n - 1] if n >= 1 else 0) + (base[n - 2] if n >= 2 else 0)
                for n in range(N + 1)]
    return [base[n] - (base[n - 2] if n >= 2 else 0) for n in range(N + 1)]


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


def chi(lam, m, k):
    t = 2 * m + 2
    d = gcd(k, t)
    q = t // d
    h = h_closed(q, d, 3 * (max(lam) + m) + 8)
    H = lambda n: (h[n] if 0 <= n < len(h) else 0)
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    return det_int(M) // 2


def raices_C(m):
    R = []
    for i in range(m):
        for j in range(i + 1, m):
            R.append(tuple(1 if x == i else (-1 if x == j else 0) for x in range(m)))
            R.append(tuple(1 if x in (i, j) else 0 for x in range(m)))
    for i in range(m):
        R.append(tuple(2 if x == i else 0 for x in range(m)))
    return R


def datos(vec, m, q):
    R = raices_C(m)
    vals = [sum(vec[i] * a[i] for i in range(m)) for a in R]
    N = sum(1 for v in vals if v % q == 0)
    Z = 1
    for v in vals:
        if v % q == 0:
            Z *= v
    ocup = {}
    for x in vec:
        c = x % q
        c = min(c, q - c)
        ocup[c] = ocup.get(c, 0) + 1
    return N, Z, tuple(sorted(ocup.items())), sorted(v % q for v in vals)


print("=" * 88)
print("(a) INFORME 1.3:  Sp(6), q=4.  rho=(3,2,1) y l=(5,4,3) -- mismo dato de raices?")
m, k = 3, 2
t = 2 * m + 2; d = gcd(k, t); q = t // d
rho = (3, 2, 1); ell = (5, 4, 3)
Nr, Zr, or_, resr = datos(rho, m, q)
Ne, Ze, oe, rese = datos(ell, m, q)
print("    t=%d k=%d -> d=%d q=%d" % (t, k, d, q))
print("    residuos por raiz iguales: %s" % (resr == rese))
print("    N_q(rho)=%d  N_q(l)=%d   (r_q=%d, los dos sobreviven: %s)"
      % (Nr, Ne, Nr, Nr == Ne))
print("    ocupaciones plegadas   rho %s   l %s   -> DISTINTAS: %s" % (or_, oe, or_ != oe))
print("    Z_q(l)/Z_q(rho) = %d/%d = %s" % (Ze, Zr, Ze // Zr if Zr else "-"))
print("    caracteres exactos:  chi_0 = %d    chi_(2,2,2) = %d"
      % (chi((0, 0, 0), m, k), chi((2, 2, 2), m, k)))

print()
print("=" * 88)
print("(b) INFORME 1.3 bis:  el testigo minimo del signo en Sp(4), k=2")
m, k = 2, 2
for lam in [(1, 0), (1, 1), (10, 10), (10, 9)]:
    ell = tuple(lam[i] + (m - i) for i in range(m))
    t = 2 * m + 2; d = gcd(k, t); q = t // d
    N, Z, oc, _ = datos(ell, m, q)
    print("    lambda=%-9s l=%-9s ocupaciones %-14s chi = %+d"
          % (str(lam), str(ell), str(oc), chi(lam, m, k)))
print("    -> (1,0) y (1,1) tienen las MISMAS ocupaciones y signos opuestos, con |lambda| 1 y 2.")
print("       El articulo llama minimo al par (10,10)/(10,9), de |lambda| 20 y 19.")

print()
print("=" * 88)
print("(c) INFORME 1.2:  la definicion literal del resumen aplicada a lambda=0")
m, k = 3, 2
t = 2 * m + 2; d = gcd(k, t); q = t // d
R = raices_C(m)
cero = (0,) * m
rho = tuple(m - i for i in range(m))
N_lit = sum(1 for a in R if sum(cero[i] * a[i] for i in range(m)) % q == 0)
r_q = sum(1 for a in R if sum(rho[i] * a[i] for i in range(m)) % q == 0)
print("    Sp(6), q=%d:  N_q(0) leido literalmente = %d  (= todas las %d raices positivas)"
      % (q, N_lit, len(R)))
print("    r_q = N_q(rho) = %d.   El resumen pide N_q(lambda)=r_q, o sea %d=%d: FALSO"
      % (r_q, N_lit, r_q))
print("    pero chi_0(a_P^2) = %d, que no es cero.  La convencion del Teorema si funciona:"
      % chi(cero, m, k))
print("    N_q(lambda+rho)=%d = r_q=%d." % (r_q, r_q))

print()
print("=" * 88)
print("(d) INFORME 1.1:  S(q)=S(q-1) en el regimen regular, por rango")


def S(m, q):
    f = 1 if q % 2 else 2
    a = (q - f) // 2
    u, v = divmod(m, q)
    def fact(n):
        r = 1
        for i in range(2, n + 1):
            r *= i
        return r
    def binom(n, r):
        if r < 0 or r > n:
            return 0
        return fact(n) // (fact(r) * fact(n - r))
    if v <= a:
        num = fact(m) * 2 ** (2 * u * a + v) * binom(a, v)
        den = fact(2 * u + 1) ** v * fact(2 * u) ** (a - v) * fact(u) ** f
    else:
        num = fact(m) * 2 ** (a * (2 * u + 1)) * binom(a + f, v - a)
        den = fact(2 * u + 1) ** a * fact(u) ** f * (u + 1) ** (v - a)
    return num // den if den and num % den == 0 else num / den


for m in range(2, 7):
    fila = []
    for q in range(2 * m + 2, 2 * m + 12, 2):
        fila.append("q=%d: %s" % (q, "SI" if S(m, q) == S(m, q - 1) else "NO"))
    print("    C_%d (regimen regular, q>2m=%d):  %s" % (m, 2 * m, "  ".join(fila)))
print("    -> se cumple en TODO rango pasado el rango singular.")
for m in (2, 3):
    print("    C_%d dentro del rango singular: S(%d)=%s  S(%d)=%s"
          % (m, 2 * m, S(m, 2 * m), 2 * m - 1, S(m, 2 * m - 1)))
