# Copyright (c) 2026 Carles Marin. Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Nota del conductor, §8: (a) Psi_q = Psi_{q'}^E mod p (capas en el anillo real);
# (b) rangos por orbitas gcd(c,15) para q'=15, p=7 (el soporte de caracteres es una union, no una suma).
# Uso: python psi_y_orbitas.py
# Comprobaciones: (a) Psi_q = Psi_{q'}^E mod p ; (b) rangos por orbitas q'=15, p=7
from math import gcd
import sympy as sp

X = sp.symbols("X")


def psi(n):
    """polinomio minimo de zeta_n + zeta_n^{-1} (n >= 3), via resultante."""
    y = sp.symbols("y")
    phi = sp.cyclotomic_poly(n, y)
    r = sp.resultant(phi, y**2 - X * y + 1, y)  # = psi(X)^2 salvo constante
    f = sp.factor_list(sp.Poly(r, X))[1]
    cands = [g for g, e in f]
    assert len(cands) == 1, cands
    return sp.Poly(cands[0], X)


def phi_(n):
    return sp.totient(n)


print("(a) Psi_q == Psi_{q'}^E (mod p)")
for (q, p) in [(35, 5), (45, 3), (56, 2), (63, 3), (175, 5), (99, 3), (40, 2), (147, 7)]:
    v = 0
    qq = q
    while qq % p == 0:
        qq //= p
        v += 1
    E = phi_(p**v)
    a = psi(q)
    b = psi(qq) ** E
    diff = (a - b).set_modulus(p) if False else sp.Poly([c % p for c in (a - b).all_coeffs()], X)
    ok = all(c % p == 0 for c in (a - b).all_coeffs())
    print("   q=%d p=%d q'=%d E=%d  deg %d = %d  congruentes: %s" % (q, p, qq, E, a.degree(), b.degree(), ok))


def rango_mod_p(filas, p):
    M = [list(f) for f in filas]
    r = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], -1, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


print("(b) q'=15, p=7, S(c)=2c-15 (c!=0): rangos por orbita gcd(c,15) y total")
n, p = 15, 7
units = [u for u in range(1, n) if gcd(u, n) == 1]
S = lambda c: 0 if c % n == 0 else (2 * (c % n) - n)
cols_all = list(range(1, n))
for g in (1, 3, 5):
    cols = [c for c in cols_all if gcd(c, n) == g]
    print("   gcd=%d  rango %d" % (g, rango_mod_p([[S(pow(u, -1, n) * c) for c in cols] for u in units], p)))
print("   total  rango %d" % rango_mod_p([[S(pow(u, -1, n) * c) for c in cols_all] for u in units], p))
