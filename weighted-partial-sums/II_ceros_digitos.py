# -*- coding: utf-8 -*-
r"""II_ceros_digitos.py -- la recursion por digitos de S(m) = sum n chi(n) y sus reglas de ceros.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

chi = chi_{q,s}: Legendre mod q (q primo impar) en los n primos con q, chi(q) = s = +-1.
C(r) = sum_{r'<=r} chi(r'), T(r) = sum_{r'<=r} r' chi(r'), T_q = T(q-1)  (0 <= r < q).
(D)  S(qk + r) = s q S(k) + k (T_q + q C(r)) + T(r)          [a comprobar en todo qk+r <= MMAX]
Reglas que se siguen de (D):
 (R1) si S(k) = 0 y r esta en R = {r : T_q + q C(r) = 0, T(r) = 0}, entonces S(qk + r) = 0.
      Par (q = 1 mod 4): T_q = 0, R = {C = T = 0} contiene 0 y q-1.
      Impar (q = 3 mod 4, q > 3): T_q = -q h(-q), R = {C(r) = h, T(r) = 0}.
 (R2) pendientes, r = 0:  S(qk)/(qk) = s S(k)/k + T_q/q   (impar: s S(k)/k - h).
      Con S(k) = 0 y h = 1: S(q(k+1)) = q (k+1)(s chi(k+1) - 1) = 0 si chi(k+1) = s.
Se mide: (D) exacta; |R| por q; ceros 'primitivos' = ceros no alcanzados desde ceros menores por R1
(ni por R2 con h=1); crecimiento del numero de ceros frente a X^{log|R|/log q}.
Uso: python II_ceros_digitos.py [MMAX]
"""
import sys, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 7


def legendre(q):
    t = np.zeros(q, dtype=np.int64)
    for r in range(1, q):
        t[r * r % q] = 1
    t = np.where(t == 1, 1, -1); t[0] = 0
    return t


def S_de(q, s, N):
    L = legendre(q)
    x = np.arange(N + 1, dtype=np.int64)
    eps = np.ones(N + 1, dtype=np.int64)
    while True:
        div = (x % q == 0) & (x > 0)
        if not div.any():
            break
        x[div] //= q; eps[div] *= s
    eps = eps * L[x % q]; eps[0] = 0
    return np.cumsum(np.arange(N + 1, dtype=np.int64) * eps)


for q in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 163):
    L = legendre(q)
    C = np.cumsum(L); T = np.cumsum(np.arange(q) * L); Tq = int(T[q - 1])
    h = -Tq // q if q % 4 == 3 else 0
    R = [r for r in range(q) if Tq + q * int(C[r]) == 0 and int(T[r]) == 0]
    for s in (1, -1):
        S = S_de(q, s, MMAX)
        m = np.arange(MMAX + 1)
        k, r = m // q, m % q
        D = s * q * S[k] + k * (Tq + q * C[r]) + T[r]
        okD = bool(np.all(D[q:] == S[q:]))
        Z = np.nonzero(S[1:] == 0)[0] + 1
        Zs = set(int(z) for z in Z)
        # R1: cierre; primitivos
        hijos = set()
        for z in Zs:
            for rr in R:
                c = q * z + rr
                if c <= MMAX:
                    assert c in Zs, ("R1 falla", q, s, z, rr)
                    hijos.add(c)
            if h == 1 and (z + 1) <= MMAX // q and int(L[(z + 1) % q] if (z + 1) % q else 0) != 0:
                pass
        # R2 con h = 1: S(k)=0, chi(k+1)=s  =>  S(q(k+1)) = 0
        r2 = set()
        if h == 1:
            for z in Zs:
                c = q * (z + 1)
                if c <= MMAX:
                    # chi(z+1) calculado de S
                    chi = (int(S[z + 1]) - int(S[z])) // (z + 1)
                    if chi == s:
                        assert c in Zs, ("R2 falla", q, s, z)
                        r2.add(c)
        prim = sorted(Zs - hijos - r2)
        exp = math.log(len(R)) / math.log(q) if len(R) > 1 else 0
        print("q=%3d q%%4=%d h=%d s=%+d  (D) %s  R=%s  ceros=%6d  primitivos=%4d %s  X^(log|R|/log q)=%.1f"
              % (q, q % 4, h, s, "OK" if okD else "FALLA", R, len(Zs), len(prim), prim[:10],
                 MMAX ** exp if exp else 0))
        sys.stdout.flush()
