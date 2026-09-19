# -*- coding: utf-8 -*-
r"""II_bloques_ceros.py -- ceros infinitos sin R_q: el lema de concatenacion de bloques, y R_q en q = 1 mod 4.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

chi = chi_{q,+1}, q = 3 mod 4, h = h(-q), D(r) = C(r) - h.  De (D):  S(qk + r) = q S(k) + q k D(r) + T(r).
Lema (bloques): si B tiene digitos d_{L-1..0} en base q con sum D(d_i) = 0 y S(B) = h B, entonces
S(k) = h k  =>  S(q^L k + B) = h (q^L k + B); y S(q m) = q (S(m) - h m).  Asi, con k_t = B (1 + Q + ...
+ Q^{t-1}), Q = q^L, los m_t = q k_t son ceros.
Se comprueba DIRECTAMENTE (suma explicita, sin (D)):
 (1) q = 59, h = 3, B = (1,22,21,0)_59 = 283200: sum D = 0, S(B) = 3B, y S(59 B) = S(16708800) = 0;
     el cero 4591439 = (22,21,0,0)_59.
 (2) q = 347, h = 5: el bloque de 13 digitos (40,49,226,119,164,210,117,28,302,58,264,0,0): sum D = 0 y
     S(B) = 5B, por la recursion (D) (ya comprobada exacta hasta 10^7 y formalizada en Lean).
 (3) R_103 = {47, 51, 55}  (en q = 7 mod 8, |R_q| puede ser mayor que 1).
 (4) q = 5 mod 24, q <= QMAX:  (q-5)/6 en R_q;  q = 1 mod 4: |R_q| par y >= 2.
Uso: python II_bloques_ceros.py [QMAX]
"""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000


def leg(q):
    t = np.full(q, -1, dtype=np.int64); t[0] = 0
    r = np.arange(1, q, dtype=np.int64)
    t[(r * r) % q] = 1
    return t


def tablas(q):
    L = leg(q); r = np.arange(q, dtype=np.int64)
    return L, np.cumsum(L), np.cumsum(r * L)


def S_directo(q, s, N):
    L = leg(q)
    x = np.arange(N + 1, dtype=np.int64)
    eps = np.ones(N + 1, dtype=np.int64)
    while True:
        div = (x % q == 0) & (x > 0)
        if not div.any():
            break
        x[div] //= q; eps[div] *= s
    eps = eps * L[x % q]; eps[0] = 0
    return np.cumsum(np.arange(N + 1, dtype=np.int64) * eps)


def S_rec(q, s, m, C, T):
    """S por la recursion (D), entero exacto de Python (para m enormes)."""
    if m == 0:
        return 0
    k, r = divmod(m, q)
    Tq = int(T[q - 1])
    return s * q * S_rec(q, s, k, C, T) + k * (Tq + q * int(C[r])) + int(T[r])


def digitos(m, q):
    d = []
    while m:
        d.append(m % q); m //= q
    return d[::-1]


# (1) q = 59
q, h = 59, 3
L, C, T = tablas(q)
B = 1 * 59 ** 3 + 22 * 59 ** 2 + 21 * 59 + 0
sD = sum(int(C[d]) - h for d in (1, 22, 21, 0))
S = S_directo(q, 1, 59 * B)
print("(1) q=59: B=%d, sum D=%d, S(B)=%d (3B=%d), S(59B)=S(%d)=%d, S(4591439)=%d, digitos(4591439)=%s"
      % (B, sD, int(S[B]), 3 * B, 59 * B, int(S[59 * B]), int(S[4591439]), digitos(4591439, 59)))
Qb = 59 ** 4
for t in (2, 3):
    kt = B * sum(Qb ** i for i in range(t))
    print("    t=%d: S(k_t)-3k_t = %d ; S(59 k_t) = %d  (por la recursion)"
          % (t, S_rec(q, 1, kt, C, T) - 3 * kt, S_rec(q, 1, 59 * kt, C, T)))

# (2) q = 347
q, h = 347, 5
L, C, T = tablas(q)
dig = (40, 49, 226, 119, 164, 210, 117, 28, 302, 58, 264, 0, 0)
B = 0
for d in dig:
    B = B * q + d
print("(2) q=347: B=%d, sum D=%d, S(B)-5B=%d, S(347B)=%d"
      % (B, sum(int(C[d]) - h for d in dig), S_rec(q, 1, B, C, T) - 5 * B, S_rec(q, 1, q * B, C, T)))

# (3) R_103
q = 103
L, C, T = tablas(q); Tq = int(T[-1])
R = [r for r in range(q) if Tq + q * int(C[r]) == 0 and int(T[r]) == 0]
print("(3) R_103 =", R)

# (4)
es = np.ones(QMAX + 1, dtype=bool); es[:2] = False
for p in range(2, int(QMAX ** .5) + 1):
    if es[p]:
        es[p * p::p] = False
m24 = ok24 = m14 = ok14 = 0
for q in np.nonzero(es)[0]:
    q = int(q)
    if q < 5 or q % 4 != 1:
        continue
    L, C, T = tablas(q)
    enR = (C == 0) & (T == 0)
    k = int(enR.sum())
    m14 += 1; ok14 += (k % 2 == 0 and k >= 2)
    if q % 24 == 5:
        m24 += 1; ok24 += bool(enR[(q - 5) // 6])
print("(4) q=5 mod 24 <= %d: (q-5)/6 en R_q en %d de %d ;  q=1 mod 4: |R_q| par y >=2 en %d de %d"
      % (QMAX, ok24, m24, ok14, m14))
