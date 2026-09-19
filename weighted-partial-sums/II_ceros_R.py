# -*- coding: utf-8 -*-
r"""II_ceros_R.py -- el conjunto de digitos R_q que propaga ceros, para todo primo q <= QMAX.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Con (D) (II_ceros_digitos.py):  S(qk + r) = s q S(k) + k (T_q + q C(r)) + T(r).  S(0) = 0, asi que
todo m > 0 cuyos digitos en base q estan en R_q = {r : T_q + q C(r) = 0, T(r) = 0} es un cero de S,
para los dos signos s = chi(q).
Se comprueba para cada primo impar q <= QMAX:
 (1) semisuma ponderada (q = 3 mod 4, q > 3):  A_q = sum_{r<q/2} r chi(r) = q (1 - chi(2)) h / 2,
     con h = -T_q/q (Dirichlet) y la semisuma sin peso B_q = (2 - chi(2)) h.
 (2) q = 7 mod 8  =>  (q-1)/2 esta en R_q;  q = 3 mod 8  =>  R_q vacio?;  q = 1 mod 4 => 0, q-1 en R_q.
 (3) R_q completo y su tamano; simetria r <-> q-1-r en el caso par.
 (4) control: con los digitos de R_q, S calculado DIRECTAMENTE (suma explicita) en los tres primeros
     numeros de digitos en R_q de cada longitud <= 3 que no pasen de 2*10^6, es 0.
Uso: python II_ceros_R.py [QMAX]
"""
import sys
from collections import Counter
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000


def es_primo(n):
    return n > 1 and all(n % d for d in range(2, int(n ** .5) + 1))


def leg(q):
    t = [-1] * q; t[0] = 0
    for r in range(1, q):
        t[r * r % q] = 1
    return t


def chi_fun(q, s, L):
    def chi(n):
        v = 1
        while n % q == 0:
            n //= q; v *= s
        return v * L[n % q]
    return chi


resumen = Counter(); fallos = []
tam = Counter()
for q in range(3, QMAX + 1):
    if not es_primo(q):
        continue
    L = leg(q)
    C = []; T = []; c = t = 0
    for r in range(q):
        c += L[r]; t += r * L[r]; C.append(c); T.append(t)
    Tq = T[q - 1]
    R = [r for r in range(q) if Tq + q * C[r] == 0 and T[r] == 0]
    if q % 4 == 3 and q > 3:
        h = -Tq // q
        c2 = L[2]
        A = sum(r * L[r] for r in range(1, (q + 1) // 2))
        B = sum(L[r] for r in range(1, (q + 1) // 2))
        if A * 2 != q * (1 - c2) * h or B != (2 - c2) * h:
            fallos.append(("semisuma", q))
        if q % 8 == 7:
            resumen["7mod8, (q-1)/2 en R"] += ((q - 1) // 2 in R)
            resumen["7mod8 total"] += 1
            tam["7mod8 |R|=%d" % len(R)] += 1
        else:
            resumen["3mod8, R vacio"] += (len(R) == 0)
            resumen["3mod8 total"] += 1
            if R:
                fallos.append(("3mod8 R no vacio", q, R))
    elif q % 4 == 1:
        resumen["1mod4, 0 y q-1 en R"] += (0 in R and q - 1 in R)
        resumen["1mod4 total"] += 1
        resumen["1mod4 R simetrico"] += all((q - 1 - r) in R for r in R)
        tam["1mod4 |R|=%d" % len(R)] += 1
    # control directo
    for s in (1, -1):
        chi = chi_fun(q, s, L)
        ds = [r for r in R if r > 0] or []
        cands = []
        for d1 in R:
            for d0 in R:
                m = d1 * q + d0
                if 0 < m <= 2 * 10 ** 6:
                    cands.append(m)
        for d2 in R[:2]:
            for d1 in R[:2]:
                for d0 in R[:2]:
                    m = (d2 * q + d1) * q + d0
                    if 0 < m <= 2 * 10 ** 6:
                        cands.append(m)
        for m in sorted(set(cands))[:6]:
            if sum(n * chi(n) for n in range(1, m + 1)) != 0:
                fallos.append(("control", q, s, m))
print("primos impares <= %d" % QMAX)
for k in sorted(resumen):
    print("  %-28s %d" % (k, resumen[k]))
print("tamanos de R:", dict(sorted(tam.items())))
print("fallos:", fallos[:10], "(total %d)" % len(fallos))
