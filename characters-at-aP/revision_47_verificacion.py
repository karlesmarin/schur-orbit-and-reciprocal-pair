# -*- coding: utf-8 -*-
# VERIFICACION DE LA RESPUESTA A LA CONSULTA 47.   25 de agosto de 2026.
#
# POR QUE.  La respuesta afirma cinco cosas, y dos de ellas contradicen lo que la nota dice.  Nada
# se toca hasta comprobarlas.  Se miden POR SEPARADO, y cada una con su senuelo cuando cabe.
#
#   R1  EL LEMA DE FILAS.  D_{P(S)h}(b) = (prod_i P(T_i)) D_h(b), con T_i el desplazamiento de b_i.
#       Si es cierto, la "obstruccion" que la nota declara en su §7 NO EXISTE: el desplazamiento S
#       sobre el indice de h es el mismo desplazamiento de b_i en los DOS sumandos h_{b_i+j} y
#       h_{b_i-j}, porque desplaza el indice, no j.
#
#   R2  LA ESPECIALIZACION PRINCIPAL (*):
#          chi_lambda(a_P(z)) = z^{-(lambda,rho)} prod_{alpha>0} (1-z^{(l,alpha)})/(1-z^{(rho,alpha)})
#       con l = lambda+rho y el estadistico de RAICES (no de corraices).
#
#   R3  EL CRITERIO (+):  chi != 0  <=>  N_q(lambda) = r_q,  con
#          N_q(lambda) = #{alpha>0 : q | (l,alpha)},   r_q = N_q(0) = |R_P(q)^+|.
#       Se contrasta contra los 21630 pesos ya medidos.
#
#   R4  EL CONTEO (++):  N_q(lambda) = sum_{c no fija} C(n_c,2) + sum_{f fija} n_f^2.
#
#   R5  G2: los seis valores (l,alpha) = u, v, u+v, 2u+v, 3u+v, 3u+2v con u=a+1, v=3(b+1);
#       y el criterio (+) contra el evaluador de Freudenthal.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python revision_47_verificacion.py > revision_47_verificacion_OUT.txt 2>&1

import cmath
import itertools
import random
import sys
from collections import Counter
from fractions import Fraction
from math import comb, gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 94
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


def det_f(M):
    n = len(M)
    M = [r[:] for r in M]
    d = 1.0
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-14:
            return 0.0
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return d


# ------------------------------------------------------------------ R1
print(SEP)
print("R1 -- EL LEMA DE FILAS:  D_{P(S)h}(b) = (prod_i P(T_i)) D_h(b)")
print("     h aleatoria; b aleatorio; P(S) = (1-S)^2 y (1-S^2).  Si esto vale, la obstruccion que")
print("     la nota declara en su §7 no existe.")
random.seed(3)
NH = 80


def D(h, b):
    m = len(b)
    H = lambda n: (h[n] if 0 <= n < len(h) else 0.0)
    return det_f([[H(b[i] + j) + H(b[i] - j) for j in range(1, m + 1)] for i in range(m)])


def aplica(h, coef):
    """(sum_r coef[r] S^r) h,  con (S h)_n = h_{n-1}."""
    return [sum(coef[r] * (h[n - r] if 0 <= n - r < len(h) else 0.0) for r in range(len(coef)))
            for n in range(len(h))]


for (nombre, coef) in (("(1-S)^2", [1.0, -2.0, 1.0]), ("(1-S^2)", [1.0, 0.0, -1.0])):
    peor = 0.0
    for _ in range(60):
        m = random.choice([2, 3, 4])
        h = [random.uniform(-1, 1) for _ in range(NH)]
        b = sorted(random.sample(range(6, 30), m), reverse=True)
        izq = D(aplica(h, coef), b)
        # derecha: prod_i P(T_i) aplicado a D, o sea suma sobre eps in {0..len(coef)-1}^m
        der = 0.0
        for eps in itertools.product(range(len(coef)), repeat=m):
            c = 1.0
            for e in eps:
                c *= coef[e]
            if c:
                der += c * D(h, [b[i] - eps[i] for i in range(m)])
        peor = max(peor, abs(izq - der) / max(1.0, abs(izq)))
    ok(peor < 1e-9, f"R1 con P = {nombre}", f"peor error relativo {peor:.2e}")

# senuelo: si el lema fuera falso, una variante con signos cambiados no deberia funcionar
peor_sen = 0.0
for _ in range(30):
    m = random.choice([2, 3])
    h = [random.uniform(-1, 1) for _ in range(NH)]
    b = sorted(random.sample(range(6, 30), m), reverse=True)
    izq = D(aplica(h, [1.0, -2.0, 1.0]), b)
    der = 0.0
    for eps in itertools.product(range(3), repeat=m):
        c = 1.0
        for e in eps:
            c *= [1.0, 2.0, 1.0][e]              # <- senuelo: +2 en vez de -2
        der += c * D(h, [b[i] - eps[i] for i in range(m)])
    peor_sen = max(peor_sen, abs(izq - der) / max(1.0, abs(izq)))
ok(peor_sen > 1e-3, "SENUELO: con los coeficientes cambiados NO funciona",
   f"error del senuelo {peor_sen:.2e}")

# ------------------------------------------------------------------ R2
print(SEP)
print("R2 -- LA ESPECIALIZACION PRINCIPAL EN TIPO C, en la recta a_P(z) = (z, z^2, ..., z^m)")


def sp_bialt_z(lam, z, m):
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    xs = [z ** (i + 1) for i in range(m)]
    num = det_f([[x ** e - x ** (-e) for x in xs] for e in ell])
    den = det_f([[x ** e - x ** (-e) for x in xs] for e in rho])
    return num / den


def producto(lam, z, m):
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    pares_l, pares_r = [], []
    for i in range(m):
        pares_l.append(2 * ell[i]); pares_r.append(2 * rho[i])
        for j in range(i + 1, m):
            pares_l += [ell[i] - ell[j], ell[i] + ell[j]]
            pares_r += [rho[i] - rho[j], rho[i] + rho[j]]
    p = 1.0
    for a, b in zip(pares_l, pares_r):
        p *= (1 - z ** a) / (1 - z ** b)
    lr = sum(lam[i] * rho[i] for i in range(m))
    return z ** (-lr) * p


print(f"{'m':>3} {'lambda':>14} {'z':>6} {'bialternante':>14} {'producto (*)':>14} {'err rel':>10}")
peor = 0.0
for m in (2, 3, 4):
    for lam in ([3, 1] if m == 2 else [3, 2, 0] if m == 3 else [4, 2, 1, 0],
                [5, 2] if m == 2 else [5, 1, 1] if m == 3 else [6, 3, 2, 1]):
        if len(lam) != m:
            continue
        for z in (1.37, 0.82):
            a = sp_bialt_z(lam, z, m)
            b = producto(lam, z, m)
            e = abs(a - b) / max(1.0, abs(a))
            peor = max(peor, e)
            print(f"{m:>3} {str(lam):>14} {z:>6} {a:>14.6f} {b:>14.6f} {e:>10.1e}")
ok(peor < 1e-8, "R2: la especializacion principal (*) se cumple", f"peor {peor:.2e}")

# ------------------------------------------------------------------ R3 / R4
print(SEP)
print("R3/R4 -- EL CRITERIO (+) Y EL CONTEO (++), contra el evaluador exacto")


def h_aP(m, k):
    t = 2 * m + 2
    d = gcd(k, t)
    q = t // d
    def c(n):
        return 0 if (n < 0 or n % q) else comb(n // q + d - 1, d - 1)
    if d % 2 == 0:
        return lambda n: 0 if n < 0 else c(n) - 2 * c(n - 1) + c(n - 2)
    return lambda n: 0 if n < 0 else c(n) - c(n - 2)


def det_int(M):
    n = len(M)
    if n == 0:
        return 1
    M = [r[:] for r in M]
    sg, prev = 1, 1
    for cc in range(n - 1):
        if M[cc][cc] == 0:
            p = next((r for r in range(cc + 1, n) if M[r][cc] != 0), None)
            if p is None:
                return 0
            M[cc], M[p] = M[p], M[cc]
            sg = -sg
        for r in range(cc + 1, n):
            for k2 in range(cc + 1, n):
                M[r][k2] = (M[r][k2] * M[cc][cc] - M[r][cc] * M[cc][k2]) // prev
            M[r][cc] = 0
        prev = M[cc][cc]
    return sg * M[n - 1][n - 1]


def chi_exacto(lam, m, H):
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    v = det_int(M)
    assert v % 2 == 0 or v == 0
    return v // 2


def Nq(ell, q, m):
    n = 0
    for i in range(m):
        n += (2 * ell[i]) % q == 0
        for j in range(i + 1, m):
            n += (ell[i] - ell[j]) % q == 0
            n += (ell[i] + ell[j]) % q == 0
    return n


def conteo_clases(ell, q, m):
    cls = lambda x: min(x % q, (-x) % q)
    fijas = {0} | ({q // 2} if q % 2 == 0 else set())
    cnt = Counter(cls(e) for e in ell)
    s = 0
    for c, n in cnt.items():
        s += n * n if c in fijas else comb(n, 2)
    return s


print(f"{'m':>3} {'k':>3} {'d':>3} {'q':>4} {'r_q':>4} {'pesos':>6} "
      f"{'(+) acierta':>12} {'(++) acierta':>13}")
for m in range(2, 8):
    t = 2 * m + 2
    for k in range(1, 7):
        d = gcd(k, t)
        q = t // d
        if q < 3:
            continue
        H = h_aP(m, k)
        rho = [m - i for i in range(m)]
        rq = Nq(rho, q, m)
        top = 7 if m <= 4 else 5
        a3 = a4 = tot = 0
        for lam in itertools.product(range(top + 1), repeat=m):
            if any(lam[i] < lam[i + 1] for i in range(m - 1)):
                continue
            lam = list(lam)
            ell = [lam[i] + rho[i] for i in range(m)]
            v = chi_exacto(lam, m, H)
            N = Nq(ell, q, m)
            tot += 1
            a3 += ((N == rq) == (v != 0))
            a4 += (N == conteo_clases(ell, q, m))
        marca = "" if (a3 == tot and a4 == tot) else "   <-- FALLA"
        print(f"{m:>3} {k:>3} {d:>3} {q:>4} {rq:>4} {tot:>6} {str(a3)+'/'+str(tot):>12} "
              f"{str(a4)+'/'+str(tot):>13}{marca}")
        if a3 != tot:
            fallos.append(f"R3 m={m} k={k}")
        if a4 != tot:
            fallos.append(f"R4 m={m} k={k}")
ok(not [f for f in fallos if f.startswith("R3")], "R3: el criterio N_q(lambda) = r_q se cumple")
ok(not [f for f in fallos if f.startswith("R4")], "R4: el conteo por clases se cumple")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {sorted(set(fallos))[:8]}")
    sys.exit(1)
print("RESULTADO: R1-R4 verificados.")
