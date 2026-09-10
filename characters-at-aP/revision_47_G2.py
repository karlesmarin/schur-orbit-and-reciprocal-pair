# -*- coding: utf-8 -*-
# VERIFICACION DE LA PARTE G2 DE LA RESPUESTA A LA CONSULTA 47.   25 de agosto de 2026.
#
#   R5  los seis valores (lambda+rho, alpha) en G2 son  u, v, u+v, 2u+v, 3u+v, 3u+2v
#       con u = a+1 y v = 3(b+1),  para lambda = a om_1 + b om_2.
#   R6  el criterio (+):  chi != 0  <=>  N_q(lambda) = r_q,  contra el evaluador de Freudenthal.
#   R7  las dos simplificaciones:  q=3 -> chi=0 <=> a = 2 mod 3;   q=2 -> chi=0 <=> a=b=1 mod 2.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python revision_47_G2.py > revision_47_G2_OUT.txt 2>&1

import cmath
import sys
from collections import Counter
from fractions import Fraction
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


GRAM = [[Fraction(2), Fraction(-3)], [Fraction(-3), Fraction(6)]]
ip = lambda u, v: sum(u[i] * v[j] * GRAM[i][j] for i in range(2) for j in range(2))
POS = [(1, 0), (0, 1), (1, 1), (2, 1), (3, 1), (3, 2)]
RHO = (5, 3)
OM = [(2, 1), (3, 2)]


def freudenthal(lam):
    lr = (lam[0] + RHO[0], lam[1] + RHO[1])
    cl = ip(lr, lr)
    N1, N2 = 2 * lam[0], 2 * lam[1]
    mult = {}
    for (n1, n2) in sorted(((x, y) for x in range(N1 + 1) for y in range(N2 + 1)),
                           key=lambda n: (n[0] + n[1], n)):
        mu = (lam[0] - n1, lam[1] - n2)
        if (n1, n2) == (0, 0):
            mult[mu] = Fraction(1)
            continue
        s = Fraction(0)
        for al in POS:
            k = 1
            while True:
                nu = (mu[0] + k * al[0], mu[1] + k * al[1])
                if nu not in mult:
                    break
                if mult[nu]:
                    s += mult[nu] * ip(nu, al)
                k += 1
        mr = (mu[0] + RHO[0], mu[1] + RHO[1])
        den = cl - ip(mr, mr)
        mult[mu] = Fraction(0) if den == 0 else 2 * s / den
    return {m_: int(v) for m_, v in mult.items() if v > 0}


def chi(lam, k, W):
    z = cmath.exp(2j * cmath.pi / 12)
    s = sum(m * z ** ((k * (mu[0] + 3 * mu[1])) % 12) for mu, m in W.items())
    return None if (abs(s.imag) > 1e-6 or abs(s.real - round(s.real)) > 1e-6) else int(round(s.real))


print(SEP)
print("R5 -- LOS SEIS VALORES (lambda+rho, alpha)  =  u, v, u+v, 2u+v, 3u+v, 3u+2v")
print("     con u = a+1,  v = 3(b+1)")
mal = 0
for a in range(6):
    for b in range(6):
        lam = (a * OM[0][0] + b * OM[1][0], a * OM[0][1] + b * OM[1][1])
        lr = (lam[0] + RHO[0], lam[1] + RHO[1])
        obs = sorted(int(ip(lr, al)) for al in POS)
        u, v = a + 1, 3 * (b + 1)
        pred = sorted([u, v, u + v, 2 * u + v, 3 * u + v, 3 * u + 2 * v])
        mal += (obs != pred)
        if (a, b) in ((0, 0), (1, 0), (0, 1), (2, 3)):
            print(f"   a={a} b={b}:  u={u} v={v}   observado {obs}   predicho {pred}")
ok(mal == 0, "R5: los seis valores son los de la formula", f"{mal} fallos en 36 pesos")

print(SEP)
print("R6 -- EL CRITERIO (+) EN G2, contra Freudenthal")
print(f"{'k':>3} {'d':>3} {'q':>4} {'r_q':>4} {'pesos':>6} {'(+) acierta':>12} {'nulos':>6}")
DOM = [(a, b) for a in range(7) for b in range(7)]
cache = {}
for (a, b) in DOM:
    cache[(a, b)] = freudenthal((a * OM[0][0] + b * OM[1][0], a * OM[0][1] + b * OM[1][1]))


def Nq(a, b, q):
    u, v = a + 1, 3 * (b + 1)
    return sum(1 for L in (u, v, u + v, 2 * u + v, 3 * u + v, 3 * u + 2 * v) if L % q == 0)


for k in range(1, 7):
    d = gcd(k, 12)
    q = 12 // d
    if q < 2:
        continue
    rq = Nq(0, 0, q)
    acc = tot = nul = 0
    for (a, b) in DOM:
        lam = (a * OM[0][0] + b * OM[1][0], a * OM[0][1] + b * OM[1][1])
        v = chi(lam, k, cache[(a, b)])
        if v is None:
            continue
        tot += 1
        nul += (v == 0)
        acc += ((Nq(a, b, q) == rq) == (v != 0))
    marca = "" if acc == tot else "   <-- FALLA"
    print(f"{k:>3} {d:>3} {q:>4} {rq:>4} {tot:>6} {str(acc)+'/'+str(tot):>12} {nul:>6}{marca}")
    if acc != tot:
        fallos.append(f"R6 k={k}")
ok(not [f for f in fallos if f.startswith("R6")], "R6: el criterio (+) cierra en G2")

print(SEP)
print("R7 -- LAS DOS SIMPLIFICACIONES")
for (k, q, nombre, regla) in ((4, 3, "q=3: chi=0 <=> a = 2 mod 3", lambda a, b: a % 3 == 2),
                              (6, 2, "q=2: chi=0 <=> a=b=1 mod 2", lambda a, b: a % 2 == 1 and b % 2 == 1)):
    acc = tot = 0
    for (a, b) in DOM:
        lam = (a * OM[0][0] + b * OM[1][0], a * OM[0][1] + b * OM[1][1])
        v = chi(lam, k, cache[(a, b)])
        if v is None:
            continue
        tot += 1
        acc += (regla(a, b) == (v == 0))
    ok(acc == tot, nombre, f"{acc}/{tot}")

print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {sorted(set(fallos))}")
    sys.exit(1)
print("RESULTADO: R5-R7 verificados.  El criterio de G2 esta cerrado.")
