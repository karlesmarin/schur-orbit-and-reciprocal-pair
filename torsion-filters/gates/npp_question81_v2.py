# -*- coding: utf-8 -*-
# QUESTION 8.1, CON EL INSTRUMENTO YA CALIBRADO --- Y EL EJE QUE EL G2 DE NPP NO PUEDE VER.
# 20 de agosto de 2026.
#
# LO QUE EL CALIBRADO DEJO CLARO (npp_calibrado_G2.py):
#   - la dicotomia de NPP, medida dentro del G2 de su §7, es una EQUIVALENCIA: 169/169, las dos
#     mitades, con dim Z_G(C2) = 6 saliendo de rho y el denominador 8 de su Prop 7.1 reproducido;
#   - PERO el senuelo --- emparejar con la RAIZ en vez de con la CORRAIZ --- no discrimina alli, y
#     no por casualidad: ese ejemplo vive en d = 2, y raiz y corraiz de G2 difieren por un factor 3,
#     que es IMPAR.  A d = 2 los dos emparejamientos tienen la misma paridad SIEMPRE.
#     El G2 calibra todo menos ese eje.
#
#   lectura A (corraices de C_m = raices del dual B_m):  <l+r, e_i>  y  <l+r, e_i +- e_j>
#   lectura B (raices de C_m, el propio grupo):          <l+r, 2e_i> y  <l+r, e_i +- e_j>
#
# QUE SE MIDE
#   Q1  las dos lecturas contra sp_lambda(x_m) exacto, en rango ANCHO, casilla por casilla.
#   Q2  el punto donde difieren, exhibido: la lambda mas pequena y su cuenta.
#   Q3  control de que las dos lecturas NO son la misma funcion en este rango --- si coincidieran,
#       la medida no decidiria nada, que es lo que paso en G2.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python npp_question81_v2.py > npp_question81_v2_OUT.txt 2>&1

import itertools
import json
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

MS = [1, 2, 3]
LMAX = 9


def ciclotomico(t):
    poly = [-1] + [0] * (t - 1) + [1]
    for d in range(1, t):
        if t % d == 0:
            poly = _div(poly, ciclotomico(d))
    return poly


def _div(a, b):
    a = a[:]
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] // b[-1]
        q[i] = c
        for j in range(len(b)):
            a[i + j] -= c * b[j]
    assert all(x == 0 for x in a)
    return q


def _red(a, phi):
    a = a[:]
    n = len(phi) - 1
    for i in range(len(a) - 1, n - 1, -1):
        c = a[i]
        if c:
            a[i] = 0
            for j in range(n):
                a[i - n + j] -= c * phi[j]
    while len(a) > n:
        a.pop()
    while len(a) < n:
        a.append(0)
    return a


def _mul(a, b, phi):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return _red(out, phi)


def xip(k, t, phi):
    k %= t
    n = len(phi) - 1
    if k < n:
        v = [0] * n
        v[k] = 1
        return v
    e = [0] * n
    e[1] = 1
    r = [1] + [0] * (n - 1)
    for _ in range(k):
        r = _mul(r, e, phi)
    return r


def det(M, phi):
    n = len(M)
    if n == 1:
        return M[0][0][:]
    acc = [0] * (len(phi) - 1)
    for j in range(n):
        men = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        term = _mul(M[0][j], det(men, phi), phi)
        if j % 2:
            term = [-x for x in term]
        acc = [u + v for u, v in zip(acc, term)]
    return acc


def sp_valor(lam, m, a, t, phi):
    rho = [m - i for i in range(m)]
    lr = [lam[i] + rho[i] for i in range(m)]

    def bialt(ex):
        M = [[[u - v for u, v in zip(xip(ex[i] * a[j], t, phi), xip(-ex[i] * a[j], t, phi))]
              for j in range(m)] for i in range(m)]
        return det(M, phi)

    num, den = bialt(lr), bialt(rho)
    if all(x == 0 for x in den):
        return None
    for c in range(-6, 7):
        if [c * x for x in den] == num:
            return c
    return "otro"


def hits(lr, m, d, lectura):
    """#{raices con d | <lambda+rho, .>}, contadas en los dos signos."""
    n = 0
    for i in range(m):
        v = lr[i] if lectura == "A" else 2 * lr[i]
        if v % d == 0:
            n += 2
    for i, j in itertools.combinations(range(m), 2):
        for s in (1, -1):
            if (lr[i] + s * lr[j]) % d == 0:
                n += 2
    return n


print("=" * 100)
print("QUESTION 8.1 CON EL INSTRUMENTO CALIBRADO --- LAS DOS LECTURAS")
print("=" * 100)

resumen = {}
for m in MS:
    t = 2 * m + 2
    phi = ciclotomico(t)
    a = list(range(1, m + 1))
    d = None
    for k in range(1, t + 1):
        ex = [(k * ai) % t for ai in a]
        if all(e == 0 for e in ex) or (t % 2 == 0 and all(e == t // 2 for e in ex)):
            d = k
            break

    lams = [l for l in itertools.product(range(0, LMAX + 1), repeat=m)
            if all(l[i] >= l[i + 1] for i in range(m - 1))]
    cA, cB = Counter(), Counter()
    difieren = 0
    primera = None
    malA, malB = [], []
    for lam in lams:
        v = sp_valor(list(lam), m, a, t, phi)
        if v is None or isinstance(v, str):
            continue
        lr = [lam[i] + (m - i) for i in range(m)]
        hA, hB = hits(lr, m, d, "A"), hits(lr, m, d, "B")
        pA, pB = hA > 0, hB > 0
        cA[(pA, v == 0)] += 1
        cB[(pB, v == 0)] += 1
        if pA != pB:
            difieren += 1
            if primera is None:
                primera = (lam, v, hA, hB)
        if pA != (v == 0) and len(malA) < 4:
            malA.append((lam, v, hA))
        if pB != (v == 0) and len(malB) < 4:
            malB.append((lam, v, hB))

    n = sum(cA.values())
    acA = cA[(True, True)] + cA[(False, False)]
    acB = cB[(True, True)] + cB[(False, False)]
    print("")
    print("-" * 100)
    print("m=%d  Sp(%d)  t=%d  d=%d   %d particiones (lambda_1 <= %d)" % (m, 2 * m, t, d, n, LMAX))
    print("-" * 100)
    print("  lectura A  CORRAICES de C_m  (= raices del dual B_m) : %d de %d   %s"
          % (acA, n, "EQUIVALENCIA" if acA == n else "falla en %d, p.ej. %s" % (n - acA, malA[:3])))
    print("  lectura B  RAICES de C_m     (el propio grupo)       : %d de %d   %s"
          % (acB, n, "EQUIVALENCIA" if acB == n else "falla en %d, p.ej. %s" % (n - acB, malB[:3])))
    print("  Q3  lambda donde las dos lecturas DIFIEREN : %d   %s"
          % (difieren, "*** ninguna: el rango no decide ***" if difieren == 0
             else "el eje ES visible aqui"))
    if primera:
        print("      la menor: lambda=%s  sp=%s  hits(corraiz)=%d  hits(raiz)=%d"
              % (primera[0], primera[1], primera[2], primera[3]))
    resumen["m=%d" % m] = dict(t=t, d=d, n=int(n), A=int(acA), B=int(acB), difieren=int(difieren),
                               primera=[list(primera[0]), primera[1], primera[2], primera[3]]
                               if primera else None)

json.dump(resumen, open("npp_question81_v2_DUMP.json", "w"), indent=1)
print("")
print("=" * 100)
print("DONE")
