# -*- coding: utf-8 -*-
# LA QUESTION 8.1 DE NPP, PROBADA CONTRA NUESTRA FAMILIA.   20 de agosto de 2026.
#
# POR QUE.  Sus autores plantean la Question 8.1 y la dejan ABIERTA.  Tenemos los dos lados que
# hacen falta para medirla sobre nuestra familia: los elementos x_m y el valor del caracter en
# ellos.  Falta traducir el enunciado y medirlo.
#
# EL ENUNCIADO, LITERAL (arXiv:2504.14684, Question 8.1), leido del PDF:
#   «Let x0 be a torsion element of G, of order d in the adjoint group of G.  Treat lambda.rho as the
#    cocharacter lambda.rho : C^x -> That.  Then if the dimension of the centralizer of
#    (lambda.rho)(e^{2 pi i/d}) in Ghat is GREATER than the dimension of the centralizer of x0 in G,
#    then is it true that Theta_lambda(x0) = 0?  If ... EQUAL ..., then Theta_lambda(x0) is up to a
#    constant, the dimension of an irreducible representation of a reductive group whose dual group is
#    the connected component of identity of the centralizer of (lambda.rho)(e^{2 pi i/d}) in Ghat.»
#
# LA TRADUCCION, y es donde se puede fallar --- por eso va escrita antes de medir:
#   G = Sp(2m), simplemente conexo.  Ghat = SO(2m+1), tipo B_m, rango m.
#   x_m regular  =>  Z_G(x_m) = toro maximal  =>  dim Z_G(x_m) = m.
#   lambda.rho es multiplicativo en pesos, o sea lambda + rho aditivo: (lambda_i + m - i + 1).
#   y = (lambda.rho)(e^{2 pi i/d})  =>  alpha(y) = 1  <=>  d | <alpha, lambda+rho>.
#   dim Z_Ghat(y) = rango + #{alpha en Phi(B_m) : alpha(y) = 1}.
#   Luego su dicotomia, para NUESTRO x_m, es exactamente:
#
#       sp_lambda(x_m) = 0   <=>   algun alpha de B_m cumple  d | <alpha, lambda+rho>
#
#   y en el caso de igualdad el centralizador es un toro, cuya representacion irreducible tiene
#   dimension 1, luego el valor debe ser +-1 (por la constante).  Eso tambien se mide.
#
# QUE SE MIDE
#   C0  CALIBRACION, y va primero: d = orden de x_m en el grupo ADJUNTO Sp(2m)/{+-I}, calculado y no
#       supuesto.  Si d se equivoca, todo lo demas mide otra cosa.
#   C1  FATAL: la equivalencia de arriba, casilla por casilla, contra sp_lambda(x_m) exacto.
#   C2  en los no nulos, .es el valor +-1 como predice la mitad de igualdad?
#   D1  SENUELO: la misma cuenta con las raices de C_m (el propio grupo) en vez de las de B_m (el
#       dual).  Si NO discrimina, la medida no distingue grupo de dual y no vale.
#   D2  SENUELO: usar d = t en vez del orden en el adjunto.  Tiene que fallar donde difieran.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python npp_question81.py > npp_question81_OUT.txt 2>&1

import itertools
import json
import math
import sys
from collections import Counter

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

MS = [1, 2, 3]
LMAX = 6


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
    lam = list(lam) + [0] * (m - len(lam))
    rho = [m - i for i in range(m)]
    lr = [lam[i] + rho[i] for i in range(m)]

    def bialt(ex):
        M = [[[u - v for u, v in zip(xip(ex[i] * a[j], t, phi), xip(-ex[i] * a[j], t, phi))]
              for j in range(m)] for i in range(m)]
        return det(M, phi)

    num, den = bialt(lr), bialt(rho)
    if all(x == 0 for x in den):
        return None
    for c in range(-4, 5):
        if [c * x for x in den] == num:
            return c
    return "otro"


def raices_B(m):
    R = []
    for i in range(m):
        R += [tuple(1 if k == i else 0 for k in range(m)),
              tuple(-1 if k == i else 0 for k in range(m))]
    for i, j in itertools.combinations(range(m), 2):
        for si in (1, -1):
            for sj in (1, -1):
                R.append(tuple(si if k == i else (sj if k == j else 0) for k in range(m)))
    return R


def raices_C(m):
    R = []
    for i in range(m):
        R += [tuple(2 if k == i else 0 for k in range(m)),
              tuple(-2 if k == i else 0 for k in range(m))]
    for i, j in itertools.combinations(range(m), 2):
        for si in (1, -1):
            for sj in (1, -1):
                R.append(tuple(si if k == i else (sj if k == j else 0) for k in range(m)))
    return R


print("=" * 100)
print("LA QUESTION 8.1 DE NPP, MEDIDA CONTRA NUESTRA FAMILIA")
print("=" * 100)

resumen = {}
for m in MS:
    t = 2 * m + 2
    phi = ciclotomico(t)
    a = list(range(1, m + 1))

    # ---- C0  el orden en el grupo ADJUNTO Sp(2m)/{+-I}: menor k>0 con x^k = +-I
    d = None
    for k in range(1, t + 1):
        ex = [(k * ai) % t for ai in a]
        # x^k = +-I  <=>  todos los xi^{k a_i} iguales a 1, o todos iguales a -1
        if all(e == 0 for e in ex):
            d = k
            break
        if t % 2 == 0 and all(e == t // 2 for e in ex):
            d = k
            break
    print("")
    print("-" * 100)
    print("m=%d  Sp(%d)   t = %d   d (orden en el adjunto) = %s   rango = %d" % (m, 2 * m, t, d, m))
    print("-" * 100)

    lams = [l for l in itertools.product(range(0, LMAX + 1), repeat=m)
            if all(l[i] >= l[i + 1] for i in range(m - 1))]
    RB, RC = raices_B(m), raices_C(m)

    c1 = Counter()
    c2 = Counter()
    d1 = Counter()
    d2 = Counter()
    mal = []
    for lam in lams:
        v = sp_valor(lam, m, a, t, phi)
        if v is None or isinstance(v, str):
            continue
        lr = [lam[i] + (m - i) for i in range(m)]
        hitB = sum(1 for r in RB if sum(r[k] * lr[k] for k in range(m)) % d == 0)
        hitC = sum(1 for r in RC if sum(r[k] * lr[k] for k in range(m)) % d == 0)
        hitB_t = sum(1 for r in RB if sum(r[k] * lr[k] for k in range(m)) % t == 0)

        pred_cero = hitB > 0                 # dim Z_Ghat > dim Z_G = rango
        c1[(pred_cero, v == 0)] += 1
        if pred_cero != (v == 0) and len(mal) < 6:
            mal.append((lam, v, hitB))
        if v != 0:
            c2[abs(v)] += 1
        d1[(hitC > 0, v == 0)] += 1
        d2[(hitB_t > 0, v == 0)] += 1

    n = sum(c1.values())
    ac1 = c1[(True, True)] + c1[(False, False)]
    ad1 = d1[(True, True)] + d1[(False, False)]
    ad2 = d2[(True, True)] + d2[(False, False)]
    print("  C1  FATAL  su dicotomia:  sp_lambda(x) = 0  <=>  algun alpha de B_%d con d | <alpha,lambda+rho>" % m)
    print("        predice cero y ES cero : %4d      predice cero y NO lo es : %4d"
          % (c1[(True, True)], c1[(True, False)]))
    print("        predice no-cero y es 0 : %4d      predice no-cero y no lo es: %4d"
          % (c1[(False, True)], c1[(False, False)]))
    print("        acierta %d de %d   %s" % (ac1, n, "EQUIVALENCIA" if ac1 == n
                                             else "*** falla en %d: %s ***" % (n - ac1, mal[:4])))
    print("  C2  en los no nulos, |valor| : %s   %s"
          % (dict(c2), "todos 1, como pide la mitad de igualdad" if set(c2) <= {1} else "*** hay otros ***"))
    print("  D1  SENUELO raices de C_%d (el grupo) en vez de B_%d (el dual) : %d de %d   %s"
          % (m, m, ad1, n, "NO DISCRIMINA" if ad1 == n else "falla en %d (bien)" % (n - ad1)))
    print("  D2  SENUELO d = t = %d en vez del orden en el adjunto : %d de %d   %s"
          % (t, ad2, n, "no distingue (d = t aqui)" if ad2 == n else "falla en %d (bien)" % (n - ad2)))
    resumen["m=%d" % m] = dict(t=t, d=d, n=int(n), C1={str(k): int(v) for k, v in c1.items()},
                               C1_ac=int(ac1), C2={str(k): int(v) for k, v in c2.items()},
                               D1=int(ad1), D2=int(ad2),
                               mal=[[list(x[0]), x[1], x[2]] for x in mal])

json.dump(resumen, open("npp_question81_DUMP.json", "w"), indent=1)
print("")
print("=" * 100)
print("DONE")
