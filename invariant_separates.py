# -*- coding: utf-8 -*-
# Authors: Carles Marin, Claude (AI assistant).
# OJO: importa de thm31_audit, que se ejecuta al importarse; por eso su _OUT.txt lleva delante la
# salida de aquel.  Es a proposito: las dos auditorias comparten el lector del enunciado.
"""Prop. 3.11 auditada aparte: (multiconjunto {d1,d2,d3}, eps) SEPARA los valores no nulos.

Dos direcciones:
  (=>) mismo dato  => mismo valor            (trivial desde la formula, se comprueba igual)
  (<=) mismo valor => mismo dato             (esto es lo que hay que romper si es falso)

El valor se representa por su vector de evaluaciones en varios z racionales, a 50 digitos.
"""
import mpmath as mp
from thm_main_independent import closed_form, partitions

mp.mp.dps = 50
ZS = [mp.mpf(3) / 2, mp.mpf(7) / 5, mp.mpf(11) / 4, mp.mpf(9) / 7, mp.mpf(13) / 6]


def datum(beta, t):
    N = len(beta)
    cls = {}
    for j, b in enumerate(beta):
        cls.setdefault(b % t, []).append((b, j))
    if len(cls) < t:
        return None
    exc = sorted(i for i in cls if len(cls[i]) >= 2)
    if len(exc) == 2:
        A = sorted(cls[exc[0]], reverse=True)
        B = sorted(cls[exc[1]], reverse=True)
    else:
        p, q, r = sorted(cls[exc[0]], reverse=True)
        A, B = [p, q], [q, r]
    d1 = A[0][0] - A[1][0]
    d2 = B[0][0] - B[1][0]
    dt3 = A[0][0] + A[1][0] - B[0][0] - B[1][0]
    if dt3 == 0:
        return None
    return tuple(sorted((d1, d2, abs(dt3))))


print("  t  |lam|<=  no nulas  datos distintos  |  mismo dato/valor distinto  |  mismo valor/dato distinto")
tot_a = tot_b = tot_n = 0
for t, MS in ((2, 14), (3, 13), (4, 12), (5, 11)):
    N = t + 2
    by_datum, by_value = {}, {}
    n = 0
    a_fail = b_fail = 0
    for lam in partitions(MS, N):
        L = list(lam) + [0] * (N - len(lam))
        beta = [L[j] + N - 1 - j for j in range(N)]
        d = datum(beta, t)
        if d is None:
            continue
        vals = []
        for z in ZS:
            v, tag = closed_form(beta, t, z)
            vals.append(v)
        if tag != "formula":
            continue
        n += 1
        eps_sign = 1 if vals[0] > 0 else -1
        key = (d, eps_sign)
        vkey = tuple(mp.nstr(v, 30) for v in vals)
        if key in by_datum and by_datum[key] != vkey:
            a_fail += 1
        by_datum.setdefault(key, vkey)
        if vkey in by_value and by_value[vkey] != key:
            b_fail += 1
            if b_fail <= 1:
                print("      COLISION: %s  vs  %s" % (by_value[vkey], key))
        by_value.setdefault(vkey, key)
    tot_a += a_fail
    tot_n += n
    tot_b += b_fail
    print("  %2d %6d %9d %16d  | %26d | %26d" % (t, MS, n, len(by_datum), a_fail, b_fail))
print()
print("TOTAL: %d formas no nulas, 4 configuraciones" % tot_n)
print("TOTAL: mismo dato con valor distinto = %d ;  mismo valor con dato distinto = %d" % (tot_a, tot_b))

# control: quitar el signo del dato TIENE que producir colisiones
print()
coll = 0
for t, MS in ((3, 13), (4, 12)):
    N = t + 2
    by_value = {}
    for lam in partitions(MS, N):
        L = list(lam) + [0] * (N - len(lam))
        beta = [L[j] + N - 1 - j for j in range(N)]
        d = datum(beta, t)
        if d is None:
            continue
        v, tag = closed_form(beta, t, ZS[0])
        if tag != "formula":
            continue
        vkey = mp.nstr(v, 30)
        if vkey in by_value and by_value[vkey] != d:
            pass
        by_value.setdefault(vkey, d)
    # ahora al reves: mismo multiconjunto SIN signo, valores opuestos
    seen = {}
    for lam in partitions(MS, N):
        L = list(lam) + [0] * (N - len(lam))
        beta = [L[j] + N - 1 - j for j in range(N)]
        d = datum(beta, t)
        if d is None:
            continue
        v, tag = closed_form(beta, t, ZS[0])
        if tag != "formula":
            continue
        if d in seen and mp.nstr(seen[d], 30) != mp.nstr(v, 30):
            coll += 1
        seen.setdefault(d, v)
print("CONTROL: dejando el signo FUERA del dato hay %d colisiones (si fuera 0, el signo sobraria)" % coll)
