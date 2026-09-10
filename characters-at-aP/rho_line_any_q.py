# -*- coding: utf-8 -*-
# LA RECTA DE rho EN CUALQUIER RAIZ DE LA UNIDAD.   8 de septiembre de 2026.
#
# POR QUE.  El Teorema 3.1 de note_prasad2 se enuncia para las POTENCIAS de a_P, o sea para
# rho(xi_q) con q = t/gcd(k,t) y t = 2m+2: solo los q que DIVIDEN a t.  Pero al releer su
# demostracion, el unico dato que usa es que z sea una raiz q-esima primitiva de la unidad y que
# (rho, alpha) sea entero.  Nada obliga a q | t.  Si eso es asi, el teorema calcula la anulacion
# en rho(xi) para TODA raiz de la unidad xi, y las potencias de a_P son solo el subcaso q | t.
# Esta gate decide si esa lectura es correcta o si me estoy pasando.
#
# QUE SE MIDE.
#   C0  CONTROL DE INSTRUMENTO.  En un punto REGULAR (donde el bialternante si vale) Jacobi-Trudi
#       y bialternante tienen que coincidir.  Si no, no se mide nada mas.
#   G1  EL CRITERIO, en todo q.   chi_lambda(rho(xi_q)) != 0  <=>  N_q(lambda) = r_q,
#       con N_q(lambda) = #{alpha>0 : q | (lambda+rho, alpha)}  y  r_q = N_q(0).
#       Se separa la columna "q | t" de la columna "q no divide a t", que es lo nuevo.
#   D1  SENUELO 1, el criterio ingenuo: "algun alpha>0 con q | (ell,alpha)" (lo que vale en a_P
#       mismo, donde r_q = 0).  Tiene que FALLAR fuera de los q con r_q = 0.
#   D2  SENUELO 2, la REGLA DE CAPACIDADES de tipo C con d := t/q.  Su derivacion usa m = dq/2 - 1,
#       o sea q | t.  Prediccion: acierta donde q | t y falla donde no.  Si acertara siempre, la
#       distincion que la nota hace entre el criterio y su traduccion no diria nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python rho_line_any_q.py > rho_line_any_q_OUT.txt 2>&1

import sys
from collections import Counter
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

try:
    from mpmath import mp, mpc, exp, pi, mpf
except Exception:
    print("FALTA mpmath.  pip install mpmath")
    sys.exit(2)

mp.dps = 40
SEP = "=" * 92


def dominants(m, top):
    """lambda_1 >= ... >= lambda_m >= 0, con lambda_1 <= top."""
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for rest in rec(k - 1, v):
                yield (v,) + rest
    return rec(m, top)


def roots_C(m):
    """raices positivas de C_m en coordenadas e_i, como vectores enteros."""
    R = []
    for i in range(m):
        for j in range(i + 1, m):
            a = [0] * m; a[i] = 1; a[j] = -1; R.append(tuple(a))
            b = [0] * m; b[i] = 1; b[j] = 1; R.append(tuple(b))
    for i in range(m):
        c = [0] * m; c[i] = 2; R.append(tuple(c))
    return R


def N_q(vec, R, q):
    return sum(1 for a in R if sum(x * y for x, y in zip(vec, a)) % q == 0)


def h_series(xs, N):
    """h_n, n=0..N, del alfabeto {x, 1/x : x in xs}."""
    h = [mpc(0)] * (N + 1)
    h[0] = mpc(1)
    for x in xs:
        for y in (x, 1 / x):
            acc = mpc(0)
            new = [mpc(0)] * (N + 1)
            for n in range(N + 1):
                acc = acc * y + h[n]
                new[n] = acc
            h = new
    return h


def det_c(M):
    n = len(M)
    A = [row[:] for row in M]
    det = mpc(1)
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(A[r][c]))
        if abs(A[p][c]) == 0:
            return mpc(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] * inv
            if abs(f) == 0:
                continue
            for k in range(c, n):
                A[r][k] -= f * A[c][k]
    return det


def chi(lam, m, h):
    """Weyl / Koike-Terada para Sp(2m), variante 'B' de aP_potencias.py."""
    H = lambda n: (h[n] if 0 <= n < len(h) else mpc(0))
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    return det_c(M) / 2


def bialt(lam, xs):
    m = len(lam)
    rho = [m - i for i in range(m)]
    ell = [lam[i] + rho[i] for i in range(m)]
    num = [[x ** e - x ** (-e) for x in xs] for e in ell]
    den = [[x ** e - x ** (-e) for x in xs] for e in rho]
    return det_c(num) / det_c(den)


def cls(x, q):
    return min(x % q, (-x) % q)


# --------------------------------------------------------------------------- C0
print(SEP)
print("C0 -- CONTROL DE INSTRUMENTO: Jacobi-Trudi = bialternante en un punto REGULAR")
peor = mpf(0)
for m in (2, 3, 4):
    xs = [exp(2j * pi * mpf(j) / mpf(97)) for j in range(1, m + 1)]   # 97 primo: regular
    h = h_series(xs, 3 * (8 + m))
    for lam in dominants(m, 4):
        a, b = chi(lam, m, h), bialt(lam, xs)
        peor = max(peor, abs(a - b))
print(f"   OK   peor discrepancia {mp.nstr(peor, 5)}" if peor < mpf("1e-25")
      else f"  FALLA peor discrepancia {mp.nstr(peor, 5)}")

# --------------------------------------------------------------------------- G1 + senuelos
print(SEP)
print("G1 -- EL CRITERIO  N_q(lambda) = r_q  EN TODO q,  y no solo en los q que dividen a t")
print("      D1 = senuelo del criterio ingenuo;  D2 = senuelo de la regla de capacidades")
print()
print(f"{'m':>3} {'q':>4} {'q|t':>4} {'r_q':>4} {'pesos':>7} {'criterio':>9} "
      f"{'D1 ingenuo':>11} {'D2 capac.':>10} {'amb':>4}")

tot = {"ok": 0, "bad": 0, "ok_nodiv": 0, "bad_nodiv": 0,
       "d1": 0, "d1n": 0, "d2": 0, "d2n": 0, "amb": 0}
filas = []
for m in range(2, 7):
    t = 2 * m + 2
    R = roots_C(m)
    rho = [m - i for i in range(m)]
    top = 8 if m <= 4 else 6
    for q in range(2, 15):
        xi = exp(2j * pi / mpf(q))
        xs = [xi ** j for j in range(1, m + 1)]
        h = h_series(xs, 3 * (top + m) + 6)
        r_q = N_q(rho, R, q)
        ok = bad = amb = 0
        d1ok = d1n = d2ok = d2n = 0
        divide = (t % q == 0)
        d = t // q if divide else None
        fixed = {0} | ({q // 2} if q % 2 == 0 else set())
        for lam in dominants(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            v = abs(chi(lam, m, h))
            if mpf("1e-18") < v < mpf("1e-6"):
                amb += 1
                continue
            nz = v >= mpf("1e-6")
            pred = (N_q(ell, R, q) == r_q)
            ok, bad = (ok + 1, bad) if pred == nz else (ok, bad + 1)
            # D1: criterio ingenuo
            p1 = not any(sum(x * y for x, y in zip(ell, a)) % q == 0 for a in R)
            d1ok, d1n = (d1ok + 1, d1n) if p1 == nz else (d1ok, d1n + 1)
            # D2: capacidades con d = t/q (solo tiene sentido escribirlo; se mide igual)
            if divide:
                cnt = Counter(cls(e, q) for e in ell)
                p2 = all(n <= ((d // 2) if c in fixed else d) for c, n in cnt.items())
                d2ok, d2n = (d2ok + 1, d2n) if p2 == nz else (d2ok, d2n + 1)
        tot["ok"] += ok; tot["bad"] += bad; tot["amb"] += amb
        tot["d1"] += d1ok; tot["d1n"] += d1n
        tot["d2"] += d2ok; tot["d2n"] += d2n
        if not divide:
            tot["ok_nodiv"] += ok; tot["bad_nodiv"] += bad
        n = ok + bad
        filas.append((m, q, divide, r_q, n, ok, bad))
        d2s = f"{d2ok}/{d2ok+d2n}" if divide else "  --  "
        print(f"{m:>3} {q:>4} {'si' if divide else 'no':>4} {r_q:>4} {n:>7} "
              f"{ok:>4}/{n:<4} {d1ok:>5}/{d1ok+d1n:<5} {d2s:>10} {amb:>4}")

print()
print(f"     CRITERIO, total ............ {tot['ok']} aciertos, {tot['bad']} fallos")
print(f"     CRITERIO, solo q NO divide a t  {tot['ok_nodiv']} aciertos, {tot['bad_nodiv']} fallos")
print(f"     D1 senuelo ingenuo ......... {tot['d1']} aciertos, {tot['d1n']} fallos"
      f"   (tiene que fallar: {'SI falla' if tot['d1n'] else 'NO FALLA -> no discrimina'})")
print(f"     D2 senuelo capacidades ..... {tot['d2']} aciertos, {tot['d2n']} fallos"
      f"   (medido solo donde q | t)")
print(f"     casos ambiguos numericamente  {tot['amb']}")
print()
if tot["bad"] == 0 and tot["bad_nodiv"] == 0 and tot["d1n"] > 0 and tot["amb"] == 0:
    print("RESULTADO: el criterio vale en TODA raiz de la unidad de la recta de rho, no solo en las")
    print("           potencias de a_P.  El senuelo ingenuo falla, luego la medida discrimina.")
else:
    print("RESULTADO: NO -- revisar antes de escribir nada en el articulo.")
