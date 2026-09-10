# -*- coding: utf-8 -*-
"""EL LEMA DE RESIDUOS, ENUNCIADO SIN TIPO:

   para todo lambda superviviente (N_q(ell) = r_q),   R_ell == R_rho   en Z[X]/(X^q - 1),

donde R_x(X) = sum_{alpha in Phi} X^{(x,alpha)} = char(adjunto)(x) - rango.

El articulo lo demuestra SOLO en tipo C y declara abierto el resto.  Aqui se mide en B, C y G2.
"""
import sys
from collections import Counter
from itertools import product
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def raices_pos(tipo, m):
    """raices POSITIVAS, en coordenadas con (,) estandar (raiz corta de cuadrado 2)."""
    R = []
    if tipo in ("B", "C"):
        for i in range(m):
            for j in range(i + 1, m):
                v = [0] * m; v[i] = 1; v[j] = -1; R.append(tuple(v))
                v = [0] * m; v[i] = 1; v[j] = 1; R.append(tuple(v))
        for i in range(m):
            v = [0] * m
            v[i] = 2 if tipo == "C" else 1
            R.append(tuple(v))
    return R


def R_de(vals, q):
    """multiconjunto de (x,alpha) mod q sobre TODAS las raices (positivas y sus negativas)."""
    c = Counter()
    for v in vals:
        c[v % q] += 1
        c[(-v) % q] += 1
    return c


def prueba(tipo, m, rho, Rpos, red, top, etq):
    t = {"C": 2 * m + 2, "B": 4 * m - 2}[tipo]
    filas = []
    for d in range(1, t + 1):
        if t % d:
            continue
        q = t // d
        if q < 2:
            continue
        vr = [sum(rho[i] * a[i] for i in range(m)) for a in Rpos]
        rq = sum(1 for v in vr if v % q == 0)
        Rrho = R_de(vr, q)
        viv = malos = 0
        for lam in red(m, top):
            ell = [lam[i] + rho[i] for i in range(m)]
            ve = [sum(ell[i] * a[i] for i in range(m)) for a in Rpos]
            if sum(1 for v in ve if v % q == 0) != rq:
                continue
            viv += 1
            if R_de(ve, q) != Rrho:
                malos += 1
        filas.append((q, d, rq, viv, malos))
    print("  %s" % etq)
    for q, d, rq, viv, malos in filas:
        print("     q=%-3d d=%-3d r_q=%-3d  supervivientes %-6d  fallos del lema: %d %s"
              % (q, d, rq, viv, malos, "" if malos == 0 else "  <<<<<< SE CAE"))


def dom_C(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for r in rec(k - 1, v):
                yield (v,) + r
    return rec(m, top)


print("=" * 92)
print("TIPO C  (donde el articulo lo demuestra)")
print("=" * 92)
for m in (2, 3, 4):
    rho = [m - i for i in range(m)]
    prueba("C", m, rho, raices_pos("C", m), dom_C, 8, "C_%d" % m)

print()
print("=" * 92)
print("TIPO B  (que el articulo declara ABIERTO).  rho semientero: se trabaja con y = 2(lambda+rho)")
print("=" * 92)
for m in (2, 3, 4):
    # en B_m, rho = (m-1/2, ..., 1/2).  Duplicamos todo: y_i = 2 ell_i, y las raices se leen igual.
    rho2 = [2 * m - 1 - 2 * i for i in range(m)]       # 2*rho
    Rpos = raices_pos("B", m)
    t = 4 * m - 2
    for d in [dd for dd in range(1, t + 1) if t % dd == 0]:
        q = t // d
        if q < 2:
            continue
        vr = [sum(rho2[i] * a[i] for i in range(m)) for a in Rpos]
        rq = sum(1 for v in vr if v % q == 0)
        Rrho = R_de(vr, q)
        viv = malos = 0
        # pesos: lambda entero O semientero (espinoriales) -> 2*lambda entero con paridad fija
        for par in (0, 1):
            for lam2 in product(range(par, 2 * 8 + 1, 2), repeat=m):
                if any(lam2[i] < lam2[i + 1] for i in range(m - 1)):
                    continue
                y = [lam2[i] + rho2[i] for i in range(m)]
                ve = [sum(y[i] * a[i] for i in range(m)) for a in Rpos]
                if sum(1 for v in ve if v % q == 0) != rq:
                    continue
                viv += 1
                if R_de(ve, q) != Rrho:
                    malos += 1
        print("  B_%d  q=%-3d d=%-3d r_q=%-3d  supervivientes %-6d  fallos: %d %s"
              % (m, q, d, rq, viv, malos, "" if malos == 0 else "  <<<<<< SE CAE"))

print()
print("=" * 92)
print("G2.  Seis raices positivas con (x,alpha) = u, v, u+v, 2u+v, 3u+v, 3u+2v;  t=12")
print("=" * 92)
for d in (1, 2, 3, 4, 6, 12):
    q = 12 // d
    if q < 2:
        continue
    vr = [1, 3, 4, 5, 6, 9]
    rq = sum(1 for v in vr if v % q == 0)
    Rrho = R_de(vr, q)
    viv = malos = 0
    for a in range(0, 3 * q + 2):
        for b in range(0, 3 * q + 2):
            u, v = a + 1, 3 * (b + 1)
            ve = [u, v, u + v, 2 * u + v, 3 * u + v, 3 * u + 2 * v]
            if sum(1 for w in ve if w % q == 0) != rq:
                continue
            viv += 1
            if R_de(ve, q) != Rrho:
                malos += 1
    print("  G2   q=%-3d d=%-3d r_q=%-3d  supervivientes %-6d  fallos: %d %s"
          % (q, d, rq, viv, malos, "" if malos == 0 else "  <<<<<< SE CAE"))
