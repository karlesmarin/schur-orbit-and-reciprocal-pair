# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Prediccion escrita antes de los datos -- el reparto con c_p = prod P^2.
#   R1  b (reticulo) == beta = dim_Fp imagen de A en F_p[y]/(Phi_n(y)^2), y = zeta_q  (sin reticulo)
#   R2  (un solo P) beta - 1 in {1, f-1}
#   R3  beta - 1 == 1  <=>  S(pc) = lambda.S(c) mod p para todo c (2c != 0), lambda en F_p^*
# Filas: reparto.leer (LOGS_RETICULO: celdas_m6/7/8, ley_conductor_*) + h2_OUT.txt, r2_OUT.txt.  Las 11 (3,2) que sugirieron la
# prediccion (m in {3,6,7,8}, q <= 72, vistas antes de la prediccion) se cuentan APARTE.
#
# Uso: python reparto_p2.py   (desde este directorio)
#
# Authors: Carles Marin, Claude (AI assistant).
import contextlib
import io
import os
import sys
from math import gcd
from collections import Counter, defaultdict
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_argv = sys.argv
sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    from reparto import leer  # noqa: E402
sys.argv = _argv
from mc2_combinatorio import primos, factor_p  # noqa: E402

# Los logs del reticulo que se leen, por NOMBRE y no por glob: en este directorio hay otros *_OUT.txt
# (las salidas de los analisis) y un glob los recogeria.  Estos cinco, con h2_OUT.txt y r2_OUT.txt
# (filas H2|, que se leen aparte), dan todas las filas contadas, y quitar cualquiera de los siete
# baja el recuento.  celdas_m2/m3_OUT.txt (que lee datos_indice.py) repiten filas de los logs de
# ley_conductor con los mismos valores: anadirlos no cambia nada.
AQUI = os.path.dirname(os.path.abspath(__file__))
LOGS_RETICULO = [os.path.join(AQUI, n) for n in (
    "celdas_m6_OUT.txt", "celdas_m7_OUT.txt", "celdas_m8_OUT.txt",
    "ley_conductor_ancho_OUT.txt", "ley_conductor_repro_OUT.txt")]
LOGS_H2 = [os.path.join(AQUI, n) for n in ("h2_OUT.txt", "r2_OUT.txt")]

SUGIRIERON = {(3, 21), (3, 28), (3, 35), (3, 56), (3, 63), (6, 28), (6, 56), (7, 28), (7, 56), (8, 36), (8, 72)}


def phi(n):
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


_cyc = {}


def cyclo(n):
    if n in _cyc:
        return _cyc[n]
    p = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            b = cyclo(d)
            a = p[:]
            out = [0] * (len(a) - len(b) + 1)
            for k in range(len(out) - 1, -1, -1):
                c = a[k + len(b) - 1]
                out[k] = c
                for j, yv in enumerate(b):
                    a[k + j] -= c * yv
            assert all(v == 0 for v in a)
            p = out
    _cyc[n] = p
    return p


def polmul(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, yv in enumerate(b):
                out[i + j] = (out[i + j] + x * yv) % p
    return out


def beta(q, p, m, n):
    """dim_Fp de F_p[e_1..e_m] dentro de F_p[y]/(Phi_n(y)^2), y = zeta_q."""
    mod = polmul([c % p for c in cyclo(n)], [c % p for c in cyclo(n)], p)
    D = len(mod) - 1
    # matriz de multiplicacion por y
    Y = np.zeros((D, D), dtype=np.int64)
    for k in range(D):
        if k + 1 < D:
            Y[k + 1, k] = 1
        else:
            for t in range(D):
                Y[t, k] = (-mod[t]) % p
    pots = [np.eye(D, dtype=np.int64)]
    for _ in range(1, q):
        pots.append((Y @ pots[-1]) % p)
    one = np.zeros(D, dtype=np.int64)
    one[0] = 1
    coef = [one]  # coeficientes de prod (T - eta_j), como vectores
    for j in range(1, m + 1):
        Eta = (pots[j % q] + pots[(-j) % q]) % p
        new = [np.zeros(D, dtype=np.int64) for _ in range(len(coef) + 1)]
        for i, c in enumerate(coef):
            new[i + 1] = (new[i + 1] + c) % p
            new[i] = (new[i] - Eta @ c) % p
        coef = new
    gens = coef[:-1]

    def matof(v):
        M = np.zeros((D, D), dtype=np.int64)
        for k in range(D):
            if v[k]:
                M = (M + int(v[k]) * pots[k]) % p
        return M
    inv = [0] + [pow(x, -1, p) for x in range(1, p)]
    piv = {}

    def add(v):
        v = v % p
        for col, row in piv.items():
            if v[col]:
                v = (v - v[col] * row) % p
        nz = np.nonzero(v)[0]
        if len(nz) == 0:
            return None
        col = int(nz[0])
        v = (v * inv[int(v[col])]) % p
        for c2 in list(piv):
            r = piv[c2]
            if r[col]:
                piv[c2] = (r - r[col] * v) % p
        piv[col] = v
        return v
    add(one)
    mats = []
    for g in gens:
        if add(g) is not None:
            mats.append(matof(g))
    for g in gens:
        mats.append(matof(g))
    cola = [r.copy() for r in piv.values()]
    while cola:
        v = cola.pop()
        for M in mats:
            w = add(M @ v)
            if w is not None:
                cola.append(w)
    return len(piv)


def S_mod(m, n, p):
    out = {}
    for c in range(1, n):
        if (2 * c) % n == 0:
            continue
        s = sum(j for j in range(1, m + 1) if j % n == c) - sum(j for j in range(1, m + 1) if j % n == (n - c) % n)
        out[c] = s % p
    return out


def autovector(m, n, p):
    S = S_mod(m, n, p)
    if all(v == 0 for v in S.values()):
        return "S=0"
    lam = None
    for c, s in S.items():
        t = S[(p * c) % n]
        if s == 0:
            if t != 0:
                return False
            continue
        l = (t * pow(s, -1, p)) % p
        if lam is None:
            lam = l
        elif l != lam:
            return False
    return True


def main():
    filas = {}
    for f in leer(LOGS_RETICULO):
        for p in set(t[0] for t in f["fac"]):
            forma = [(ff, ee) for (pp, ff, ee) in f["fac"] if pp == p]
            a = ZZv(f["idx"], p)
            b = ZZv(f["acx"], p)
            filas[(f["m"], f["q"], p)] = (a, b, forma)
    for ruta in LOGS_H2:
      try:
        for ln in open(ruta, encoding="utf-8", errors="replace"):
            ln = ln.strip()
            while ln.startswith("sage:") or ln.startswith("....:"):
                ln = ln[5:].strip()
            if ln.startswith("H2|"):
                _, m, q, p, a, b, dm, forma, c1, c2 = ln.split("|")
                fo = [tuple(int(x) for x in t.split(".")) for t in forma.split("+")]
                filas[(int(m), int(q), int(p))] = (int(a), int(b), fo)
      except FileNotFoundError:
        pass
    r1 = Counter()
    r2 = Counter()
    r3 = Counter()
    detalle = []
    fuera = Counter()
    for (m, q, p), (a, b, forma) in sorted(filas.items()):
        v, n = factor_p(q, p)
        if n < 3:
            continue
        E = phi(p ** v)
        todas_e2 = all(ee == 2 for (ff, ee) in forma)
        av = autovector(m, n, p)
        if not (todas_e2 and E >= 2):
            if forma and av is True:
                fuera[tuple(sorted(forma))] += 1
            continue
        grupo = "sugirieron" if (m, q) in SUGIRIERON and forma == [(3, 2)] else "nuevas"
        bt = beta(q, p, m, n)
        r1[(grupo, b == bt)] += 1
        f = forma[0][0]
        un_P = len(forma) == 1
        dV = bt - 1
        if un_P:
            r2[(grupo, dV in (1, f - 1))] += 1
            r3[(grupo, (dV == 1) == (av is True), av)] += 1
        detalle.append((grupo, m, q, p, v, n, E, forma, a, b, bt, av))
    print("R1  b(reticulo) == beta(F_p[y]/Phi_n^2):", dict(r1))
    print("R2  un solo P, dim V in {1, f-1}:", dict(r2))
    print("R3  (dim V == 1) <=> S autovector de xp:", dict(r3))
    print("detalle (grupo, m, q, p, v, q', E, forma, a, b, beta, S autovector):")
    for t in detalle:
        print("   ", t)
    print("casillas con S autovector FUERA de la hipotesis c_p = P^2 (por forma):", dict(fuera))


def ZZv(x, p):
    v = 0
    while x and x % p == 0:
        x //= p
        v += 1
    return v


if __name__ == "__main__":
    main()
