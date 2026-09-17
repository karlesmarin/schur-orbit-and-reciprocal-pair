# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# PERFIL DE PROFUNDIDAD de A en p:  B_k = imagen de A en F_p[y]/(Phi_n(y)^k), y = zeta_q, k <= E = phi(p^v)
#   W_i = dim B_{i+1} - dim B_i      (W_0 = 1 por H2', W_1 = dim F_p[<p>].S por el teorema del reparto)
# y lo cruza con el reticulo:  e (exponente de c en p, todos los P iguales) y b = v_p([A:c]).
#   H-e:  e = primer k tal que W_i = phi(n)/2 (total) para todo i en [k, E-1]      (si e <= E-1)
#   H-b:  b = sum_{i<e} W_i                                                          (teorema si e <= E)
# Casillas con e > E o E*phi(n) grande: se cuentan aparte, no se miden.
#
# Uso: python perfil_W.py [DMAX]   (DMAX = tope de E*phi(n), por coste)
#
# Authors: Carles Marin, Claude (AI assistant).
import contextlib
import io
import os
import sys
from collections import Counter
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_argv = sys.argv
sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    from reparto import leer  # noqa: E402
sys.argv = _argv
from mc2_combinatorio import factor_p  # noqa: E402
from reparto_p2 import cyclo, polmul, phi, S_mod, LOGS_RETICULO, LOGS_H2  # noqa: E402


def dims_B(q, p, m, n, K):
    """[dim B_1, ..., dim B_K] con B_k en F_p[y]/(Phi_n^k).  Se calcula en k = K y se proyecta."""
    phin = [c % p for c in cyclo(n)]
    mod = [1]
    for _ in range(K):
        mod = polmul(mod, phin, p)
    D = len(mod) - 1
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
    coef = [one]
    for j in range(1, m + 1):
        Eta = (pots[j % q] + pots[(-j) % q]) % p
        new = [np.zeros(D, dtype=np.int64) for _ in range(len(coef) + 1)]
        for i, c in enumerate(coef):
            new[i + 1] = (new[i + 1] + c) % p
            new[i] = (new[i] - Eta @ c) % p
        coef = new
    gens = coef[:-1]
    inv = [0] + [pow(x, -1, p) for x in range(1, p)]
    piv = {}
    base = []

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

    def matof(v):
        M = np.zeros((D, D), dtype=np.int64)
        for k in range(D):
            if v[k]:
                M = (M + int(v[k]) * pots[k]) % p
        return M
    add(one)
    base.append(one)
    mats = [matof(g) for g in gens]
    for g in gens:
        if add(g) is not None:
            base.append(g)
    cola = list(base)
    while cola:
        v = cola.pop()
        for M in mats:
            w = M @ v % p
            if add(w) is not None:
                cola.append(w)
    B = np.array([base_v for base_v in piv.values()], dtype=np.int64) if piv else np.zeros((0, D), dtype=np.int64)
    # proyecciones a F_p[y]/(Phi_n^k): reducir cada vector de la base modulo Phi_n^k y medir rango
    out = []
    for k in range(1, K + 1):
        modk = [1]
        for _ in range(k):
            modk = polmul(modk, phin, p)
        Dk = len(modk) - 1
        filas = []
        for vrow in B:
            poly = [int(x) for x in vrow]
            # resto de poly modulo modk (monico)
            for deg in range(len(poly) - 1, Dk - 1, -1):
                c = poly[deg]
                if c:
                    for t in range(Dk + 1):
                        poly[deg - Dk + t] = (poly[deg - Dk + t] - c * modk[t]) % p
            filas.append(poly[:Dk])
        out.append(rango(filas, p))
    return out


def rango(filas, p):
    if not filas:
        return 0
    M = np.array(filas, dtype=np.int64) % p
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r] * pow(int(M[r, c]), -1, p)) % p
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c] * M[r]) % p
        r += 1
        if r == rows:
            break
    return r


def filas_reticulo():
    filas = {}
    for f in leer(LOGS_RETICULO):
        for p in set(t[0] for t in f["fac"]):
            forma = [(ff, ee) for (pp, ff, ee) in f["fac"] if pp == p]
            a = vp(f["idx"], p)
            b = vp(f["acx"], p)
            filas[(f["m"], f["q"], p)] = (a, b, forma)
    for ruta in LOGS_H2:
        try:
            for ln in open(ruta, encoding="utf-8", errors="replace"):
                ln = ln.strip()
                if ln.startswith("H2|"):
                    _, m, q, p, a, b, dm, forma, c1, c2 = ln.split("|")
                    fo = [tuple(int(x) for x in t.split(".")) for t in forma.split("+")]
                    filas[(int(m), int(q), int(p))] = (int(a), int(b), fo)
        except FileNotFoundError:
            pass
    return filas


def vp(x, p):
    v = 0
    while x and x % p == 0:
        x //= p
        v += 1
    return v


def main():
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 240
    tab = Counter()
    detalle = []
    for (m, q, p), (a, b, forma) in sorted(filas_reticulo().items()):
        v, n = factor_p(q, p)
        if n < 3:
            continue
        es = set(ee for (ff, ee) in forma)
        if len(es) != 1:
            tab["exponentes distintos entre P"] += 1
            continue
        e = es.pop()
        E = phi(p ** v)
        full = phi(n) // 2
        if e > E:
            tab["e > E (fuera del regimen de caracteristica p)"] += 1
            detalle.append(("e>E", m, q, p, v, n, E, forma, a, b, None))
            continue
        K = min(E, e + 1) if e < E else E
        if K * phi(n) > DMAX:
            tab["demasiado grande (no medida)"] += 1
            continue
        dB = dims_B(q, p, m, n, K)
        W = [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, len(dB))]
        hb = (b == sum(W[:e]))
        # H-e: con K = e+1 se ve W_e; e es el primer k con W_i total en [k, K-1] -- y W_{e-1} no total
        if e < K:
            he = (W[e] == full) and (e == 0 or W[e - 1] != full or e == 1 and W[0] == full)
        else:
            he = None
        tab[("H-b", hb)] += 1
        tab[("H-e", he)] += 1
        detalle.append(("ok" if hb and he is not False else "FALLA", m, q, p, v, n, E, forma, a, b, W))
    print("resumen:", dict(tab))
    for t in detalle:
        print("   %-5s m=%2d q=%3d p=%d v=%d q'=%2d E=%2d forma=%s a=%d b=%d  W=%s" % t)


if __name__ == "__main__":
    main()
