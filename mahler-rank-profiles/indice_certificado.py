# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: certificado de completitud por el indice, piso 1, caso semisimple.
#      H_seed = subalgebra graduada generada por los soportes exactos (orden_por_caracteres.py: capas 0..2p-1) cerrada
#      por productos. H_seed,k contenido en H_k, luego
#        delta_seed - delta_A = sum_k (dim H_k - dim H_seed,k) >= 0,   delta_A = p d - 2 + Z(pn) + Z(n)  (indice cerrado).
#      Si delta_seed = delta_A (con cola llena certificada), el semilla es TODO el graduado: afirmacion (K) probada en el
#      caso, mas fuerte que e exacto. Como las capas 0..2p-1 del semilla suman ya delta_A por construccion, el
#      certificado equivale a: el cierre llena toda capa k >= 2p. Cola: una ventana llena de longitud p-1 basta
#      (multiplicar por el trivial de la capa p-1 lleva capa llena k a k+p-1).
#      Se registra por caso: defecto Delta = delta_seed - delta_A, primera capa defectuosa, y si la cola se certifica.
#      Casos: los 90 de orden_por_caracteres (e medido) + controles compuestos (63,5), (65,7), (105,11)
#      (sin e medido; FUERA de los casos donde la formula del indice se contrasto: se marcan aparte).
# FECHA: 1-oct-2026
# Uso (docker, /work): sage -python indice_certificado.py
import re

from sage.all import GF, QQ, DirichletGroup, Integers, bernoulli, euler_phi

import orden_por_caracteres as OC  # noqa: F401  (su bloque final tambien corre: 90 casos, rapido)


def soportes(n, p, Dmax):
    """copia de OC.orden devolviendo los soportes S (misma construccion)."""
    lam = Integers(n).unit_group_exponent()
    f = 1
    while (p ** f - 1) % lam:
        f += 1
    F = GF(p ** f, 'g')
    z = F.multiplicative_generator() ** ((p ** f - 1) // lam)
    G = list(DirichletGroup(n, F, zeta=z, zeta_order=lam))
    idx = {chi: i for i, chi in enumerate(G)}
    triv = next(i for i, c in enumerate(G) if c.is_trivial())
    mult = [[idx[G[i] * G[j]] for j in range(len(G))] for i in range(len(G))]
    npar = sum(1 for c in G if c(-1) == 1)
    nimp = len(G) - npar
    S = [set() for _ in range(Dmax + 1)]
    S[0] = {triv}
    S[p - 1] = {triv}
    prim = {i: G[i].primitive_character() for i in range(len(G))}
    for j in range(p - 1):
        dd = p + j
        par = (-1) ** (j + 1)
        for i, c in enumerate(G):
            if c(-1) != par:
                continue
            if i == triv:
                if j % 2 == 1 and (j == p - 2 or F((1 - 2 ** j) * QQ(bernoulli(j + 1))) != 0):
                    S[dd].add(i)
                continue
            if (1 - F(2) ** j * c(2)) * OC.B_prim(prim[i], j + 1, F) != 0:
                S[dd].add(i)
    S[2 * p - 1] = set(S[p])
    cambio = True
    while cambio:
        cambio = False
        for a in range(p - 1, Dmax + 1):
            if not S[a]:
                continue
            for b in range(a, Dmax + 1 - a):
                if not S[b]:
                    continue
                nuevo = {mult[x][y] for x in S[a] for y in S[b]}
                if not nuevo <= S[a + b]:
                    S[a + b] |= nuevo
                    cambio = True
    llena = [len(S[k]) == (npar if k % 2 == 0 else nimp) for k in range(Dmax + 1)]
    return S, llena, npar, nimp


def certificado(n, p):
    d = euler_phi(n) // 2
    Dmax = 6 * p + 4
    S, llena, npar, nimp = soportes(n, p, Dmax)
    dim = lambda k: npar if k % 2 == 0 else nimp  # noqa: E731
    Zpn = sum(dim(p + j) - len(S[p + j]) for j in range(p - 1))
    Zn = dim(2 * p - 1) - len(S[2 * p - 1])
    delta_A = p * d - 2 + Zpn + Zn
    # cola certificada: primera ventana llena [k0, k0+p-1) con k0 >= 2p dentro de Dmax
    k0 = next((k for k in range(2 * p, Dmax - p + 2) if all(llena[k:k + p - 1])), None)
    if k0 is None:
        return d, delta_A, None, None, None
    delta_seed = sum(dim(k) - len(S[k]) for k in range(k0))
    defectuosas = [k for k in range(2 * p, k0) if not llena[k]]
    return d, delta_A, delta_seed, (defectuosas[0] if defectuosas else None), k0


PAT = re.compile(r"n=\s*(\d+) p=\s*(\d+) d=\s*(\d+) .*?e=(\d+|None)")
casos = []
for nom in ('conjetura_C_lote_OUT.txt', 'conjetura_C_p17_OUT.txt'):
    for ln in open('/work/' + nom, encoding='utf-8'):
        m = PAT.search(ln)
        if m and m[4] != 'None':
            casos.append((int(m[1]), int(m[2])))
print("=" * 30, "certificado de indice", flush=True)
tot = cert = 0
fallos = []
for n, p in casos:
    if euler_phi(n) % p == 0:
        continue
    d, dA, ds, k1, k0 = certificado(n, p)
    tot += 1
    ok = ds is not None and ds == dA
    cert += ok
    if not ok:
        fallos.append((n, p, dA, ds, k1))
    print("n=%3d p=%2d d=%2d: delta_A=%d delta_seed=%s Delta=%s  primera defectuosa>=2p: %s  cola desde %s  %s" % (
        n, p, d, dA, ds, (ds - dA) if ds is not None else '-', k1, k0, 'CERTIFICADO' if ok else 'NO'), flush=True)
print("CERTIFICADO (K) en %d de %d casos; fallos: %s" % (cert, tot, fallos), flush=True)
print("=" * 30, "controles compuestos (fuera de la zona contrastada de la formula del indice)", flush=True)
for n, p in [(63, 5), (65, 7), (105, 11)]:
    d, dA, ds, k1, k0 = certificado(n, p)
    print("n=%3d p=%2d d=%2d: delta_A=%d delta_seed=%s Delta=%s  primera defectuosa>=2p: %s  cola desde %s" % (
        n, p, d, dA, ds, (ds - dA) if ds is not None else '-', k1, k0), flush=True)
