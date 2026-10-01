# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: instrumento de Mahler v3, perfiles mod P SIN numeros de Bernoulli (masas de clase;
#      derivacion abajo). Solo da FILAS PIVOTE (perfiles mod P): Phi_2 y los escalares no se
#      transportan sin su normalizacion.
#   Objeto. P_chi(Y) = sum_{u<n} chi(u) Y^u, F_chi = P_chi/(1 - Y^n), H_chi(Y) = -2 (F_chi(Y) - chi(2) F_chi(Y^2));
#      H_chi(e^s) = sum_j beta_j s^j/j!, beta_j = 2 (1 - 2^j chi(2)) B_{j+1,chi}/(j+1).
#   Por M0, sum_i (A/tau)_{J,i} C(x,i) = [t^J] H_chi(E(t)^x). Columnas x = 0..p^a-1: el cambio C(x,i) es unitriangular
#      entero; filas: T = E(t) - 1 = t + O(t^2) con coeficientes p-enteros, cambio unitriangular. Luego
#      B_{J,x} = [T^J] H_chi((1+T)^x) tiene las MISMAS filas pivote que A (se comprueba en 'control', no se supone).
#   Masas. M = p^b >= L+1. Con Z = Y^M y S_a(Z) = sum_{r<n} chi(a+Mr) Z^r: F_chi = sum_{a<M} Y^a [S_a(Z)/(1-Z)]/D(Z),
#      D(Z) = sum_{r<n} Z^r, D(1) = n unidad, S_a(1) = 0 (M invertible mod n). En caracteristica p, Z = 1 + (Y-1)^M, asi
#      que F_chi(Y) = sum_a w_a Y^a mod (P, (Y-1)^M) con
#        (8) w_a = -S_a'(1)/n = -(1/n) sum_{r<n} r chi(a+Mr),
#        (9) w_a = -chi(M) (B_chi/n + C_chi(c)),  B_chi = sum_{u<n} u chi(u),  C_chi(c) = sum_{u<c} chi(u),  c = a/M mod n.
#      Y = (1+T)^x: B_{J,x} = -2 sum_a w_a [C(ax mod M, J) - chi(2) C(2ax mod M, J)] mod P.  Coste O(n + M + p^a M L).
#   Primo: P = (p, g(zeta_m)), g = mahler_padico.primo_sobre_p(m, p) (convencion del 2-oct), zeta_m -> raiz de g.
#   El factor global -2 chi(M) no cambia pivotes y se omite en perfil (se conserva en masas_directas / masas_cerradas).
# FECHA: 2-oct-2026
# Uso (docker, /work): sage -python mahler_masas.py control | masas | lambda_alta
import sys
import time

import numpy as np

sys.path.insert(0, '/work')
import mahler_padico as MP  # noqa: E402
from sage.all import GF, Integers, matrix, kronecker_character  # noqa: E402


def tabla_potencias(g, m):
    """coordenadas (base 1, x, ..., x^{f-1} de F_p[x]/(g)) de x^e, e = 0..m-1; y matriz de multiplicacion por x^e."""
    f = g.degree()
    X = g.parent().gen()
    pot, c = [], g.parent()(1)
    for _ in range(m):
        pot.append(([int(v) for v in c.list()] + [0] * f)[:f])
        c = (c * X) % g
    return np.array(pot, dtype=np.int64)


def mult_matriz(g, e, m):
    """matriz f x f de la multiplicacion por x^e en F_p[x]/(g) (columnas: imagenes de la base)."""
    f = g.degree()
    X = g.parent().gen()
    cols = []
    for k in range(f):
        c = (X ** (k + e % m)) % g
        cols.append(([int(v) for v in c.list()] + [0] * f)[:f])
    return np.array(cols, dtype=np.int64).T


def binomiales(M, L, p):
    """Bin[y, J] = C(y, J) mod p, y < M, J <= L."""
    Bn = np.zeros((M, L + 1), dtype=np.int64)
    Bn[0, 0] = 1
    for y in range(1, M):
        Bn[y, 0] = 1
        Bn[y, 1:] = (Bn[y - 1, 1:] + Bn[y - 1, :-1]) % p
    return Bn


def masas_cerradas(ex, m, n, p, M, g):
    """v_a = B_chi/n + C_chi(a/M mod n) en coordenadas (M x f); w_a = -chi(M) v_a. ex[u] = e con chi(u) = zeta_m^e,
    o -1 si chi(u) = 0."""
    f = g.degree()
    pot = tabla_potencias(g, m)
    ch = np.zeros((n, f), dtype=np.int64)
    for u in range(n):
        if ex[u] >= 0:
            ch[u] = pot[ex[u]]
    Bchi = (np.arange(n, dtype=np.int64)[:, None] % p * ch).sum(axis=0) % p
    ninv = pow(int(n), -1, int(p))
    C = np.zeros((n, f), dtype=np.int64)
    C[1:] = np.cumsum(ch[:-1], axis=0) % p
    Minv = pow(int(M), -1, int(n))
    c = (np.arange(M, dtype=np.int64) * Minv) % n
    return (Bchi * ninv + C[c]) % p


def masas_directas(ex, m, n, p, M, g):
    """(8): n w_a / (-1) = sum_{r<n} r chi(a+Mr) -> v'_a = (1/n) sum_r r chi(a+Mr), para comparar con chi(M) v_a."""
    f = g.degree()
    pot = tabla_potencias(g, m)
    ninv = pow(int(n), -1, int(p))
    out = np.zeros((M, f), dtype=np.int64)
    for a in range(M):
        acc = np.zeros(f, dtype=np.int64)
        for r in range(n):
            e = ex[(a + M * r) % n]
            if e >= 0:
                acc += r * pot[e]
        out[a] = acc % p * ninv % p
    return out


def matriz_B(ex, e2, m, n, p, a, L):
    """B (L+1) x p^a sobre F_q (q = p^f), salvo el factor -2 chi(M)."""
    g = MP.primo_sobre_p(m, p)
    M = modulo_M(p, L)
    v = masas_cerradas(ex, m, n, p, M, g)
    v2 = (v @ mult_matriz(g, e2, m).T) % p  # chi(2) v_a
    return matriz_desde_masas(v, v2, g, p, a, L)


def modulo_M(p, L):
    """M = p^b >= L+1: (1+T)^M = 1 mod (p, T^{L+1})."""
    M = p
    while M < L + 1:
        M *= p
    return M


def matriz_desde_masas(v, v2, g, p, a, L):
    """B_{J,x} = [T^J] sum_a (v_a (1+T)^{ax} - v2_a (1+T)^{2ax}) mod P, J <= L, x < p^a; v, v2: M x f coordenadas en
    F_p[x]/(g). Comun a los caracteres (v2 = chi(2) v) y a zeta (v2 = 0)."""
    f = g.degree()
    M = v.shape[0]
    P = p ** a
    ar = np.arange(M, dtype=np.int64)
    Bn = binomiales(M, L, p)
    # producto por BLAS en float64: exacto porque cada suma es < M (p-1)^2 < 2^53
    assert M * (p - 1) ** 2 < 2 ** 53
    Bnf = Bn.astype(np.float64)
    Bk = [np.zeros((P, L + 1), dtype=np.int64) for _ in range(f)]  # P x (L+1)
    bloque = max(1, min(P, (1 << 25) // (M * f)))  # columnas x por bloque (memoria ~ 256 MB)
    for x0 in range(0, P, bloque):
        xs = range(x0, min(P, x0 + bloque))
        mass = np.zeros((len(xs), M, f), dtype=np.int64)
        for i, x in enumerate(xs):
            np.add.at(mass[i], (ar * x) % M, v)
            np.add.at(mass[i], (2 * ar * x) % M, -v2)
        mass %= p
        for k in range(f):
            Bk[k][x0:x0 + len(xs)] = np.rint(mass[:, :, k].astype(np.float64) @ Bnf).astype(np.int64) % p
    Fq = GF(p ** f, 'y', modulus=g) if f > 1 else GF(p)
    if f == 1:
        return matrix(Fq, Bk[0].T.tolist())
    ent = []
    for J in range(L + 1):
        for x in range(P):
            ent.append(Fq([int(Bk[k][x, J]) for k in range(f)]))
    return matrix(Fq, L + 1, P, ent)


exponentes = MP.exponentes  # mismas orbitas y mismo orden que mahler_padico.perfil
orbitas = MP.orbitas


def perfil(n, p, a, L):
    """como mahler_padico.perfil: (d, [(chi(2), chi(p), tamano, filas pivote)], #primitivos)."""
    sal = []
    obs, npr = orbitas(n, p)
    for chi, orb, m, Km, chim in obs:
        ex = exponentes(chim, Km, n)
        piv = matriz_B(ex, ex[2], m, n, p, a, L).pivot_rows()
        sal.append((str(chi(2)), str(chi(p)), len(orb), piv))
    return Integers(n).unit_group_order() // 2, sal, npr


def _control():
    from mahler_rapido import CONTROL
    tot = ok = 0
    difs = 0
    for n, p, a, med in CONTROL:
        P = p ** a
        L = max(c - P for c in med)
        t0 = time.time()
        d, sal, npr = perfil(n, p, a, L)
        dt = time.time() - t0
        # fila a fila frente a mahler_padico (mismo primo, mismo orden de orbitas)
        d2, sal2, _ = MP.perfil(n, p, a, L, 'bernoulli')
        for (c2, cp, tam, piv), (_, _, _, piv2) in zip(sal, sal2):
            if tuple(piv) != tuple(piv2):
                difs += 1
                print("  DIFIERE n=%d p=%d a=%d chi(2)=%s: masas %s padico %s" % (n, p, a, c2, piv, piv2), flush=True)
        for c, w in sorted(med.items()):
            W = d - npr + sum(tam for _, _, tam, piv in sal if (c - P) in piv)
            tot += 1
            ok += (W == w)
            if W != w:
                print("FALLA n=%d p=%d capa %d: %d vs %d" % (n, p, c, W, w), flush=True)
        print("n=%d p=%d a=%d: %d orbitas, %.2f s, misma cantidad de orbitas que padico: %s" % (
            n, p, a, len(sal), dt, len(sal) == len(sal2)), flush=True)
    print("CONTROL masas: %d de %d; orbitas con pivotes distintos de mahler_padico: %d" % (ok, tot, difs), flush=True)


def _masas():
    """(a): (9) frente a (8) en casos pequenos, con el factor chi(M): w = -chi(M) v (9) y w = -v' (8)."""
    casos = fallos = 0
    for n in (3, 5, 7, 9, 11, 13, 15, 21, 23, 25, 27, 29):
        for p in (5, 7, 11, 13):
            if n % p == 0 or Integers(n).unit_group_order() % p == 0:
                continue
            obs, _ = orbitas(n, p)
            for chi, orb, m, Km, chim in obs:
                ex = exponentes(chim, Km, n)
                g = MP.primo_sobre_p(m, p)
                for M in (p, p * p, p ** 3):
                    v = masas_cerradas(ex, m, n, p, M, g)
                    vd = masas_directas(ex, m, n, p, M, g)
                    lhs = (v @ mult_matriz(g, ex[M % n], m).T) % p  # chi(M) v
                    casos += 1
                    if not np.array_equal(lhs, vd):
                        fallos += 1
                        print("  (9) != (8): n=%d p=%d M=%d" % (n, p, M), flush=True)
    assert casos > 0
    print("MASAS (9) vs (8): %d casos, %d fallos" % (casos, fallos), flush=True)


def cuadratico(n, p, a, L):
    chi = kronecker_character(-n)
    assert chi(-1) == -1 and chi.modulus() == n
    vals = chi.values()
    ex = [-1 if v == 0 else (0 if v == 1 else 1) for v in vals]
    B = matriz_B(ex, ex[2], 2, n, p, a, L)
    return B, B.pivot_rows(), chi(2), chi(p)


def _lambda_alta():
    n, p = 75619, 5
    for a, L in ((2, 48), (3, 164)):
        t0 = time.time()
        B, piv, c2, cp = cuadratico(n, p, a, L)
        aus = [J for J in range(0, L + 1, 2) if J not in piv]
        print("piso %d (J <= %d, %.1f s): chi(2) = %s, chi(p) = %s; filas pares ausentes: %s" % (
            a, L, time.time() - t0, c2, cp, aus), flush=True)
        print("   impares ausentes: %s" % [J for J in range(1, L + 1, 2) if J not in piv], flush=True)
        print("   rango filas <= 163: %s; <= 164: %s" % (B[:164].rank() if L >= 163 else '-',
                                                          B[:165].rank() if L >= 164 else '-'), flush=True)
        # certificado de dependencia: B28 + 3B26 + 4B22 + 4B18 + 3B10 + 2B6 = 0 sobre F_5, en la coordenada Y = 1+T (filas de B)
        cert = B[28] + 3 * B[26] + 4 * B[22] + 4 * B[18] + 3 * B[10] + 2 * B[6]
        print("   certificado (4) en piso %d: %s" % (a, 'CERO' if cert.is_zero() else 'NO cero'), flush=True)


if __name__ == '__main__':
    modo = sys.argv[1] if len(sys.argv) > 1 else 'control'
    if modo == 'control':
        _control()
    elif modo == 'masas':
        _masas()
    elif modo == 'lambda_alta':
        _lambda_alta()
