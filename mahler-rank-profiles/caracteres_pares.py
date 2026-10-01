# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>
#
# QUE: control de la observacion sobre caracteres PARES del paper III (seccion de la prueba del Teorema B), en casos
#      donde una rama impar b tiene orden r_{chi,b} >= 1. chi = (D/.) cuadratico par de conductor n = D (D = 1 mod 4).
#      Para (n, p, a, b):
#        (1) medida mu_chi (masas por la ec. (masas), mu_chi = -2(mu_F - chi(2)[2]_* mu_F)); paridad de mu_chi;
#            r_b = ord_T f_b para todas las ramas (nivel p^{a+1}, truncado en p^a), y para b impar la descomposicion
#            r_b = [chi(2) omega(2)^b = 1] p^{v_p(theta(2))} + ord_T g_b, g_b = int_{Z_p^x} omega^b (1+T)^theta dmu_F;
#            g_b(0) = 0 frente a B_{b+1,chi} = 0 mod p (b+1 <= p-2);
#        (2) filas: B_{J,x} = int C(xy, J) dmu_chi, J, x < p^a, filas dependientes por eliminacion mod p, frente al
#            Teorema A con esos r_b (truncados en p^{a-1});
#        (3) capas: A = F_p-algebra generada por los e_k(eta_1..eta_m), eta_j = C_j(x), 2m+1 = p^a n, dentro de
#            F_p[x]/psi_n(x)^L, L = 2p^a - p^{a-1} (las capas de grado < L del orden no dependen de v si
#            phi(p^v) >= L, porque psi_{p^v n} = psi_n^{phi(p^v)} mod p); capas G_i en la base x^k psi_n^i;
#            componente chi: pi_chi(y) = sum_c chi(c) y(omega^c + omega^-c) w_c^i, w_c = psi_n'(.)(omega^c - omega^-c),
#            calculada en F_p[omega]/Phi_n (la anulacion no depende de la componente: Galois permuta las componentes y
#            multiplica pi_chi por chi(u));
#        (4) para 0 <= j < p^a - p^{a-1}: componente chi de G_{P+j} no nula  <=>  fila j de B nueva.
#      Control que debe fallar: la prediccion del Teorema A con r_b = 0 en la rama b elegida debe discrepar de las capas.
#      Python 3 + numpy, sin Sage. Memoria: matrices de lado D = d L <= ~800 (unos MB); aborta si la estimacion
#      supera 2 GB.
# POR QUE: la observacion afirma el Teorema B para chi par; los ejemplos con todas las r_b de ramas impares nulas no
#      distinguen la regla de "fila j nueva para todo j impar". Aqui la rama b tiene r_b >= 1: por el factor de Euler
#      en 2 (n = 17, p = 7, b = 3: chi(2) = 1, omega(2)^3 = 1) y por un numero de Bernoulli (n = 37, p = 5, b = 1:
#      B_{2,chi} = 20 = 0 mod 5).
# FECHA: 1-oct-2026
# Uso: python caracteres_pares.py > caracteres_pares_OUT.txt   (los dos casos: (17,7,2,3) y (37,5,2,1); ~4 min)
#      python caracteres_pares.py n p a b          (un caso; n = D primo con D = 1 mod 4)
import os
import sys
import time
from fractions import Fraction
from math import comb, gcd

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo import cyclo, phi, psi  # noqa: E402

LIMITE_BYTES = 2 * 1024 ** 3


# ------------------------------------------------------------------ aritmetica mod p
def kron(D, x):
    """simbolo de Kronecker (D/x), x >= 0."""
    if x == 0:
        return 1 if abs(D) == 1 else 0
    res = 1
    while x % 2 == 0:
        x //= 2
        if D % 2 == 0:
            return 0
        res *= 1 if D % 8 in (1, 7) else -1
    a, n = D % x, x
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            res = -res
        a %= n
    return res if n == 1 else 0


def mm(A, B, p):
    """A @ B mod p, exacto en coma flotante (se comprueba la cota)."""
    assert A.shape[1] * (p - 1) ** 2 < 2 ** 52
    return np.rint(A.astype(np.float64) @ B.astype(np.float64)).astype(np.int64) % p


def matmul_anillo(g, f, p):
    """matriz D x D cuya fila k es x^k g mod (f, p); f monico de grado D."""
    D = len(f) - 1
    M = np.zeros((D, D), dtype=np.int64)
    r = np.array(g, dtype=np.int64) % p
    fl = np.array(f[:D], dtype=np.int64)
    for k in range(D):
        M[k] = r
        top = r[D - 1]
        r = np.concatenate(([0], r[:D - 1]))
        if top:
            r = (r - top * fl) % p
    return M


def por_x(v, f, p):
    D = len(f) - 1
    top = v[D - 1]
    r = np.concatenate(([0], v[:D - 1]))
    return (r - top * np.array(f[:D], dtype=np.int64)) % p


class Escalon:
    """base escalonada reducida mod p; pivote = primera columna no nula de cada fila."""

    def __init__(self, D, p):
        self.p, self.D = p, D
        self.B = np.zeros((0, D), dtype=np.int64)
        self.piv = []

    def reduce(self, X):
        X = np.atleast_2d(X) % self.p
        if not self.piv:
            return X
        return (X - mm(X[:, self.piv], self.B, self.p)) % self.p

    def add(self, X):
        p = self.p
        R = self.reduce(X)
        R = R[R.any(axis=1)]
        nuevas = []
        while len(R):
            r = R[0]
            c = int(np.flatnonzero(r)[0])
            r = r * pow(int(r[c]), -1, p) % p
            R = R[1:]
            if len(R):
                R = (R - np.outer(R[:, c], r)) % p
                R = R[R.any(axis=1)]
            if len(self.B):
                self.B = (self.B - np.outer(self.B[:, c], r)) % p
            self.B = np.vstack([self.B, r])
            self.piv.append(c)
            nuevas.append(r)
        return np.array(nuevas, dtype=np.int64).reshape(len(nuevas), self.D)


# ------------------------------------------------------------------ (1) medida y ramas
def masas_F(ch, n, p, M):
    r = np.arange(n, dtype=np.int64)
    mF = np.array([int((r * ch[(a + M * r) % n]).sum() % p) for a in range(M)], dtype=np.int64)
    return (-mF * pow(n, -1, p)) % p


def masas_chi(ch, n, p, M):
    mF = masas_F(ch, n, p, M)
    idx = (np.arange(M) * pow(2, -1, M)) % M
    return (-2 * (mF - int(ch[2 % n]) * mF[idx])) % p


def theta(p, N):
    mod = p ** N
    dl, c = {}, 1
    for e in range(p ** (N - 1)):
        dl[c] = e
        c = c * (1 + p) % mod
    out = {}
    for y in range(mod):
        if y % p:
            om = pow(y, p ** (N - 1), mod)
            out[y] = dl[y * pow(om, -1, mod) % mod]
    return out


def orden_T(m, p, N, b, th):
    """ord_T de int_{Z_p^x} omega^b (1+T)^theta dm, truncado: < p^{N-1}, o None (>= p^{N-1})."""
    for i in range(p ** (N - 1)):
        s = sum(int(m[y]) * pow(y % p, b, p) * comb(t, i) for y, t in th.items() if m[y]) % p
        if s:
            return i
    return None


def bernoulli_chi(k, ch, n):
    """B_{k,chi} = n^{k-1} sum_{a=1}^{n} chi(a) B_k(a/n)."""
    Bl = [Fraction(0)] * (k + 1)
    A = [Fraction(0)] * (k + 1)
    for mm_ in range(k + 1):
        A[mm_] = Fraction(1, mm_ + 1)
        for j in range(mm_, 0, -1):
            A[j - 1] = j * (A[j - 1] - A[j])
        Bl[mm_] = A[0]
    Bl[1] = Fraction(-1, 2)

    def Bpol(x):
        return sum(comb(k, i) * Bl[i] * x ** (k - i) for i in range(k + 1))
    return n ** (k - 1) * sum(int(ch[a % n]) * Bpol(Fraction(a, n)) for a in range(1, n + 1))


def mod_p(q, p):
    return q.numerator * pow(q.denominator, -1, p) % p


# ------------------------------------------------------------------ (2) filas
def matriz_B(m, p, P):
    Bn = np.zeros((P, P), dtype=np.int64)
    Bn[0, 0] = 1
    for z in range(1, P):
        Bn[z, 0] = 1
        Bn[z, 1:] = (Bn[z - 1, 1:] + Bn[z - 1, :-1]) % p
    y = np.arange(P, dtype=np.int64)
    B = np.zeros((P, P), dtype=np.int64)
    for x in range(P):
        mass = np.zeros(P, dtype=np.int64)
        np.add.at(mass, (x * y) % P, m)
        B[:, x] = (mass % p) @ Bn % p
    return B


def filas_nuevas(B, p):
    E = Escalon(B.shape[1], p)
    nuevas = []
    for J in range(B.shape[0]):
        if len(E.add(B[J:J + 1])):
            nuevas.append(J)
    return nuevas


def u_de(J, p):
    x = J + 1
    while x % p == 0:
        x //= p
    return x


def teorema_A(p, a, r, eps):
    """filas dependientes segun el Teorema A; r[b] entero o None (infinito o >= p^{a-1})."""
    dep = [] if eps else [0]
    for J in range(1, p ** a):
        rb = r[J % (p - 1)]
        if rb is None or u_de(J, p) < p * rb:
            dep.append(J)
    return dep


# ------------------------------------------------------------------ (3) capas del orden
def capas(n, p, a):
    P = p ** a
    L = 2 * P - p ** (a - 1)
    ps = [c % p for c in psi(n)]
    d = len(ps) - 1
    D = d * L
    m = (P * n - 1) // 2
    est = 8 * (6 * D * D + (m + 2) * D)
    print("   capas: d=%d L=%d D=%d m=%d; memoria estimada %.0f MB" % (d, L, D, m, est / 2 ** 20), flush=True)
    if est > LIMITE_BYTES:
        sys.exit("ABORTA: memoria estimada > 2 GB")
    f = [1]
    for _ in range(L):
        f = list(np.convolve(f, ps) % p)
    f = [int(c) for c in f]
    assert len(f) == D + 1 and f[-1] == 1
    # eta_j = C_j(x) mod f (Chebyshev) y e_k(eta_1..eta_m) por la recurrencia de prod (X - eta_j)
    x = np.zeros(D, dtype=np.int64)
    x[1] = 1
    eta_ant = np.zeros(D, dtype=np.int64)
    eta_ant[0] = 2 % p
    eta = x.copy()
    C = np.zeros((1, D), dtype=np.int64)
    C[0, 0] = 1
    for j in range(1, m + 1):
        Me = matmul_anillo(eta, f, p)
        Cn = np.zeros((len(C) + 1, D), dtype=np.int64)
        Cn[1:] += C
        Cn[:-1] -= mm(C, Me, p)
        C = Cn % p
        eta, eta_ant = (por_x(eta, f, p) - eta_ant) % p, eta
    gens = [((-1) ** k * C[m - k]) % p for k in range(1, m + 1)]
    # cierre: generadores de algebra elegidos de los e_k; semi-ingenuo
    alg = Escalon(D, p)
    uno = np.zeros(D, dtype=np.int64)
    uno[0] = 1
    alg.add(uno)
    T = []
    for g in gens:
        if not alg.reduce(g).any():
            continue
        Mt = matmul_anillo(g, f, p)
        viejo = alg.B.copy()
        T.append(Mt)
        frente = alg.add(g)
        frente = np.vstack([frente, alg.add(mm(viejo, Mt, p))])
        while len(frente):
            nuevo = []
            for M in T:
                nuevo.append(alg.add(mm(frente, M, p)))
            frente = np.vstack(nuevo)
    dimA = len(alg.B)
    # base x^k psi^i (fila i*d+k) -> monomios; matriz unitriangular inferior; inversa por duplicacion
    Cm = np.zeros((D, D), dtype=np.int64)
    pw = np.array([1], dtype=np.int64)
    for i in range(L):
        for k in range(d):
            Cm[i * d + k, k:k + len(pw)] = pw
        pw = np.convolve(pw, ps) % p
    N = (np.eye(D, dtype=np.int64) - Cm) % p
    Ci = (np.eye(D, dtype=np.int64) + N) % p
    Nk = N
    for _ in range(int(np.ceil(np.log2(D))) + 1):
        Nk = mm(Nk, Nk, p)
        if not Nk.any():
            break
        Ci = mm(Ci, (np.eye(D, dtype=np.int64) + Nk) % p, p)
    assert (mm(Cm, Ci, p) == np.eye(D, dtype=np.int64)).all()
    dig = mm(alg.B, Ci, p)
    E = Escalon(D, p)
    E.add(dig)
    G = [[] for _ in range(L)]
    for r, c in zip(E.B, E.piv):
        i = c // d
        G[i].append(r[i * d:(i + 1) * d])
    G = [np.array(g, dtype=np.int64).reshape(len(g), d) for g in G]
    return dict(P=P, L=L, d=d, D=D, m=m, G=G, dimA=dimA, nT=len(T), ps=ps)


def componentes_chi(n, p, ch, cap):
    """grados i con componente chi de G_i no nula."""
    cf = [c % p for c in cyclo(n)]
    e = len(cf) - 1
    unid = [c for c in range(1, n) if gcd(c, n) == 1]
    ypow = [np.zeros(e, dtype=np.int64)]
    ypow[0][0] = 1
    for _ in range(n):
        ypow.append(por_x(ypow[-1], cf, p))
    d, L, ps = cap['d'], cap['L'], cap['ps']
    dps = [(k * ps[k]) % p for k in range(1, d + 1)]
    out = []
    datos = {}
    for c in unid:
        xc = (ypow[c] + ypow[n - c]) % p
        Mx = matmul_anillo(xc, cf, p)
        Ec = np.zeros((d, e), dtype=np.int64)
        v = ypow[0].copy()
        for k in range(d):
            Ec[k] = v
            v = mm(v[None], Mx, p)[0]
        dval = (np.array(dps, dtype=np.int64) @ Ec) % p
        w = mm(dval[None], matmul_anillo((ypow[c] - ypow[n - c]) % p, cf, p), p)[0]
        datos[c] = (Ec, matmul_anillo(w, cf, p))
    Wpow = {c: np.eye(e, dtype=np.int64) for c in unid}
    for i in range(L):
        if i > 0:
            Wpow = {c: mm(Wpow[c], datos[c][1], p) for c in unid}
        Gi = cap['G'][i]
        if not len(Gi):
            continue
        tot = np.zeros((len(Gi), e), dtype=np.int64)
        for c in unid:
            tot = (tot + int(ch[c]) * mm(mm(Gi, datos[c][0], p), Wpow[c], p)) % p
        if tot.any():
            out.append(i)
    return out


# ------------------------------------------------------------------ un caso
def caso(n, p, a, b_sel):
    t0 = time.time()
    P = p ** a
    ch = np.array([kron(n, x) for x in range(n)], dtype=np.int64)
    print("=" * 100)
    print("(n, p, a) = (%d, %d, %d), chi = (%d/.), rama impar elegida b = %d" % (n, p, a, n, b_sel))
    assert n % 4 == 1 and (n * phi(n)) % p and p >= 5
    print("   chi(-1) = %d, chi(2) = %d, chi(p) = %d, omega(2)^b = 2^b mod p = %d"
          % (int(ch[n - 1]), int(ch[2]), int(ch[p % n]), pow(2, b_sel, p)))
    # (1)
    Nn = a + 1
    mN = masas_chi(ch, n, p, p ** Nn)
    mFN = masas_F(ch, n, p, p ** Nn)
    impar = all((mN[(-y) % p ** Nn] + mN[y]) % p == 0 for y in range(p ** Nn))
    th = theta(p, Nn)
    v2 = 0 if th[2] % p else 1
    r = [orden_T(mN, p, Nn, b, th) for b in range(p - 1)]
    print("   (1) mu_chi impar (mu(-X) = -mu(X)): %s; mu_chi(Z_p) = %d mod p; v_p(theta(2)) %s 0"
          % (impar, int(mN.sum() % p), '=' if v2 == 0 else '>'))
    print("       r_b (nivel p^%d, None = >= %d): %s" % (Nn, p ** (Nn - 1), r))
    ok_R = True
    for b in range(1, p - 1, 2):
        og = orden_T(mFN, p, Nn, b, th)
        eu = 1 if (int(ch[2]) * pow(2, b, p)) % p == 1 else 0
        pred = None if (og is None or r[b] is None) else eu + og
        ok_R &= (pred == r[b])
        linea = "       b=%d: [chi(2)omega(2)^b=1] = %d, ord g_b = %s, suma = %s, r_b = %s" % (b, eu, og, pred, r[b])
        if b + 1 <= p - 2:
            Bk = bernoulli_chi(b + 1, ch, n)
            linea += "; B_{%d,chi} = %s = %d mod p, g_b(0) = 0: %s" % (b + 1, Bk, mod_p(Bk, p), og is not None and og >= 1)
            ok_R &= ((mod_p(Bk, p) == 0) == (og is not None and og >= 1))
        print(linea)
    print("       ramas impares: r_b = factor de Euler + ord g_b (y g_b(0) = 0 <=> B_{b+1,chi} = 0): %s" % ok_R)
    print("       rama elegida b = %d: r_b = %s %s" % (b_sel, r[b_sel], '>= 1' if (r[b_sel] is None or r[b_sel] >= 1) else '= 0 (el caso NO sirve)'))
    # (2)
    m = masas_chi(ch, n, p, P)
    Bm = matriz_B(m, p, P)
    nuevas = filas_nuevas(Bm, p)
    dep = [J for J in range(P) if J not in nuevas]
    eps = int(m.sum() % p) != 0
    rt = [None if (x is None or x >= p ** (a - 1)) else x for x in r]
    predA = teorema_A(p, a, rt, eps)
    jw = P - p ** (a - 1)
    print("   (2) filas de B (%d x %d): rango %d; perfil medido = Teorema A: %s" % (P, P, len(nuevas), dep == predA))
    print("       filas dependientes j < %d con j impar: %s" % (jw, [j for j in dep if j < jw and j % 2]))
    # (3)
    cap = capas(n, p, a)
    print("       dim A mod psi^L = %d, generadores de algebra usados: %d, %.1f s" % (cap['dimA'], cap['nT'], time.time() - t0))
    nz = componentes_chi(n, p, ch, cap)
    nzw = [i for i in nz if P <= i < P + jw]
    print("   (3) grados P+j, 0 <= j < %d, con componente chi no nula: %s" % (jw, nzw))
    print("       (ninguno de grado impar: %s)" % all(i % 2 == 0 for i in nz))
    # (4)
    linea, fallos = [], 0
    for j in range(jw):
        capa = (P + j) in nzw
        fila = j in nuevas
        fallos += capa != fila
        linea.append('.' if capa and fila else ('x' if not capa and not fila else '!'))
    print("   (4) capa vs fila (. ambas no nulas/nuevas, x ambas nulas/dependientes, ! discrepancia), j = 0..%d:" % (jw - 1))
    print("       %s" % ''.join(linea))
    print("       j impares con capa nula: %s" % [j for j in range(1, jw, 2) if (P + j) not in nzw])
    print("       discrepancias capa/fila: %d de %d" % (fallos, jw))
    r0 = list(rt)
    r0[b_sel] = 0
    pred0 = teorema_A(p, a, r0, eps)
    dis = sum(((P + j) in nzw) != (j not in pred0) for j in range(jw))
    print("   control (debe fallar): Teorema A con r_%d = 0 frente a las capas: %d discrepancias (en j = %s)"
          % (b_sel, dis, [j for j in range(jw) if ((P + j) in nzw) != (j not in pred0)]))
    ok = impar and ok_R and dep == predA and fallos == 0 and dis > 0 and (r[b_sel] is None or r[b_sel] >= 1)
    print("   VEREDICTO (n,p,a,b) = (%d,%d,%d,%d): %s  (%.1f s)" % (n, p, a, b_sel, 'ACUERDO' if ok else 'DESACUERDO', time.time() - t0),
          flush=True)
    return ok


if __name__ == '__main__':
    if len(sys.argv) > 4:
        casos = [tuple(int(x) for x in sys.argv[1:5])]
    else:
        casos = [(17, 7, 2, 3), (37, 5, 2, 1)]
    res = [caso(*c) for c in casos]
    print("=" * 100)
    print("RESULTADO: %d/%d casos en acuerdo" % (sum(res), len(res)))
