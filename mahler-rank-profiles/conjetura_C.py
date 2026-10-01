# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: falsacion de la conjetura (C): para 2m+1 = pn, modulo Q = p^v n, n >= 3 impar, p >= 5, p no divide n phi(n),
#      y M_n = [s_n(u^-1 c)] de rango lleno phi(n)/2 mod p, el conductor de A = Z[e_k(eta_1..eta_m)] en p cumple
#      e <= 2p (r^{2p} en A_p).  Es la familia del descenso (Teorema B) con a = 1, M = 1, q' = n.
#      Instrumento (Sage, GF(p)): las capas de A_p en grado i < L solo ven la imagen de A en
#        O/(p, r^L) = F_p[x]/(psi_n(x)^L)       (valido si L <= E = phi(p^v): psi_Q = psi_n^E mod p, se comprueba),
#      y los e_k son polinomios en theta independientes de Q, asi que el resultado vale para TODO v con E >= L.
#      Se calcula el algebra generada por los e_k en ese cociente (cierre B <- B + B*B_1, podando valoraciones >= L),
#      en la base x^k psi^i (digitos de la expansion en base psi), con pivote = menor valoracion; de ahi
#        G_i = gr_i(A) subset k = F_p[x]/psi_n,   W_i = dim G_i,   Q_i = G_i / sum_{0<j<i} G_j G_{i-j}  (generadores nuevos),
#      y el exponente e = min{i : G_j = k para todo i <= j < 2i}.  CERTIFICADO DE COLA: si G_j = k en [e, 2e), todo
#      grado >= e es suma de grados de [e, 2e) y el producto de capas llenas es lleno (gr(O) = k[u]), luego r^e en A_p;
#      y G_{e-1} != k da la agudeza.  Si 2e > L no hay certificado y se dice.
#      Controles: (i) calibracion contra nucleo.local_completo (reticulo, sin truncar) en casos pequenos;
#      (ii) near-miss: casos con M_n de rango deficiente (fuera de la hipotesis) se miden igual, para ver si el
#      instrumento puede dar e > 2p.
# POR QUE: SEGUIR_AQUI_2026-09-29 punto 3: falsar antes de probar; lote estratificado con p = 17 reservado.
# FECHA: 29-sep-2026
# Uso (Git Bash):
#   MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/note_conductor/ext_clasicos:/work" \
#       --entrypoint bash -w /work sage-normaliz:local -c "sage -python conjetura_C.py <modo> > conjetura_C_<modo>_OUT.txt"
#   modo = calibra | lote | p17
import sys
import time
from math import gcd

from sage.all import GF, PolynomialRing, matrix, vector

sys.path.insert(0, '/work')
from nucleo import cheb_polys, matriz_momento, phi, psi, rango_mod, sierra, sierra_desplazada  # noqa: E402


def es_primo(n):
    return n > 1 and all(n % q for q in range(2, int(n ** 0.5) + 1))


class Caso:
    def __init__(self, n, p, L, a=1):
        self.n, self.p, self.L = n, p, L
        F = GF(p)
        self.F = F
        R = PolynomialRing(F, 'x')
        self.R = R
        x = R.gen()
        self.x = x
        self.psi = R([int(c) for c in psi(n)])
        self.d = self.psi.degree()
        self.f = self.psi ** L
        self.D = self.d * L
        self.m = (p ** a * n - 1) // 2  # a = 1: la familia de (C); a = 2: autoprueba

    def digitos(self, g):
        """coordenadas de g (grado < D) en la base x^k psi^i, i < L, k < d."""
        out = []
        for _ in range(self.L):
            g, r = g.quo_rem(self.psi)
            c = r.list()
            out.extend(c + [0] * (self.d - len(c)))
        assert g == 0
        return out

    def poli(self, row):
        g = self.R(0)
        P = self.R(1)
        d = self.d
        for i in range(self.L):
            blk = row[i * d:(i + 1) * d]
            if any(blk):
                g += self.R(list(blk)) * P
            P *= self.psi
        return g

    def generadores(self):
        C = cheb_polys(self.m)
        etas = [self.R([int(c) for c in C[a]]) % self.f for a in range(1, self.m + 1)]
        coef = [self.R(1)]
        for eta in etas:
            new = [self.R(0)] * (len(coef) + 1)
            for i, c in enumerate(coef):
                new[i + 1] += c
                new[i] -= (eta * c) % self.f
            coef = new
        m = self.m
        return [((-1) ** k * coef[m - k]) % self.f for k in range(1, m + 1)]

    def val(self, row):
        for j, c in enumerate(row):
            if c:
                return j // self.d
        return None

    def echelon(self, filas):
        M = matrix(self.F, filas).echelon_form()
        return [r for r in M.rows() if not r.is_zero()]

    def cierre(self):
        gens = self.generadores()
        B1 = self.echelon([self.digitos(self.R(1))] + [self.digitos(g) for g in gens])
        B1c = [(r, self.val(r), self.poli(r)) for r in B1 if self.val(r) > 0]
        B = B1
        while True:
            Bc = [(r, self.val(r), self.poli(r)) for r in B if self.val(r) > 0]
            nuevas = []
            for r, vr, gr in Bc:
                for s, vs, gs in B1c:
                    if vr + vs < self.L:
                        nuevas.append(self.digitos((gr * gs) % self.f))
            B2 = self.echelon([list(r) for r in B] + nuevas) if nuevas else B
            if len(B2) == len(B):
                return B
            B = B2

    def capas(self, B):
        d, L = self.d, self.L
        G = [[] for _ in range(L)]
        for r in B:
            i = self.val(r)
            G[i].append(list(r[i * d:(i + 1) * d]))
        return G

    def producto_k(self, U, V):
        """span de U*V en k = F_p[x]/psi."""
        out = []
        for u in U:
            for w in V:
                c = ((self.R(u) * self.R(w)) % self.psi).list()
                out.append(c + [0] * (self.d - len(c)))
        return out

    def rango(self, filas):
        return matrix(self.F, filas).rank() if filas else 0


def medir(n, p, L=None, a=1):
    if L is None:
        L = 4 * p + 2
    t0 = time.time()
    c = Caso(n, p, L, a)
    B = c.cierre()
    G = c.capas(B)
    d = c.d
    W = [len(g) for g in G]
    Q = []
    for i in range(L):
        prods = []
        for j in range(1, i):
            if G[j] and G[i - j]:
                prods += c.producto_k(G[j], G[i - j])
        rp = c.rango(prods)
        Q.append(W[i] - rp if W[i] else 0)
        assert rp <= W[i] or i == 0, (n, p, i, rp, W[i])  # los productos caen dentro de G_i
    e = None
    for i in range(1, L):
        if 2 * i > L:
            break
        if all(W[j] == d for j in range(i, 2 * i)) and W[i - 1] < d:
            e = i
            break
    rk = rango_mod(matriz_momento(sierra(n), n), p)
    rs = rango_mod(matriz_momento(sierra_desplazada(n), n), p)  # residuo con signo (hipotesis del descenso)
    return dict(n=n, p=p, d=d, L=L, W=W, Q=Q, e=e, rango=rk, lleno=(rk == d), rango_signo=rs, t=time.time() - t0, dimA=len(B))


def fila(r):
    p, e = r['p'], r['e']
    Wtxt = ','.join('d' if w == r['d'] else str(w) for w in r['W'][:min(r['L'], 2 * p + 4)])  # d = capa llena
    Qnz = {i: q for i, q in enumerate(r['Q']) if q and i > 0}
    ver = ('SIN CERTIFICADO (2e > L)' if e is None else ('e<=2p' if e <= 2 * p else 'e>2p  <<< FALLA'))
    return (f"n={r['n']:3d} p={p:2d} d={r['d']:2d} rango pura={r['rango']:2d} signo={r['rango_signo']:2d}{'' if r['lleno'] else ' DEFICIENTE'} "
            f"e={e} 2p={2 * p} {ver} | W[0..]={Wtxt} | Q_i(i>0)={Qnz} ({r['t']:.1f}s)")


def casos(ps, nmax):
    out = []
    for p in ps:
        for n in range(3, nmax + 1, 2):
            if (n * phi(n)) % p:
                out.append((n, p))
    return out


if __name__ == '__main__':
    modo = sys.argv[1] if len(sys.argv) > 1 else 'calibra'
    if modo == 'calibra':
        # psi_Q = psi_n^E mod p, en los Q usados
        for n, p, v in [(3, 5, 2), (7, 5, 2), (3, 7, 2), (5, 7, 2)]:
            R = PolynomialRing(GF(p), 'x')
            Q = p ** v * n
            E = phi(p ** v)
            assert R([int(c) for c in psi(Q)]) == R([int(c) for c in psi(n)]) ** E, (n, p, v)
            print(f"psi_{Q} = psi_{n}^{E} mod {p}: ok")
        from nucleo import local_completo, segmento_lleno
        for n, p, v in [(3, 5, 2), (7, 5, 2), (3, 7, 2), (5, 7, 2), (9, 5, 2)]:
            Q = p ** v * n
            E = phi(p ** v)
            L = min(E, 4 * p + 2)
            r = medir(n, p, L)
            a, e_ref, Wref = local_completo(Q, segmento_lleno((p * n - 1) // 2), p)
            k = min(len(Wref), L)
            print(fila(r))
            print(f"    nucleo.local_completo Q={Q}: a={a} e={e_ref} W={Wref[:k]}")
            print(f"    W coincide en i < {k}: {r['W'][:k] == Wref[:k]}   e coincide: {r['e'] == e_ref}"
                  f"   a = sum(d - W_i, i<e) coincide: {a == sum(r['d'] - w for w in r['W'][:e_ref])}")
            sys.stdout.flush()
    elif modo == 'lote':
        todos = casos([5, 7, 11, 13], 45)
        print(f"{len(todos)} casos (p en 5,7,11,13; n impar 3..45; p no divide n phi(n))")
        for n, p in todos:
            print(fila(medir(n, p)))
            sys.stdout.flush()
    elif modo == 'control':
        # near-miss: p | phi(n) (fuera de la hipotesis de (C)), p no divide n; puede el instrumento dar e > 2p?
        todos = [(n, p) for p in (5, 7, 11) for n in range(3, 70, 2) if n % p and phi(n) % p == 0]
        print(f"{len(todos)} casos de control (p | phi(n), p no divide n; n impar < 70, p en 5,7,11)")
        for n, p in todos:
            print(fila(medir(n, p)))
            sys.stdout.flush()
    elif modo == 'autoprueba':
        # el instrumento debe poder dar e > 2p: familia con a = 2 (2m+1 = p^2 n), escala natural 2P = 2p^2.
        # Validez: L <= E = phi(p^v) con v grande (las capas < L no dependen de v).
        for n, p in [(3, 5), (7, 5), (3, 7)]:
            L = 4 * p * p + 2
            r = medir(n, p, L, a=2)
            print(f"a=2 (2m+1 = {p * p * n}), L={L}: " + fila(r))
            sys.stdout.flush()
    elif modo == 'a2_extra':
        # dos casos a = 2: (11,5) (predicho e = 48, W_25 = 5, W_47 = 4) y el caso de control (31,5)
        # (5 | phi(31), r_signo = 12 < d = 15: la ley e <= 2P con igualdad sii r < d fuerza e = 50).
        # L = 100 = phi(5^3): valido para v >= 3; el certificado necesita 2e <= L.
        for n, p in [(11, 5), (31, 5)]:
            r = medir(n, p, 100, a=2)
            print(f"a=2 (2m+1 = {p * p * n}), L=100: " + fila(r))
            print(f"    W completo = {r['W']}")
            sys.stdout.flush()
    elif modo == 'p17':
        todos = casos([17], 45)
        print(f"{len(todos)} casos (p = 17, bloque reservado)")
        for n, p in todos:
            print(fila(medir(n, p)))
            sys.stdout.flush()
