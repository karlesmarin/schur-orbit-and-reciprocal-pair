# -*- coding: utf-8 -*-
# inst_tipoA.sage -- el instrumento de capas de `inst.sage`, con el SEGMENTO como parametro.
#
# `inst.sage` calcula la imagen de A_m(q) = Z[e_1..e_m] en Z[xi]/(p, pi^K) = F_p[y]/(Phi_n^K)
# con generadores e_k(eta_1..eta_m), eta_j = y^j + y^{-j}: eso es TIPO C y esta cableado.
# Aqui lo unico que cambia es de donde salen los generadores: se pasa la lista de EXPONENTES del
# elemento y se toman los e_k de los autovalores y^{exp}.  Asi el mismo instrumento sirve para
#
#   tipo C:  exps = [+-1, ..., +-m]        sobre M = q
#   tipo A SU: exps = [(n+1-2j)/2]         sobre M = q     (n impar)
#              exps = [n+1-2j]             sobre M = 2q    (n par)
#   tipo A GL: exps = [0, 1, ..., n-1]     sobre M = q
#
# El algebra es identica: arbol de productos empaquetado (Kronecker) para el polinomio
# caracteristico, cierre multiplicativo del F_p-espacio, y rangos modulo Phi^k.
#
# Authors: Carles Marin, Claude (AI assistant).


class ImgSeg:
    def __init__(self, M, p, exps, r, K):
        v = valuation(M, p)
        n = M // p ** v
        self.M, self.p, self.exps, self.r, self.K = M, p, list(exps), r, K
        self.v, self.n = v, n
        self.E = euler_phi(p ** v) if v > 0 else 1
        assert 1 <= K <= self.E, "K debe ser <= E = phi(p^v) = %d" % self.E
        F = GF(p); self.F = F
        R = PolynomialRing(F, 'y'); self.R = R
        y = R.gen()
        self.phin = R(cyclotomic_polynomial(n))
        assert R(cyclotomic_polynomial(M)) == self.phin ** self.E, \
            "Phi_M no es Phi_n^E mod p: M=%d p=%d" % (M, p)
        Mod = self.phin ** K; self.Mod = Mod
        D = Mod.degree(); self.D = D
        facs = []
        for e in self.exps:
            facs.append([-(y ** (e % M)) % Mod, R(1)])
        coefs = self._tree(facs) if facs else [R(1)]
        nn = len(self.exps)
        assert len(coefs) == nn + 1 and coefs[nn] == 1
        # e_k = (-1)^k * coef de T^{nn-k};  se toman k = 1..r
        gens = [((-1) ** k * coefs[nn - k]) % Mod for k in range(1, r + 1)]
        self.gens = [g for g in gens if g != 0]
        self.B = self._closure(self.gens)

    def _mul(self, A, B):
        R, Mod, D = self.R, self.Mod, self.D
        if len(A) * len(B) <= 64:
            out = [R(0)] * (len(A) + len(B) - 1)
            for i, a in enumerate(A):
                if a:
                    for j, b in enumerate(B):
                        out[i + j] += a * b
            return [c % Mod for c in out]
        S = 2 * D - 1
        def pack(L):
            c = [0] * (len(L) * S)
            for i, a in enumerate(L):
                l = a.list()
                c[i * S:i * S + len(l)] = [int(x) for x in l]
            return R(c)
        P = pack(A) * pack(B)
        l = P.list()
        L = len(A) + len(B) - 1
        l = l + [0] * (L * S - len(l))
        return [R(l[i * S:(i + 1) * S]) % Mod for i in range(L)]

    def _tree(self, facs):
        while len(facs) > 1:
            nxt = []
            for i in range(0, len(facs) - 1, 2):
                nxt.append(self._mul(facs[i], facs[i + 1]))
            if len(facs) % 2:
                nxt.append(facs[-1])
            facs = nxt
        return facs[0]

    def vec(self, f, D=None):
        D = D or self.D
        l = f.list()
        return l + [0] * (D - len(l))

    def _closure(self, gens):
        F, R, Mod = self.F, self.R, self.Mod
        G = matrix(F, [self.vec(R(1))] + [self.vec(g) for g in gens]).echelon_form()
        r = G.rank()
        G = G.matrix_from_rows(range(r))
        gb = [R(list(row)) for row in G.rows()]
        B = G
        while True:
            bb = [R(list(row)) for row in B.rows()]
            prods = [self.vec((a * g) % Mod) for a in bb for g in gb]
            S = B.stack(matrix(F, prods)).echelon_form()
            r2 = S.rank()
            if r2 == B.nrows():
                return B
            B = S.matrix_from_rows(range(r2))

    def dims(self, K=None):
        K = K or self.K
        out = []
        for k in range(1, K + 1):
            Mk = self.phin ** k
            Dk = Mk.degree()
            rows = [self.vec(self.R(list(r)) % Mk, Dk) for r in self.B.rows()]
            out.append(matrix(self.F, rows).rank())
        return out

    def W(self):
        dB = self.dims()
        return [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, len(dB))]


def exps_su(n):
    if n % 2 == 1:
        return [(n + 1 - 2 * j) // 2 for j in range(1, n + 1)]
    return [(n + 1 - 2 * j) for j in range(1, n + 1)]


def exps_C(m):
    out = []
    for j in range(1, m + 1):
        out += [j, -j]
    return out
