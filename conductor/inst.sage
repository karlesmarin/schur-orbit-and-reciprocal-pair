# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# inst.sage -- fast layer instrument: image of A_m(q)=Z[e_1..e_m] in Z[xi]/(p, pi^K) = F_p[y]/(Phi_n(y)^K), K<=E.
# What: e_k(eta_1..eta_m), eta_j=y^j+y^-j, by a Kronecker-packed product tree (FLINT); ring closure; dims mod Phi_n^k.
# Run: load("inst.sage") from another sage script.
import sys, time

def split_q(q, p):
    v = valuation(q, p)
    return v, q // p**v

class Img:
    def __init__(self, q, p, m, K):
        v, n = split_q(q, p)
        self.q, self.p, self.m, self.K, self.v, self.n = q, p, m, K, v, n
        self.E = euler_phi(p**v)
        assert 1 <= K <= self.E, "K must be <= E"
        F = GF(p); self.F = F
        R = PolynomialRing(F, 'y'); self.R = R; y = R.gen()
        self.phin = R(cyclotomic_polynomial(n))
        assert R(cyclotomic_polynomial(q)) == self.phin**self.E
        Mod = self.phin**K; self.Mod = Mod; D = Mod.degree(); self.D = D
        yinv = y.inverse_mod(Mod)
        # eta_j mod Mod, via y^j (period of y in the ring divides the order, fine to just power)
        facs = []
        yj = R(1); yij = R(1)
        for j in range(1, m + 1):
            yj = (yj * y) % Mod; yij = (yij * yinv) % Mod
            facs.append([-(yj + yij) % Mod, R(1)])
        coefs = self._tree(facs) if facs else [R(1)]
        assert len(coefs) == m + 1 and coefs[m] == 1
        gens = [coefs[k] for k in range(m)]   # e_{m-k} up to sign
        self.B = self._closure(gens)

    def _mul(self, A, B):
        """product of two Z-polys with coefficients in R_{<D} (lists low->high), reduced mod Mod."""
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
        l = f.list(); return l + [0] * (D - len(l))

    def _closure(self, gens):
        F, R, Mod, D = self.F, self.R, self.Mod, self.D
        G = matrix(F, [self.vec(R(1))] + [self.vec(g) for g in gens]).echelon_form()
        r = G.rank(); G = G.matrix_from_rows(range(r))
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
            Mk = self.phin**k; Dk = Mk.degree()
            rows = [self.vec(self.R(list(r)) % Mk, Dk) for r in self.B.rows()]
            out.append(matrix(self.F, rows).rank())
        return out

    def W(self):
        dB = self.dims()
        return [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, len(dB))]

    def max_ideal_square_zero(self):
        """is m^2 = 0 in the image mod Phi^K ?  (m = kernel of reduction mod Phi)"""
        F, R = self.F, self.R
        P1 = matrix(F, [self.vec(R(list(r)) % self.phin, self.phin.degree()) for r in self.B.rows()])
        ker = P1.left_kernel().basis_matrix()
        mm = [R(list(row)) for row in (ker * self.B).rows()]
        for i in range(len(mm)):
            for j in range(i, len(mm)):
                if (mm[i] * mm[j]) % self.Mod != 0:
                    return False, len(mm)
        return True, len(mm)

def moment_rank(n, mprime, p):
    S = [0] * n
    for j in range(1, mprime + 1):
        S[j % n] += j
        S[(-j) % n] -= j
    units = [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]
    rows = []
    for u in units:
        ui = inverse_mod(u, n)
        rows.append([S[(ui * c) % n] for c in range(1, n)])
    return matrix(GF(p), rows).rank(), len(units)
