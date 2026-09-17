# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# capas.sage -- layers W_i, kappa and F_p-type of A_m(q) at p, by linear algebra in O/(p, r^K) = F_p[X]/(Psi_{q'}^K), K<=E.
# What: independent Sage re-implementation of the layer instrument (X = alpha = xi+xi^-1, Dickson recursion for eta_j).
# Run: load("capas.sage") from another sage script.

def val_p(x, p):
    x = abs(Integer(x))
    return x.valuation(p) if x else Infinity

def split_q(q, p):
    v = valuation(q, p)
    return v, q // p**v

def Psi_ZZ(n):
    """minimal polynomial of zeta_n + zeta_n^-1 over ZZ (n>=3)."""
    R = PolynomialRing(ZZ, 'X'); X = R.gen()
    Phi = cyclotomic_polynomial(n)
    c = Phi.list(); d = Phi.degree() // 2
    V = [R(2), X]
    for k in range(2, d + 1):
        V.append(X * V[-1] - V[-2])
    return R(c[d]) + sum(c[d + k] * V[k] for k in range(1, d + 1))

class Capas:
    def __init__(self, q, p, m, K):
        self.q, self.p, self.m, self.K = q, p, m, K
        v, n = split_q(q, p)
        self.v, self.n = v, n
        self.E = euler_phi(p**v)
        assert K <= self.E, "K must be <= E"
        self.d = euler_phi(n) // 2
        F = GF(p); self.F = F
        R = PolynomialRing(F, 'X'); self.R = R; X = R.gen(); self.X = X
        self.psi = R(Psi_ZZ(n))
        # sanity: Psi_q == Psi_n^E mod p
        assert R(Psi_ZZ(q)) == self.psi**self.E
        Mod = self.psi**K; self.Mod = Mod
        D = Mod.degree(); self.D = D
        eta_prev, eta = R(2), X
        coef = [R(1)]
        for j in range(1, m + 1):
            new = [R(0)] * (len(coef) + 1)
            for i, c in enumerate(coef):
                new[i + 1] += c
                new[i] = (new[i] - eta * c) % Mod
            coef = new
            eta_prev, eta = eta, (X * eta - eta_prev) % Mod
        gens = coef[:-1]
        V = VectorSpace(F, D)
        def vec(f):
            l = f.list(); return vector(F, l + [0] * (D - len(l)))
        self.vec = vec
        G = matrix(F, [vec(R(1))] + [vec(g) for g in gens]).echelon_form()
        G = G.matrix_from_rows([i for i in range(G.nrows()) if G.row(i) != 0])
        gb = [R(list(r)) for r in G.rows()]
        mats = []
        for g in gb:
            rows = []; h = g
            for i in range(D):
                rows.append(vec(h)); h = (X * h) % Mod
            mats.append(matrix(F, rows))
        B = matrix(F, [vec(R(1))])
        while True:
            S = B.stack(G)
            for Mg in mats:
                S = S.stack(B * Mg)
            S = S.echelon_form()
            r = S.rank()
            S = S.matrix_from_rows(range(r))
            if r == B.nrows():
                break
            B = S
        self.B = B  # basis of image of A in F_p[X]/Psi^K
        self.mats_basis = None

    def poly(self, row):
        return self.R(list(row))

    def dims(self):
        out = []
        for k in range(1, self.K + 1):
            Mk = self.psi**k
            rows = [self.vec_k(self.poly(r) % Mk, k) for r in self.B.rows()]
            out.append(matrix(self.F, rows).rank())
        return out

    def vec_k(self, f, k):
        Dk = k * self.d
        l = f.list(); return vector(self.F, l + [0] * (Dk - len(l)))

    def W(self):
        dB = self.dims()
        return [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, len(dB))]

    def image_mod(self, k):
        """basis (rows) of the image of A in F_p[X]/Psi^k, as matrix over F of size r x kd."""
        Mk = self.psi**k
        M = matrix(self.F, [self.vec_k(self.poly(r) % Mk, k) for r in self.B.rows()]).echelon_form()
        return M.matrix_from_rows(range(M.rank()))

    def mult_matrix(self, f, k):
        Mk = self.psi**k; X = self.X
        rows = []; h = f % Mk
        for i in range(k * self.d):
            rows.append(self.vec_k(h, k)); h = (X * h) % Mk
        return matrix(self.F, rows)

    def kappa(self):
        """returns (W1, kappa): V = gr_1(A) in Rbar = F_p[X]/psi; kappa = dim {a : aV in V}."""
        assert self.K >= 2 and self.E >= 2
        F, R, psi, d = self.F, self.R, self.psi, self.d
        B2 = self.image_mod(2)
        # kernel of projection B2 -> F_p[X]/psi
        P = matrix(F, [self.vec_k(self.poly(r) % psi, 1) for r in B2.rows()])
        ker = P.left_kernel().basis_matrix()
        U = ker * B2  # elements of B2 in (psi)/(psi^2)
        Vl = []
        for u in U.rows():
            f = self.poly(u)
            qq, rr = f.quo_rem(psi)
            assert rr == 0
            Vl.append(self.vec_k(qq % psi, 1))
        if not Vl:
            return 0, d
        Vm = matrix(F, Vl).echelon_form(); W1 = Vm.rank(); Vm = Vm.matrix_from_rows(range(W1))
        if W1 == 0:
            return 0, d
        Z = Vm.right_kernel().basis_matrix().transpose()  # Vm*Z = 0
        blocks = []
        for vrow in Vm.rows():
            fv = self.poly(vrow)
            # rows: a = X^k -> (X^k * v mod psi) * Z
            Mk = matrix(F, [self.vec_k((self.X**k * fv) % psi, 1) for k in range(d)])
            blocks.append(Mk * Z if Z.ncols() else matrix(F, d, 0))
        C = block_matrix(F, [blocks], subdivide=False) if Z.ncols() else matrix(F, d, 0)
        kap = d - C.rank()
        return W1, kap

    def type_Fp(self, e):
        """CM type of A_p assuming f O_p = r^e with e <= K."""
        assert e <= self.K
        F = self.F
        Be = self.image_mod(e)
        psi = self.psi
        P = matrix(F, [self.vec_k(self.poly(r) % psi, 1) for r in Be.rows()])
        ker = P.left_kernel().basis_matrix()
        mm = ker * Be
        n = e * self.d
        Z = Be.right_kernel().basis_matrix().transpose()
        if Z.ncols() == 0:
            return 0
        if mm.nrows() == 0:
            return n - Be.nrows()
        blocks = [self.mult_matrix(self.poly(b), e) * Z for b in mm.rows()]
        C = block_matrix(F, [blocks], subdivide=False)
        Tdim = n - C.rank()
        return Tdim - Be.nrows()
