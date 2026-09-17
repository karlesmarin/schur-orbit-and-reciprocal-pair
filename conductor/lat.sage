# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# lat.sage -- integer-lattice computation of A_m(q) in O = Z[alpha]: index, conductor exponent e at p (r^i in A_p test
#             by lattice sum), layers W_i by intersections, Cohen-Macaulay type by colon lattice; kappa and F_p type by capas.
# What: independent of gor.sage (Dickson recursion in Z[X]/Psi_q, HNF closure); cross-checks e, W, type.
# Run (in this directory): sage lat.sage m q p [colon]
import sys, time
load("capas.sage")

def hnf_rows(rows):
    M = matrix(ZZ, rows).echelon_form(include_zero_rows=False)
    return M

def build(m, q):
    R = PolynomialRing(ZZ, 'X'); X = R.gen()
    P = Psi_ZZ(q); D = P.degree()
    eta_prev, eta = R(2), X
    coef = [R(1)]
    for j in range(1, m + 1):
        new = [R(0)] * (len(coef) + 1)
        for i, c in enumerate(coef):
            new[i + 1] += c
            new[i] = (new[i] - eta * c) % P
        coef = new
        eta_prev, eta = eta, (X * eta - eta_prev) % P
    def vec(f):
        l = (f % P).list(); return l + [0] * (D - len(l))
    G = hnf_rows([vec(R(1))] + [vec(c) for c in coef[:-1]])
    gb = [R(list(r)) for r in G.rows()]
    def mulmat(g):
        rows = []; h = g % P
        for i in range(D):
            rows.append(vec(h)); h = (X * h) % P
        return matrix(ZZ, rows)
    mats = [mulmat(g) for g in gb]
    L = hnf_rows([vec(R(1))])
    it = 0
    while True:
        it += 1
        rows = [L]
        S = L.stack(G)
        for Mg in mats:
            S = S.stack(L * Mg)
        L2 = S.echelon_form(include_zero_rows=False)
        if L2 == L:
            break
        L = L2
    assert L.nrows() == D
    return R, X, P, D, L, vec, mulmat, it

def analyse(m, q, p, colon=False):
    res = {}
    t0 = time.time()
    R, X, P, D, A, vec, mulmat, it = build(m, q)
    v, n = split_q(q, p); E = euler_phi(p**v); d = euler_phi(n) // 2
    idx = abs(A.det())
    print("A_%d(%d) D=%d closure_iters=%d index=%s  [%.1fs]" % (m, q, D, it, factor(idx), time.time() - t0)); sys.stdout.flush()
    psi = R(Psi_ZZ(n))
    vA = val_p(idx, p)
    def rad_rows(i):
        rows = []
        for s in range(i + 1):
            base = (p**s) * psi**(i - s) % P
            h = base
            for j in range(D):
                rows.append(vec(h)); h = (X * h) % P
        return rows
    e = None
    for i in range(0, 4 * E * d + 5):
        S = hnf_rows(list(A.rows()) + rad_rows(i))
        if val_p(S.det(), p) == vA:
            e = i; break
    print("  p=%d E=%d d=%d v_p[O:A]=%d  e(lattice sum test)=%s  [%.1fs]" % (p, E, d, vA, e, time.time() - t0)); sys.stdout.flush()
    # sanity: r^(e-1) not in A_p (by construction) ; layers by intersection
    ZD = ZZ**D
    LA = ZD.span(A.rows())
    Ai = []
    for i in range(0, e + 2):
        Ri = ZD.span(hnf_rows(rad_rows(i)).rows()) if i > 0 else ZD
        Ai.append(LA.intersection(Ri))
    W = [val_p(Ai[i + 1].basis_matrix().det(), p) - val_p(Ai[i].basis_matrix().det(), p) for i in range(e + 1)]
    lOA = vA; lAf = e * d - vA
    print("  W(lattice)=%s  l(O/A)=%d  l(A/f)=%d  sum(d-W_i,i<e)=%d sum(W_i,i<e)=%d  [%.1fs]" % (W, lOA, lAf, sum(d - w for w in W[:e]), sum(W[:e]), time.time() - t0)); sys.stdout.flush()
    res.update(dict(D=D, E=E, d=d, e=e, W=W, lOA=lOA, lAf=lAf))
    if e <= E and E >= 2:
        C = Capas(q, p, m, max(2, e))
        Wc = C.W()
        W1, kap = C.kappa()
        tF = C.type_Fp(e)
        print("  capas: W=%s W1=%d kappa=%d type_Fp=%d  pred(d-W1+kappa-1)=%d  [%.1fs]" % (Wc, W1, kap, tF, d - W1 + kap - 1, time.time() - t0)); sys.stdout.flush()
        res.update(dict(W1=W1, kappa=kap, typeF=tF))
    if colon:
        mx = Ai[1]
        Hinv = A.change_ring(QQ).inverse()
        T = ZD
        for b in mx.basis():
            bb = R(list(b))
            Mb = mulmat(bb).change_ring(QQ)
            N = Mb * Hinv
            den = lcm([x.denominator() for x in N.list()])
            # x in T iff x*N integral
            Lb = (QQ**D).span((N.inverse()).rows(), ZZ) if N.det() != 0 else None
            T = T.intersection(Lb)
        typ = val_p(idx, p) - val_p(T.basis_matrix().det(), p)
        print("  type(colon lattice)=%d  HS12.2.3: l(O/A)-l(A/f)=%d >= type-1=%d : %s  [%.1fs]" % (typ, lOA - lAf, typ - 1, lOA - lAf >= typ - 1, time.time() - t0))
        res["typeC"] = typ
    return res

if sys.argv and sys.argv[0].endswith("lat.sage.py"):
    a_ = sys.argv[1:]
    analyse(int(a_[0]), int(a_[1]), int(a_[2]), len(a_) > 3 and a_[3] == "colon")
