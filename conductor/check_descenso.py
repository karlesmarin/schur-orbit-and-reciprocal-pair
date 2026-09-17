# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Control del enunciado de descenso: W_P(A_m(q)) = W_1(A_{m0}(q0)), pi^P notin B, B subset R_rho + pi^P V0.
import sys
import numpy as np

def polmul(a, b, p):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] = (r[i+j] + x*y) % p
    return r

def pdiv(num, den):
    num = num[:]; out = [0]*(len(num)-len(den)+1)
    for i in range(len(num)-len(den), -1, -1):
        c = num[i+len(den)-1] // den[-1]
        out[i] = c
        for j, d in enumerate(den):
            num[i+j] -= c*d
    assert all(x == 0 for x in num[:len(den)-1])
    return out

def cyclo(n):
    phis = {}
    for d in range(1, n+1):
        if n % d == 0:
            f = [-1] + [0]*(d-1) + [1]
            for e in range(1, d):
                if d % e == 0:
                    f = pdiv(f, phis[e])
            phis[d] = f
    return phis[n]

class Ring:
    def __init__(self, p, n, K):
        self.p = p
        mod = [1]
        ph = [c % p for c in cyclo(n)]
        for _ in range(K):
            mod = polmul(mod, ph, p)
        self.D = D = len(mod)-1
        Y = np.zeros((D, D), dtype=np.int64)
        for k in range(D):
            if k+1 < D:
                Y[k+1, k] = 1
            else:
                for t in range(D):
                    Y[t, k] = (-mod[t]) % p
        self.Y = Y
    def one(self):
        v = np.zeros(self.D, dtype=np.int64); v[0] = 1
        return v
    def ypow(self, e, q):
        e %= q
        r = np.eye(self.D, dtype=np.int64); M = self.Y.copy()
        while e:
            if e & 1:
                r = (M @ r) % self.p
            M = (M @ M) % self.p; e >>= 1
        return (r @ self.one()) % self.p
    def mat(self, v):
        cols = []; w = v.copy()
        for k in range(self.D):
            cols.append(w); w = (self.Y @ w) % self.p
        return np.array(cols).T % self.p

def rank(rows, p):
    if len(rows) == 0:
        return 0
    M = np.array(rows, dtype=np.int64) % p
    r = 0; R, C = M.shape
    for c in range(C):
        piv = next((i for i in range(r, R) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        M[r] = (M[r]*pow(int(M[r, c]), -1, p)) % p
        for i in range(R):
            if i != r and M[i, c]:
                M[i] = (M[i] - M[i, c]*M[r]) % p
        r += 1
        if r == R:
            break
    return r

def algebra(R, gens):
    p = R.p
    basis = []
    for v in [R.one()] + list(gens):
        if rank(basis + [v % p], p) > len(basis):
            basis.append(v % p)
    mats = [R.mat(g) for g in gens]
    queue = list(basis)
    while queue:
        v = queue.pop()
        for Mt in mats:
            w = (Mt @ v) % p
            if rank(basis + [w], p) > len(basis):
                basis.append(w); queue.append(w)
    return basis

def gens_A(R, m, q):
    p = R.p
    coef = [R.one()]
    for j in range(1, m+1):
        Me = R.mat((R.ypow(j, q) + R.ypow(-j, q)) % p)
        new = [np.zeros(R.D, dtype=np.int64) for _ in range(len(coef)+1)]
        for i, c in enumerate(coef):
            new[i+1] = (new[i+1] + c) % p
            new[i] = (new[i] - Me @ c) % p
        coef = new
    return coef[:-1]

def reduce_to(vecs, p, n, k):
    mod = [1]; ph = [c % p for c in cyclo(n)]
    for _ in range(k):
        mod = polmul(mod, ph, p)
    Dk = len(mod)-1
    out = []
    for v in vecs:
        poly = [int(x) for x in v]
        for deg in range(len(poly)-1, Dk-1, -1):
            c = poly[deg]
            if c:
                for t in range(Dk+1):
                    poly[deg-Dk+t] = (poly[deg-Dk+t] - c*mod[t]) % p
        out.append(poly[:Dk])
    return out

def layers(B, p, n, K):
    dims = [rank(reduce_to(B, p, n, k), p) for k in range(1, K+1)]
    return [dims[0]] + [dims[i]-dims[i-1] for i in range(1, K)]

def crt(a1, m1, a2, m2):
    for x in range(m1*m2):
        if x % m1 == a1 and x % m2 == a2:
            return x

def signed(r, n):
    r %= n
    return r if r <= n//2 else r-n

def caso(p, a, n, M, extra_v=0):
    P = p**a
    m = (P*M*n-1)//2
    v = a+1
    while p**v*n <= 2*m:
        v += 1
    v += extra_v
    q = p**v*n
    assert 2*m < q and v > a and M % p and n % p and (M*n) % 2
    K = P+1
    R = Ring(p, n, K)
    B = algebra(R, gens_A(R, m, q))
    W = layers(B, p, n, K)
    erho = crt(0, n, 1, p**v); eom = crt(1, n, 0, p**v)
    rho = R.ypow(erho, q); om = R.ypow(eom, q)
    one = R.one()
    Mpi = R.mat((rho - one) % p)
    powpi = [one]
    for _ in range(P):
        powpi.append((Mpi @ powpi[-1]) % p)
    piP = powpi[P]
    Mom = R.mat(om)
    ompow = [one]
    for _ in range(n-1):
        ompow.append((Mom @ ompow[-1]) % p)
    G = []
    for j in range(n):
        g = np.zeros(R.D, dtype=np.int64)
        for r in range(n):
            g = (g + (M*signed(r, n)) * ompow[(P*r*j) % n]) % p
        G.append(g)
    MpiP = R.mat(piP)
    piPG = [(MpiP @ g) % p for g in G]
    w = rank(piPG, p)
    dimB = len(B)
    piP_in_B = rank(B + [piP], p) == dimB
    V0_in_B = rank(B + piPG, p) == dimB
    big = rank(B + powpi + piPG, p)
    ctrl = rank(powpi + piPG, p)
    m0 = (M*n-1)//2; q0 = q//P
    R0 = Ring(p, n, 2)
    W0 = layers(algebra(R0, gens_A(R0, m0, q0)), p, n, 2)
    reps = []
    for u in range(1, n):
        if np.gcd(u, n) == 1 and (n-u) not in reps:
            reps.append(u)
    Smat = [[(M*signed(pow(u, -1, n)*c, n)) % p for c in range(1, n)] for u in reps]
    rS = rank(Smat, p)
    lad = [i for i in range(1, P) if W[i]]
    ladder_ok = lad == sorted(P-p**j for j in range(a)) and all(W[i] == 1 for i in lad)
    ok = (W[P] == W0[1] == w == rS) and (not piP_in_B) and V0_in_B and big == ctrl == P+1+w and ladder_ok and dimB == 1 + a + w
    print("p=%d a=%d n=%d M=%d m=%d q=%d v=%d D=%d | W=%s | W1(A0)=%d rankS=%d dimV0=%d | piP in B: %s | piP.V0 in B: %s | dim(B+Rrho+piPV0)=%d ctrl=%d | dimB=%d 1+a+w=%d | ladder %s | %s" % (p, a, n, M, m, q, v, R.D, W, W0[1], rS, w, piP_in_B, V0_in_B, big, ctrl, dimB, 1+a+w, ladder_ok, "OK" if ok else "FALLA"), flush=True)
    return ok

if __name__ == "__main__":
    casos = [tuple(int(x) for x in s.split(",")) for s in sys.argv[1:]]
    res = []
    for c in casos:
        try:
            res.append((c, caso(*c)))
        except Exception as e:
            print(c, "ERROR", repr(e), flush=True); res.append((c, None))
    print("RESUMEN", sum(1 for _, r in res if r), "/", len(res), "OK; no OK:", [c for c, r in res if not r])
