# ONE PASS over every load-bearing formula of the session, re-derived from scratch.
#
# Nothing here reuses an earlier script's output. Each claim is rebuilt from its own definition and
# checked against an independently computed object. A claim that is not checked here is not claimed.
#
# Carles Marin + Claude (AI assistant).  2026-07-30.

FAIL = {}
def rec(name, bad, n):
    FAIL[name] = (bad, n)
    print("  %-58s %6d checks  %s" % (name, n, "OK" if bad == 0 else "FAIL x%d" % bad))

# ---------- the object, defined once, from the bialternant ----------
def obj(lam, r, ring=None, alpha=None):
    N = 2 * r + 2
    if ring is None:
        ring = LaurentPolynomialRing(QQ, ['z%d' % i for i in range(1, r + 1)])
        alpha = [ring(1), ring(-1)]
        for z in ring.gens():
            alpha += [z, z**-1]
    L = list(lam) + [0] * max(0, N - len(lam))
    if len([x for x in lam if x > 0]) > N:
        return ring(0)
    L = L[:N]
    D = max(L[0], 1) + N + 1
    PK = PowerSeriesRing(ring, 'u', default_prec=D + 2); u = PK.gen()
    G = prod(1 / (1 - a * u) for a in alpha)
    h = [G[k] for k in range(D + 1)]
    H = lambda k: ring(0) if k < 0 else h[k]
    return matrix(ring, N, N, lambda i, j: H(L[i] - i + j)).det()

print("FORMULA AUDIT -- every load-bearing claim, rebuilt and rechecked\n")

# ---------- 1. the PUBLISHED r=1 closed form: eps * chi_a chi_b chi_c ----------
# beta = lam + rho, E = even part, O = odd part.  e=o=2: E={p>q}, O={s>t}
#   spins a=(p-q)/2-1, b=(s-t)/2-1, c=|p+q-s-t|/2-1 ; chi_k(z) = (z^{k+1}-z^{-k-1})/(z-1/z)
R1 = LaurentPolynomialRing(QQ, 'z'); z1 = R1.gen()
A1 = [R1(1), R1(-1), z1, z1**-1]
def chi(k):
    if k < 0:
        return R1(0)
    return sum(z1**(k - 2 * j) for j in range(k + 1))
bad = n = 0
for S in range(0, 22):
    for lam in Partitions(S, max_length=4):
        L = list(lam) + [0] * (4 - len(lam))
        beta = [L[i] + (3 - i) for i in range(4)]
        E = sorted([b for b in beta if b % 2 == 0], reverse=True)
        O = sorted([b for b in beta if b % 2 == 1], reverse=True)
        if len(E) != 2 or len(O) != 2:
            continue
        p, q = E; ss, t = O
        a = (p - q) // 2 - 1; b = (ss - t) // 2 - 1; c = abs(p + q - ss - t) // 2 - 1
        pred = chi(a) * chi(b) * chi(c)
        got = obj(lam, 1, R1, A1)
        n += 1
        if got != pred and got != -pred:
            bad += 1
rec("1. published r=1 closed form  eps*chi_a chi_b chi_c (2-2 split)", bad, n)

# ---------- 2. branch (a): all beta one parity => 0, via the squared alphabet ----------
bad = n = 0
for r in [1, 2, 3]:
    N = 2 * r + 2
    for S in range(0, 15):
        for lam in Partitions(S, max_length=N):
            L = list(lam) + [0] * (N - len(lam))
            beta = [L[i] + (N - 1 - i) for i in range(N)]
            if len(set(x % 2 for x in beta)) != 1:
                continue
            n += 1
            if obj(lam, r) != 0:
                bad += 1
rec("2. branch (a): constant-parity beta => object = 0", bad, n)

# ---------- 3. branch (b): lam self-complementary of ODD width => 0 ; EVEN width => != 0 ----------
bad = n = badc = nc = 0
for r in [1, 2, 3]:
    N = 2 * r + 2
    for S in range(0, 19 - 2 * r):
        for lam in Partitions(S, max_length=N):
            L = list(lam) + [0] * (N - len(lam))
            if len(L) != N:
                continue
            w = L[0] + L[N - 1]
            if not all(L[i] + L[N - 1 - i] == w for i in range(N)):
                continue
            beta = [L[i] + (N - 1 - i) for i in range(N)]
            if len(set(x % 2 for x in beta)) == 1:
                continue                      # branch (a), covered above
            if w % 2 == 1:
                n += 1
                if obj(lam, r) != 0:
                    bad += 1
            else:
                nc += 1
                if obj(lam, r) == 0:
                    badc += 1
rec("3. branch (b): self-complementary, ODD width  => object = 0", bad, n)
rec("3b. CONTROL: self-complementary, EVEN width  => object != 0", badc, nc)

# ---------- 4. the criterion, BOTH directions ----------
import random
random.seed(int(7))
bad = n = 0
for r, NMAX in [(1, 20), (2, 16), (3, 13)]:
    N = 2 * r + 2
    for S in range(0, NMAX + 1):
        for lam in Partitions(S, max_length=N):
            L = list(lam) + [0] * (N - len(lam))
            beta = [L[i] + (N - 1 - i) for i in range(N)]
            ba = len(set(x % 2 for x in beta)) == 1
            w = L[0] + L[N - 1] if len(L) == N else None
            sc = (len(L) == N and w % 2 == 1 and
                  all(L[i] + L[N - 1 - i] == w for i in range(N)))
            n += 1
            if (ba or sc) != (obj(lam, r) == 0):
                bad += 1
rec("4. criterion  (a) or (b)  <=>  object = 0,  r=1,2,3", bad, n)

# ---------- 5. the complementation identity  s_{lam-hat}(A) = (-1)^c * (-1) * s_lam(A) ----------
bad = n = 0
for r in [1, 2]:
    N = 2 * r + 2
    for S in range(0, 13 - 2 * r):
        for lam in Partitions(S, max_length=N):
            L = list(lam) + [0] * (N - len(lam))
            if len(L) != N:
                continue
            v = obj(lam, r)
            for w in range(L[0], L[0] + 3):
                if any(x > w for x in L):
                    continue
                hh = [w - L[N - 1 - i] for i in range(N)]
                c = w + N - 1
                n += 1
                if obj(hh, r) != (-1)**c * (-1) * v:
                    bad += 1
rec("5. complementation  s_{lam^}(A) = (-1)^c (-1) s_lam(A)", bad, n)

# ---------- 6. Phi_A anti-self-reciprocal, middle coefficient 0 ----------
bad = n = 0
for r in [1, 2, 3]:
    N = 2 * r + 2
    F = FractionField(PolynomialRing(QQ, ['z%d' % i for i in range(1, r + 1)]))
    Px = PolynomialRing(F, 'x'); x = Px.gen()
    A = [F(1), F(-1)]
    for zz in F.gens():
        A += [zz, 1 / zz]
    Phi = prod(x - a for a in A)
    co = [Phi[k] for k in range(N + 1)]
    n += 1
    if not (all(co[k] == -co[N - k] for k in range(N + 1)) and co[N // 2] == 0):
        bad += 1
rec("6. Phi_A anti-self-reciprocal, middle coefficient zero", bad, n)

# ---------- 7. (SP): object = 0  <=>  a sparse multiple of Phi_A supported in beta exists -------
bad = n = 0
for r in [1, 2]:
    N = 2 * r + 2
    F = FractionField(PolynomialRing(QQ, ['z%d' % i for i in range(1, r + 1)]))
    Px = PolynomialRing(F, 'x'); x = Px.gen()
    A = [F(1), F(-1)]
    for zz in F.gens():
        A += [zz, 1 / zz]
    Phi = prod(x - a for a in A)
    co = [Phi[k] for k in range(N + 1)]
    for S in range(0, 13 - 2 * r):
        for lam in Partitions(S, max_length=N):
            L = list(lam) + [0] * (N - len(lam))
            beta = [L[i] + (N - 1 - i) for i in range(N)]
            b1 = beta[0]; D = b1 - N
            if D < 0:
                continue
            offs = [e for e in range(0, b1 + 1) if e not in beta]
            M = matrix(F, len(offs), D + 1)
            for i, e in enumerate(offs):
                for t in range(D + 1):
                    k = e - t
                    M[i, t] = co[k] if 0 <= k <= N else F(0)
            n += 1
            if (M.rank() < D + 1) != (obj(lam, r) == 0):
                bad += 1
rec("7. (SP) sparse-multiple form <=> vanishing", bad, n)

# ---------- 8. d'Ocagne for the Chebyshev-like S_n ----------
Rc = PolynomialRing(QQ, 'c'); cc = Rc.gen()
_S = [Rc(0), Rc(1)]
for k in range(2, 70):
    _S.append(cc * _S[k - 1] - _S[k - 2])
Sv = lambda k: _S[k] if k >= 0 else -_S[-k]
bad = n = 0
for b in range(1, 32):
    for bp in range(0, b):
        n += 1
        if Sv(b) * Sv(bp - 1) - Sv(bp) * Sv(b - 1) != -Sv(b - bp):
            bad += 1
rec("8. d'Ocagne  S_b S_{b'-1} - S_{b'} S_{b-1} = -S_{b-b'}", bad, n)
n2 = sum(1 for k in range(1, 40))
rec("8b. deg_c S_n = n-1", sum(1 for k in range(1, 40) if Sv(k).degree() != k - 1), n2)

# ---------- 9. reflection stability of the sparse locus ----------
bad = n = 0
for r in [1, 2]:
    N = 2 * r + 2
    F = FractionField(PolynomialRing(QQ, ['z%d' % i for i in range(1, r + 1)]))
    Px = PolynomialRing(F, 'x'); x = Px.gen()
    A = [F(1), F(-1)]
    for zz in F.gens():
        A += [zz, 1 / zz]
    Phi = prod(x - a for a in A)
    for D in range(0, 5):
        M = D + N
        for _ in range(12):
            coeffs = [F(int(random.randint(-3, 3))) for _ in range(D + 1)]
            if all(c == 0 for c in coeffs):
                continue
            Rp = Px(coeffs)
            P = Rp * Phi
            Rstar = Px([Rp[D - t] for t in range(D + 1)])
            Pstar = Px([P[P.degree() - t] for t in range(P.degree() + 1)]) if P != 0 else Px(0)
            n += 1
            # (R Phi)^{*M} = R^{*D} * (-Phi)  as polynomials, when deg(R) = D
            if Rp[D] != 0 and Pstar != -(Rstar * Phi):
                bad += 1
rec("9. reflection: (R Phi)^{*M} = R^{*D} (-Phi)", bad, n)

print("\n" + "=" * 74)
tot_bad = sum(v[0] for v in FAIL.values()); tot_n = sum(v[1] for v in FAIL.values())
print("TOTAL: %d checks, %d failures" % (tot_n, tot_bad))
for k, v in FAIL.items():
    if v[0]:
        print("   FAILED:", k)
