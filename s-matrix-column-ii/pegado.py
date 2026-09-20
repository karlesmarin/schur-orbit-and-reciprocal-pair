# P3_pegado.py -- author: Carles Marin (Claude as AI assistant)
# What: the p-adic anatomy of the fusion ring R = Fus(SU(n)_k), q = n+k = p prime.
#   R  ⊆  prod_j A_j  ⊆  prod_j O_j   (j = Galois orbits of columns, A_j = ring of one column)
#   v_p disc R = sum_j v_p disc K_j + 2 sum_j v_p[O_j:A_j] + 2 v_p[prod A_j : R]      (standard)
#   Prediction (Gepner hep-th/0608140 eq.(13) D = prod S_{0mu}^{-2}, Vandermonde v_p = (n-1)(k-1)/(p-1)
#   per column; ramification (e-1)/e and defect delta/e per column, e = (p-1)/h):
#     2 v_p[prod A_j : R] = sum_mu [ (n-1)(k-1)/(p-1) - (e_mu - 1 + 2 delta_mu)/e_mu ].
#   Test: compute [prod A_j : R] DIRECTLY as a lattice index in L^J (L = Q(zeta_{nq})), with the
#   Schur values s_lambda(g_mu) by the bialternant, and compare its p-part with the prediction.
#   Control: also v_p of the full disc(R) from the Gram determinant of the trace form.
# Run: MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/w" -w /w sage-normaliz:local sage -python -u /w/P3_pegado.py
import sys, itertools
from sage.all import CyclotomicField, matrix, ZZ, QQ, gcd, euler_phi, Integer


def weights(n, k):
    for mu in itertools.product(range(k + 1), repeat=n - 1):
        if sum(mu) <= k:
            yield mu


def particion(mu):
    # Dynkin labels -> partition with n-1 rows
    n1 = len(mu)
    return [sum(mu[i:]) for i in range(n1)]


def exps(n, k, mu):
    q = n + k; m = n * q
    b = [0] * n
    for i in range(n - 2, -1, -1):
        b[i] = b[i + 1] + mu[i] + 1
    S = sum(b)
    return [(n * bj - S) % m for bj in b]


def semigrupo(T, p):
    H = [a for a in range(1, p) if sorted(a * t % p for t in T) == sorted(T)]
    h = len(H); e = (p - 1) // h
    G = sorted(set([e] + [j for j in range(1, e) if sum(pow(t, h * j, p) for t in T) % p]))
    lim = G[0] * G[-1] + 2
    alc = [False] * lim; alc[0] = True
    for x in range(1, lim):
        alc[x] = any(x >= g and alc[x - g] for g in G)
    return e, sum(1 for x in range(1, lim) if not alc[x])


def espectro(n, k, mu):
    p = n + k
    b = [0] * n
    for i in range(n - 2, -1, -1):
        b[i] = b[i + 1] + mu[i] + 1
    s = sum(b) * pow(n, -1, p) % p
    return sorted((x - s) % p for x in b)


def vp(x, p):
    x = Integer(x); v = 0
    while x % p == 0 and x != 0:
        x //= p; v += 1
    return v


def run(n, k):
    p = n + k; N = n * p
    L = CyclotomicField(N); z = L.gen(); d = L.degree()
    cols = list(weights(n, k))
    lams = cols  # integrable weights = basis of R
    parts = [particion(l) + [0] for l in lams]
    # eigenvalues of g_mu
    X = {mu: [z ** c for c in exps(n, k, mu)] for mu in cols}
    def schur(part, xs):
        num = matrix(L, n, n, lambda i, j: xs[i] ** (part[j] + n - 1 - j))
        den = matrix(L, n, n, lambda i, j: xs[i] ** (n - 1 - j))
        return num.det() / den.det()
    val = {mu: [schur(pt, X[mu]) for pt in parts] for mu in cols}
    # Galois orbits of columns: sigma_a acts on the value vector
    units = [a for a in range(1, N) if gcd(a, N) == 1]
    key = {mu: tuple(val[mu][i] for i in range(1, n)) for mu in cols}  # fundamental reps first? use all
    resto = set(cols); orbitas = []
    while resto:
        mu = min(resto)
        orb = set()
        for a in units:
            sa = L.hom([z ** a])
            v = tuple(sa(x) for x in val[mu])
            for nu in list(resto):
                if nu not in orb and tuple(val[nu]) == v:
                    orb.add(nu)
        orbitas.append((mu, orb)); resto -= orb
    reps = [mu for mu, _ in orbitas]
    def coords(x):
        return list(x.list()) + [0] * (d - len(x.list()))
    # lattice of R inside L^J, and of each A_j
    filasR = [sum((coords(val[mu][i]) for mu in reps), []) for i in range(len(lams))]
    def covol2(rows):
        M = matrix(QQ, rows)
        den = Integer(1)
        for r in rows:
            for c in r:
                den = den.lcm(QQ(c).denominator())
        Mz = (M * den).change_ring(ZZ)
        H = Mz.hermite_form(include_zero_rows=False)
        return (H * H.transpose()).det() / den ** (2 * H.nrows()), H.nrows()
    gR, rR = covol2(filasR)
    gA = QQ(1); rA = 0
    for mu in reps:
        g, r = covol2([coords(val[mu][i]) for i in range(len(lams))])
        gA *= g; rA += r
    idx2 = gR / gA
    idx = idx2.sqrt()
    # prediction
    pred2 = QQ(0)
    for mu in cols:
        e, delta = semigrupo(espectro(n, k, mu), p)
        pred2 += QQ((n - 1) * (k - 1)) / (p - 1) - QQ(e - 1 + 2 * delta) / e
    # disc(R) from trace form (control of Gepner)
    tr = [[sum(val[mu][i] * val[mu][j] for mu in cols) for j in range(len(lams))] for i in range(len(lams))]
    discR = matrix(QQ, [[QQ(x) for x in fila] for fila in tr]).det()
    print("SU(%d)_%d p=%d: columnas %d, orbitas %d, rango R %d / prod A %d" % (n, k, p, len(cols), len(reps), rR, rA))
    print("   [prod A_j : R] = %s = %s   v_p = %d   | prediccion 2*gl = %s -> gl = %s"
          % (idx, Integer(idx).factor() if idx in ZZ else idx, vp(idx, p) if idx in ZZ else -1, pred2, pred2 / 2))
    print("   v_p disc R (traza) = %d   | Gepner-Vandermonde: %s" % (vp(abs(discR.numerator()), p) - vp(discR.denominator(), p),
          QQ(len(cols) * (n - 1) * (k - 1)) / (p - 1)))
    sys.stdout.flush()


for (n, k) in [(2, 3), (3, 2), (3, 4), (4, 3), (2, 9), (3, 8)]:
    run(n, k)
