# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# item2.sage -- staircase before p^a: W_0=1, W_i = [i = P - p^j, 0<=j<a] for 1<=i<P; and W_P vs rank of moment matrix of m'=(Mn-1)/2.
# What: layers by capas.sage (F_p[X]/Psi_n^K), K=P+1<=E, over a grid of (p,a,n,M,v).
# Run (in this directory): sage item2.sage
import sys, time
load("capas.sage")

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

grid = []
for (p, a, ns, Ms) in [(3, 1, [5, 7, 11, 13, 17], [1, 5, 7]),
                       (3, 2, [5, 7, 11, 13], [1, 5, 7]),
                       (3, 3, [5, 7, 11, 13], [1, 5]),
                       (5, 1, [7, 9, 11, 13], [1, 3, 7]),
                       (5, 2, [7, 9, 11, 13], [1, 3]),
                       (7, 1, [5, 9, 11, 13, 15], [1, 3, 5])]:
    for n in ns:
        for M in Ms:
            N = p**a * M
            v0 = a + 1
            while p**v0 <= N:
                v0 += 1
            vs = [v0]
            if p**a <= 9 and n <= 13:
                vs.append(v0 + 1)
            for v in vs:
                grid.append((p, a, n, M, v))

fails = 0; failsP = 0; tot = 0
for (p, a, n, M, v) in grid:
    t0 = time.time()
    P = p**a; N = P * M; m = (N * n - 1) // 2; q = p**v * n
    E = euler_phi(p**v); d = euler_phi(n) // 2
    assert 2 * m < q and E > P
    K = P + 1
    C = Capas(q, p, m, K)
    W = C.W()
    pred = [1] + [1 if any(i == P - p**j for j in range(a)) else 0 for i in range(1, P)]
    ok = (W[:P] == pred)
    mp = (M * n - 1) // 2
    rk, nu = moment_rank(n, mp, p)
    okP = (W[P] == rk)
    tot += 1; fails += (not ok); failsP += (not okP)
    print("p=%d a=%d n=%2d M=%d v=%d q=%d m=%d d=%d | W[0..P]=%s | stair %s | W_P=%d rank(m'=%d)=%d %s [%.1fs]" % (
        p, a, n, M, v, q, m, d, W, "OK" if ok else "FAIL pred=%s" % pred, W[P], mp, rk, "OK" if okP else "FAIL", time.time() - t0))
    sys.stdout.flush()
print("TOTAL cases=%d staircase_fail=%d W_P_rank_fail=%d" % (tot, fails, failsP))
