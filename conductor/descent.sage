# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# descent.sage -- item 1: ladder below P=p^a and descent W_P(A_m(q)) = W_1(A_{m0}(q0)) = rank_Fp[S_0(u^-1 c)].
# What: 2m+1 = P*M*n, q = p^v n, q0 = q/P, m0 = (Mn-1)/2; layers by inst.sage with K=P+1; ring check dB[P]==1+a+W_1(A0), m^2=0.
# Run: sh sg.sh descent.sage SHARD NSHARDS COSTCAP     (cost proxy m*D; cases above cap are listed as SKIP)
load("inst.sage")
shard, nsh, cap = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
ncomp = [9, 15, 21, 25, 27, 33, 35, 39, 45]
nprim = [5, 7, 11, 13, 17, 19, 23, 29, 31]
grid = []
for p in [3, 5, 7, 11]:
    for a in ([1, 2, 3] if p == 3 else [1, 2]):
        P = p**a
        for n in ncomp + nprim:
            if n % p == 0: continue
            for M in [1, 3, 5, 7, 11, 13]:
                if M % p == 0: continue
                v0 = a + 1
                while p**(v0 - a) <= M:
                    v0 += 1
                for v in [v0, v0 + 1]:
                    m = (P * M * n - 1) // 2
                    D = (P + 1) * euler_phi(n)
                    cost = m * D
                    if v == v0 + 1 and cost > cap / 20:
                        continue
                    grid.append((cost, p, a, n, M, v))
grid.sort()
mine = [g for i, g in enumerate(grid) if i % nsh == shard]
tot = fl = flP = flR = flQ = sk = 0
for (cost, p, a, n, M, v) in mine:
    P = p**a; m = (P * M * n - 1) // 2; q = p**v * n; d = euler_phi(n) // 2
    if cost > cap:
        sk += 1
        print("SKIP p=%d a=%d n=%d M=%d v=%d m=%d cost=%d" % (p, a, n, M, v, m, cost)); sys.stdout.flush()
        continue
    t0 = time.time()
    assert 2 * m < q and euler_phi(p**v) >= P + 1
    I = Img(q, p, m, P + 1)
    dB = I.dims(); W = [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, len(dB))]
    pred = [1] + [1 if any(i == P - p**j for j in range(a)) else 0 for i in range(1, P)]
    ok = (W[:P] == pred)
    m0 = (M * n - 1) // 2; q0 = p**(v - a) * n
    I0 = Img(q0, p, m0, 2); W0 = I0.W()
    rk, nu = moment_rank(n, m0, p)
    okP = (W[P] == W0[1] == rk)
    okR = (dB[P] == 1 + a + W0[1])
    sq, dm = I.max_ideal_square_zero()
    tot += 1; fl += (not ok); flP += (not okP); flR += (not okR); flQ += (not sq)
    print("p=%d a=%d n=%2d M=%2d v=%d q=%d m=%d d=%d | W=%s | ladder %s | W_P=%d W1(A_%d(%d))=%d rank=%d %s | dB[P]=%d 1+a+W1=%d %s | m^2=0 %s | %.1fs" % (
        p, a, n, M, v, q, m, d, W, "OK" if ok else "FAIL", W[P], m0, q0, W0[1], rk, "OK" if okP else "FAIL",
        dB[P], 1 + a + W0[1], "OK" if okR else "FAIL", sq, time.time() - t0))
    sys.stdout.flush()
print("SHARD %d/%d TOTAL run=%d skipped=%d ladder_fail=%d descent_fail=%d ring_dim_fail=%d msq_fail=%d" % (shard, nsh, tot, sk, fl, flP, flR, flQ))
