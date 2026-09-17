# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# dos.sage -- item 2: factor at 2 for odd n (composite and prime): 2m+1=Nn, p odd, p !| n*N*h^-(Q(zeta_n)).
# What: layers W_0..W_2 (inst.sage, K=3) vs prediction from H=<2> in (Z/n)^x: -1 in H -> W1=d; else delta=d/t, W1=d-delta, W2=d;
#       rank of [S(u^-1 c)] mod p; F_p type (capas.sage) vs 2*delta-1 on a sample.   Run: sh sg.sh dos.sage
load("capas.sage")
load("inst.sage")
hm = {}
for ln in open("stickelberger_smith_OUT.txt"):
    t = ln.split()
    if t and t[0].isdigit():
        hm[int(t[0])] = Integer(t[4])
ncomp = [9, 15, 21, 25, 27, 33, 35, 39, 45, 49, 51, 55, 57, 63, 65, 69, 75, 77, 85, 87]
nprim = [n for n in prime_range(5, 90)]
tot = bad = 0; badrk = 0
typ_done = {}; typ_tot = typ_bad = 0
struct = {}
for n in ncomp + nprim:
    G = [u for u in range(1, n) if gcd(u, n) == 1]
    t = Mod(2, n).multiplicative_order()
    H = set(power_mod(2, i, n) for i in range(t))
    minus = (n - 1) in H
    d = euler_phi(n) // 2
    delta = None if minus else d // t
    assert minus or d % t == 0
    struct[n] = (t, minus, d, delta)
    print("n=%d %s t=%d -1inH=%s d=%d delta=%s h^-=%s" % (n, "comp" if n in ncomp else "prime", t, minus, d, delta, hm[n])); sys.stdout.flush()
    for p in prime_range(3, 32):
        if n % p == 0 or hm[n] % p == 0: continue
        for N in [1, 3, 5]:
            if N % p == 0: continue
            v = 1
            while euler_phi(p**v) < 3 or p**v <= N:
                v += 1
            q = p**v * n; m = (N * n - 1) // 2
            assert 2 * m < q
            I = Img(q, p, m, 3); W = I.W()
            pred = [1, d, d] if minus else [1, d - delta, d]
            rk, _ = moment_rank(n, m, p)
            ok = (W == pred); okr = (rk == W[1])
            tot += 1; bad += (not ok); badrk += (not okr)
            extra = ""
            if not minus:
                extra = " v_p[O:A]=%d v_p[A:f]=%d" % (d + delta - 1, d - delta + 1)
            line = "  p=%2d N=%d v=%d q=%d m=%d W=%s pred=%s %s rank[S]=%d %s%s" % (p, N, v, q, m, W, pred, "OK" if ok else "FAIL", rk, "OK" if okr else "FAIL", extra)
            # type sample: first N=1 case per (n) with -1 not in H, up to 2 primes, only if cheap
            key = n
            if (not minus) and N == 1 and typ_done.get(key, 0) < 2 and euler_phi(q) // 2 * 2 <= 200 and euler_phi(p**v) >= 2:
                C = Capas(q, p, m, 2)
                W1c, kap = C.kappa(); tF = C.type_Fp(2)
                okt = (tF == 2 * delta - 1 and kap == delta and W1c == d - delta)
                typ_tot += 1; typ_bad += (not okt); typ_done[key] = typ_done.get(key, 0) + 1
                line += " | capas W1=%d kappa=%d type_Fp=%d pred=%d %s" % (W1c, kap, tF, 2 * delta - 1, "OK" if okt else "FAIL")
            print(line); sys.stdout.flush()
print("TOTAL layer_cases=%d layer_fail=%d rank_fail=%d | type_cases=%d type_fail=%d" % (tot, bad, badrk, typ_tot, typ_bad))
print("t even and -1 not in H:", [(n, struct[n][0], struct[n][2], struct[n][3]) for n in ncomp + nprim if struct[n][0] % 2 == 0 and not struct[n][1]])
print("t odd:", [(n, struct[n][0], struct[n][2], struct[n][3]) for n in ncomp + nprim if struct[n][0] % 2 == 1])
print("delta>=2:", [(n, struct[n][0], struct[n][3]) for n in ncomp + nprim if struct[n][3] and struct[n][3] >= 2])
