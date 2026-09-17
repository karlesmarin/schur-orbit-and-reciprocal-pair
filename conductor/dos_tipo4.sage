# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# dos_tipo4.sage -- item 2 addendum: F_p type (capas.sage) for the delta=4 cases n=85 (composite, t=8 even), 73, 89 (prime).
# What: W1, kappa, type_Fp at e=2 vs d-delta, delta, 2*delta-1; also layers W_0..W_2.   Run: sh sg.sh dos_tipo4.sage
import sys, time
load("capas.sage")
hm = {}
for ln in open("stickelberger_smith_OUT.txt"):
    t = ln.split()
    if t and t[0].isdigit():
        hm[int(t[0])] = Integer(t[4])
for n in [85, 73, 89, 51, 63]:
    d = euler_phi(n) // 2; t = Mod(2, n).multiplicative_order(); delta = d // t
    done = 0
    for p in prime_range(3, 32):
        if n % p == 0 or hm[n] % p == 0 or done >= 2: continue
        v = 1
        while euler_phi(p**v) < 3:
            v += 1
        q = p**v * n; m = (n - 1) // 2
        t0 = time.time()
        C = Capas(q, p, m, 3)
        W = C.W(); W1, kap = C.kappa(); tF = C.type_Fp(2)
        ok = (W == [1, d - delta, d] and kap == delta and tF == 2 * delta - 1)
        print("n=%d t=%d d=%d delta=%d p=%d q=%d m=%d | W=%s kappa=%d type_Fp=%d pred type=%d %s [%.1fs]" % (n, t, d, delta, p, q, m, W, kap, tF, 2 * delta - 1, "OK" if ok else "FAIL", time.time() - t0))
        sys.stdout.flush(); done += 1
