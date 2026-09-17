# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# item3.sage -- factor at 2 without semisimplicity: n prime <= 61, p odd <= 31, p !| n h^-, 2m+1 = N n, p !| N.
# Claim: d - W_1 = 0 (t even) or d/t (t odd), t = ord_n(2); and W_2 = d (so e = 1 or 2).
# Run (in this directory): sage item3.sage
import sys, time
load("capas.sage")
hminus = {}
for ln in open("stickelberger_smith_OUT.txt"):
    f = ln.split()
    if len(f) > 5 and f[0].isdigit():
        hminus[int(f[0])] = int(f[4])
tot = 0; fail = 0; tot_div = 0; fail_div = 0
for n in prime_range(5, 62):
    h = hminus[n]; d = (n - 1) // 2
    t = Mod(2, n).multiplicative_order()
    delta = 0 if t % 2 == 0 else d // t
    for p in prime_range(3, 32):
        if p == n or h % p == 0:
            continue
        Ns = [N for N in range(1, 40, 2) if N % p != 0][:3]
        for N in Ns:
            m = (N * n - 1) // 2
            v = 1
            while p**v <= N or euler_phi(p**v) < 3:
                v += 1
            q = p**v * n
            E = euler_phi(p**v); K = min(3, E)
            t0 = time.time()
            C = Capas(q, p, m, K)
            W = C.W()
            ok1 = (d - W[1] == delta)
            ok2 = (W[2] == d)
            ok = ok1 and ok2
            div = ((n - 1) % p == 0)
            tot += 1; fail += (not ok)
            tot_div += div; fail_div += (div and not ok)
            if div or not ok or N == Ns[0]:
                print("n=%2d p=%2d N=%2d v=%d q=%d m=%d d=%d t=%d delta=%d p|n-1=%s | W=%s | d-W1=%d pred=%d W2full=%s %s  [%.1fs]" % (
                    n, p, N, v, q, m, d, t, delta, div, W, d - W[1], delta, ok2, "OK" if ok else "FAIL", time.time() - t0))
                sys.stdout.flush()
print("TOTAL cases=%d fails=%d ; with p|n-1: cases=%d fails=%d" % (tot, fail, tot_div, fail_div))
