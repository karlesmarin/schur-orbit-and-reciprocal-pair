# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# item3_ctrl.sage -- near-miss control for item 3: same formula on pairs with p | h^-(Q(zeta_n)) (hypothesis violated).
# Run (in this directory): sage item3_ctrl.sage
load("capas.sage")
hminus = {}
for ln in open("stickelberger_smith_OUT.txt"):
    f = ln.split()
    if len(f) > 5 and f[0].isdigit():
        hminus[int(f[0])] = int(f[4])
for n in prime_range(5, 62):
    h = hminus[n]; d = (n - 1) // 2
    t = Mod(2, n).multiplicative_order(); delta = 0 if t % 2 == 0 else d // t
    for p in prime_range(3, 32):
        if p == n or h % p != 0:
            continue
        v = 2 if p == 3 else 1
        q = p**v * n; m = (n - 1) // 2
        W = Capas(q, p, m, 3).W()
        print("CTRL n=%d p=%d v_p(h^-)=%d t=%d delta=%d | W=%s d-W1=%d pred=%d %s" % (n, p, valuation(h, p), t, delta, W, d - W[1], delta, "holds" if d - W[1] == delta else "DEVIATES"))
