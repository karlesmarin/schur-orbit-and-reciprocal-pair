# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# dos_ctrl.sage -- item 2 near-miss control: same prediction with p | h^-(Q(zeta_n)) (hypothesis violated), -1 not in <2>.
# What: layers W_0..W_2 by inst.sage; the prediction [1, d-delta, d] is expected to be able to fail here.  Run: sh sg.sh dos_ctrl.sage
load("inst.sage")
for (n, p) in [(75, 11), (77, 5), (51, 5), (23, 3), (31, 3), (47, 5), (49, 43)]:
    d = euler_phi(n) // 2; t = Mod(2, n).multiplicative_order(); delta = d // t
    v = 1
    while euler_phi(p**v) < 3:
        v += 1
    q = p**v * n; m = (n - 1) // 2
    W = Img(q, p, m, 3).W(); rk, _ = moment_rank(n, m, p)
    print("CTRL n=%d p=%d (p|h^-) t=%d d=%d delta=%d q=%d m=%d | W=%s rank[S]=%d pred(if p!|h^-)=%s %s" % (n, p, t, d, delta, q, m, W, rk, [1, d - delta, d], "matches" if W == [1, d - delta, d] else "DIFFERS"))
    sys.stdout.flush()
