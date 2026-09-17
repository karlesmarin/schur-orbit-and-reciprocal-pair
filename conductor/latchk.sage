# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# latchk.sage -- item 2 lattice checks with lat.sage: index, e (r^e in A_p by lattice sum), layers, F_p type, colon-lattice type.
# What: A_7(105) p=7 (n=15, t=4, -1 not in <2>, delta=1); optional A_25(357) p=7 (n=51, delta=2); A_31(315) p=5 (n=63, delta=3).
# Run: sh sg.sh latchk.sage m q p
load("lat.sage")
m, q, p = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
res = analyse(m, q, p, colon=True)
v, n = split_q(q, p); d = euler_phi(n) // 2
t = Mod(2, n).multiplicative_order(); delta = d // t
print("PRED (Cor dos, -1 not in <2>): W1=%d e=2 v_p[O:A]=%d v_p[A:f]=%d type=%d" % (d - delta, d + delta - 1, d - delta + 1, 2 * delta - 1))
print("GOT: W=%s e=%s v_p[O:A]=%d v_p[A:f]=%d typeF=%s typeC=%s" % (res["W"], res["e"], res["lOA"], res["lAf"], res.get("typeF"), res.get("typeC")))
