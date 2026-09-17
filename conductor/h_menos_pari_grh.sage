# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# What: independent h^- = h(Q(zeta_n)) / h(Q(zeta_n)^+) via PARI bnfinit (class_number), for 20<phi(n)<=32, proof=False (GRH), n != 2 mod 4, compared to the analytic formula.
# Run (in this directory): sage h_menos_pari_grh.sage
import sys, time
def hminus(n):
    G = DirichletGroup(n); K = G.base_ring(); pr = K(1)
    for chi in G:
        if chi.is_odd():
            cp = chi.primitive_character(); f = cp.modulus()
            pr *= -(sum(cp(a)*a for a in range(1, f+1)) / f)/2
    Q = 1 if len(factor(n)) == 1 else 2
    w = 2*n if n % 2 == 1 else n
    return ZZ(Q*w*QQ(pr))
allok = True
for n in range(3, 90):
    if n % 4 == 2 or (euler_phi(n) <= 20 or euler_phi(n) > 32): continue
    t0 = time.time()
    C = CyclotomicField(n)
    h = C.class_number(proof=False)
    Kp = C.maximal_totally_real_subfield()[0]
    hp = Kp.class_number(proof=False)
    ha = hminus(n)
    ok = (h % hp == 0 and h // hp == ha)
    allok = allok and ok
    print(n, euler_phi(n), "h=", h, "h+=", hp, "h/h+=", h/hp, "analytic h-=", ha, "OK" if ok else "MISMATCH", "%.1fs" % (time.time()-t0))
    sys.stdout.flush()
print("all OK:", allok)
