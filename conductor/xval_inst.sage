# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# xval_inst.sage -- cross-validate inst.sage (y-ring, product tree) against capas.sage (X-ring, Dickson recursion).
# What: W profiles on small (m,q,p,K), grid incl. non-stable segments; timing on one large case.  Run: sh sg.sh xval_inst.sage
load("capas.sage")
load("inst.sage")
tot = 0; bad = 0
cases = []
for (p, n) in [(3,5),(3,7),(5,7),(3,11),(5,9),(7,5),(3,13),(5,13),(7,15),(3,20),(5,12)]:
    for v in [1,2,3]:
        q = p**v * n; E = euler_phi(p**v)
        if E < 2 or q > 700: continue
        for m in range(2, min(q//2, 60), 3):
            cases.append((q, p, m, min(E, 5)))
for (q,p,m,K) in cases:
    a = Capas(q, p, m, K).W()
    b = Img(q, p, m, K).W()
    tot += 1
    if a != b:
        bad += 1; print("MISMATCH", q, p, m, K, a, b)
print("xval cases=%d mismatches=%d" % (tot, bad)); sys.stdout.flush()
t0 = time.time()
p=3; a=3; n=31; M=13; P=27; v=6; q=p**v*n; m=(P*M*n-1)//2
I = Img(q, p, m, P+1)
print("big m=%d D=%d dimB=%d time %.1fs" % (m, I.D, I.B.nrows(), time.time()-t0))
print("W=", I.W())
