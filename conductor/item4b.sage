# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# item4b.sage -- search for e=2 cases with kappa > 1 (discriminating test of type = d - W1 + kappa - 1), plus A_40(451) at 11.
# e=2 certified by: W1>=1 and W2=d (Prop cert (b)), or E=2 and d/2<W1<d (Prop mono (b)).
# type by F_p colon (capas, valid since e=2<=E); lattice colon (lat.sage) when D <= DLAT.
# Run (in this directory): sage item4b.sage QMAX DLAT
import sys, time
load("lat.sage")
QMAX = int(sys.argv[1]); DLAT = int(sys.argv[2])

def residue_degree(p, n):
    x = p % n; f = 1
    while x not in (1, n - 1):
        x = (x * p) % n; f += 1
    return f

# A_40(451) at p = 11 first
q, p, m = 451, 11, 40
v, n = split_q(q, p); d = euler_phi(n) // 2; f = residue_degree(p, n)
C = Capas(q, p, m, 3)
W = C.W(); W1, kap = C.kappa(); tF = C.type_Fp(2) if (W[1] >= 1 and W[2] == d) else None
print("A_40(451) p=11: q'=%d d=%d E=%d f=%d g=%d | W=%s W1=%d kappa=%d type_Fp(e=2)=%s pred=%d  2m+1=%d" % (
    n, d, euler_phi(p**v), f, d // f, W, W1, kap, tF, d - W1 + kap - 1, 2 * m + 1))
sys.stdout.flush()

tot = 0; fails = 0; big = []
from collections import Counter
hs = Counter()
for q in range(5, QMAX + 1):
    for p in prime_divisors(q):
        v, n = split_q(q, p)
        if n in (1, 2, 3, 4, 6):
            continue
        E = euler_phi(p**v)
        if E < 2:
            continue
        d = euler_phi(n) // 2
        for m in range(1, (q + 1) // 2):
            if 2 * m >= q:
                continue
            if (2 * m) % n and (2 * m + 1) % n and (2 * m + 2) % n:
                continue
            if d * min(3, E) > 120:
                continue
            K = min(3, E)
            C = Capas(q, p, m, K)
            W = C.W()
            if W[1] >= d or W[1] == 0:
                continue
            if K >= 3:
                if W[2] != d:
                    continue
            else:
                if not (2 * W[1] > d):
                    continue
            W1, kap = C.kappa()
            tF = C.type_Fp(2)
            pred = d - W1 + kap - 1
            ok = (tF == pred)
            tot += 1; fails += (not ok)
            lOA = 2*d - 1 - W1; lAf = 1 + W1
            hs[("HS_ok" if lOA - lAf >= tF - 1 else "HS_VIOLATION", "almost_gor" if lOA - lAf == tF - 1 else "strict", "gor" if tF == 1 else "nongor", "p=2" if p == 2 else "p odd")] += 1
            if kap > 1 or not ok:
                g = d // residue_degree(p, n)
                line = "q=%d m=%d p=%d q'=%d d=%d g=%d E=%d W=%s kappa=%d type_Fp=%d pred=%d %s" % (q, m, p, n, d, g, E, W, kap, tF, pred, "OK" if ok else "FAIL")
                typC = None
                D = euler_phi(q) // 2
                if D <= DLAT:
                    res = analyse(m, q, p, colon=True)
                    typC = res.get("typeC")
                    line += " | lattice: e=%s type_colon=%s %s" % (res["e"], typC, "OK" if (res["e"] == 2 and typC == pred) else "FAIL")
                    if not (res["e"] == 2 and typC == pred):
                        fails += 1
                print("KAPPA>1 " + line if kap > 1 else line)
                sys.stdout.flush()
print("HS12.2.3 on e=2 cases (lengths via Thm loc(d)):", dict(hs))
print("ITEM4b TOTAL e=2 cases (q<=%d, d*min(3,E)<=120)=%d fails=%d" % (QMAX, tot, fails))
