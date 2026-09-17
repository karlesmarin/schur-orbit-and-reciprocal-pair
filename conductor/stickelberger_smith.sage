# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# What: Smith form of the Stickelberger/moment matrix M_n vs exact h^-(Q(zeta_n)) (analytic class number formula), n in [3,90), n != 2 mod 4.
# Run (in this directory): sage stickelberger_smith.sage
import sys, time

def s_(x, n):
    x = x % n
    return 0 if x == 0 else 2*x - n

def M_of(n):
    reps = [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]
    rows = []
    for u in reps:
        ui = inverse_mod(u, n)
        rows.append([s_(ui*c, n) for c in range(1, n)])
    return matrix(ZZ, rows)

def hminus(n):
    G = DirichletGroup(n)
    K = G.base_ring()
    prod_ = K(1)
    nodd = 0
    check_sage = True
    for chi in G:
        if chi.is_odd():
            nodd += 1
            cp = chi.primitive_character()
            f = cp.modulus()
            B = sum(cp(a)*a for a in range(1, f+1)) / f
            if check_sage:
                B2 = cp.bernoulli(1)
                assert K(B2) == K(B), (n, chi)
            prod_ *= (-B/2)
    assert nodd == euler_phi(n)//2
    assert prod_ in QQ, n
    pr = QQ(prod_)
    Q = 1 if len(factor(n)) == 1 else 2
    w = 2*n if n % 2 == 1 else n
    h = Q * w * pr
    assert h in ZZ and h > 0, (n, h)
    return ZZ(h)

known = {23:3, 29:8, 31:9, 37:37, 39:2, 41:121, 43:211, 47:695, 49:43, 52:3, 53:4889, 56:2, 57:9, 59:41241,
         61:76301, 63:7, 64:17, 65:64, 67:853513, 68:8, 69:69, 71:3882809, 72:3, 73:11957417, 75:11, 76:19,
         77:1280, 79:100146415, 80:5, 83:838216959, 84:1, 85:6205, 87:1536, 88:55, 89:13379363737}
h1 = [3,4,5,7,8,9,11,12,13,15,16,17,19,20,21,24,25,27,28,32,33,35,36,40,44,45,48,60]
for m_ in h1:
    known[m_] = 1

odd_primes = [p for p in prime_range(3, 98)]
tot = 0; match = 0; mism = []
tot_semi = 0; match_semi = 0   # p | phi(n)
tot_ph = 0; match_ph = 0       # p | h^-
rank_checks = []; rank_other = []; semi_list = []
v2info = []; pn_info = []
known_res = []
t0 = time.time()
print("n phi d rank h^- Delta elemdiv_nontrivial")
for n in range(3, 90):
    if n % 4 == 2:
        continue
    phi = euler_phi(n); d = phi // 2
    M = M_of(n)
    assert M.nrows() == d and M.ncols() == n-1
    ed = M.elementary_divisors()
    nz = [e for e in ed if e != 0]
    rk = M.rank()
    Delta = prod(nz) if rk == d else 0
    h = hminus(n)
    if n in known:
        known_res.append((n, h, known[n], h == known[n]))
    nontriv = [e for e in nz if e != 1]
    print(n, phi, d, rk, h, Delta, "Delta/h^- =", (Delta/h if h else None), nontriv, "h^- fact:", factor(h), "Delta fact:", factor(Delta) if Delta else 0)
    sys.stdout.flush()
    for p in odd_primes:
        if n % p == 0:
            pn_info.append((n, p, Delta.valuation(p), h.valuation(p)))
            continue
        vD = Delta.valuation(p); vh = h.valuation(p)
        tot += 1
        ok = (vD == vh)
        match += ok
        if not ok:
            mism.append((n, p, vD, vh, phi % p == 0))
        if phi % p == 0:
            tot_semi += 1; match_semi += ok
            semi_list.append((n, p, vD, vh))
        if vh > 0:
            tot_ph += 1; match_ph += ok
            rF = matrix(GF(p), M).rank()
            if vh == 1:
                rank_checks.append((n, p, rF, d, rF == d-1))
            else:
                rank_other.append((n, p, vh, rF, d))
    v2info.append((n, Delta.valuation(2), h.valuation(2)))

print()
print("TOTAL pairs (n, odd p<=97, p!|n):", tot, "matches:", match, "mismatches:", len(mism))
print("  with p | phi(n):", tot_semi, "matches:", match_semi)
print("  with p | h^-   :", tot_ph, "matches:", match_ph)
print("MISMATCHES (n,p,vp(Delta),vp(h-),p|phi):", mism)
print("pairs p|phi(n) detail (n,p,vD,vh):", semi_list)
print("pairs p|phi(n) AND p|h^-:", [x for x in semi_list if x[3] > 0])
print("rank_Fp for vp(h-)>=2 (n,p,vp(h-),rank,d):", rank_other)
print("rank_Fp checks for vp(h-)=1 (n,p,rank,d,rank==d-1):")
for r in rank_checks: print("  ", r)
print("all rank==d-1:", all(r[4] for r in rank_checks), "count", len(rank_checks))
print("v2 info (n, v2(Delta), v2(h-)):", v2info)
print("p|n info (n,p,vp(Delta),vp(h-)):", pn_info)
print("KNOWN-values check (n, computed, table, eq):")
for r in known_res: print("  ", r)
print("known all equal:", all(r[3] for r in known_res), "count", len(known_res))
print("time", time.time()-t0)
