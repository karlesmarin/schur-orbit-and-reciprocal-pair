# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# medio.sage -- half-period families m = Br, Br-1 (B odd, r=q'/2) when 4 | q'.
# What: (A) shift identity hat s_n(c) = s_n(c+n/2)/2 for even n, hence L'_n = (1/2)L_n and equal odd
#       elementary divisors; (B) end-to-end W_1 by inst.sage vs rank mod p of M_{q'} and vs p | h^-;
#       (C) controls with p | B.   Run: sh sg.sh medio.sage
load("inst.sage")

hm = {}
for ln in open("stickelberger_smith_OUT.txt"):
    t = ln.split()
    if t and t[0].isdigit():
        hm[int(t[0])] = Integer(t[4])

def s(n, c):
    c %= n
    return 0 if c == 0 else 2 * c - n

def hs(n, c):          # signed residue, hat s_n; hat s_n(0) = hat s_n(n/2) = 0
    c %= n
    if 2 * c == n: return 0
    return c if 2 * c < n else c - n

def reps(n):
    return [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]

def Mfull(n, f):
    R = reps(n)
    return matrix(ZZ, [[f(n, inverse_mod(u, n) * c) for c in range(n)] for u in R])

def colmod(M):
    return M.transpose().row_module()

def oddpart(x):
    x = Integer(x)
    return 0 if x == 0 else x // 2**valuation(x, 2)

print("=== (A) shift identity, lattices and elementary divisors, even n ===")
tot = bad = 0
for n in range(4, 71, 2):
    d = euler_phi(n) // 2
    M = Mfull(n, s); H = Mfull(n, hs)
    sh = all(H.column(c) * 2 == M.column((c + n // 2) % n) for c in range(n))
    L = colmod(M); Lh = colmod(H)
    lat = (L == colmod(2 * H))
    edM = [Integer(x) for x in M.elementary_divisors() if x != 0]
    edH = [Integer(x) for x in H.elementary_divisors() if x != 0]
    ed = (edM == [2 * x for x in edH])
    odd = ([oddpart(x) for x in edM] == [oddpart(x) for x in edH])
    rk = all(M.change_ring(GF(p)).rank() == H.change_ring(GF(p)).rank()
             for p in [3, 5, 7, 11, 13] if n % p)
    ok = sh and lat and ed and odd and rk and len(edM) == d
    tot += 1; bad += (not ok)
    print("n=%2d (n mod 4 = %d) d=%2d | hat=shift/2 %s | L=2L' %s | ed(M)=2ed(Mhat) %s | odd parts equal %s | rank_p equal %s | %s"
          % (n, n % 4, d, sh, lat, ed, odd, rk, "OK" if ok else "FAIL"))
    sys.stdout.flush()
print("IDENTITY TOTAL=%d fail=%d" % (tot, bad))

print("")
print("=== (B) end-to-end: 4 | q', m = Br and Br-1 with B odd, p odd, p not dividing B ===")
QS = [8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60]
BIG = [(52, 3), (64, 17), (88, 5)]
cases = []
for qp in QS:
    for p in [3, 5, 7]:
        if qp % p == 0: continue
        for B in [1, 3, 5]:
            if B % p == 0: continue
            cases.append((qp, p, B))
for (qp, p) in BIG:
    cases.append((qp, p, 1))

tot = bad = 0
for (qp, p, B) in cases:
    r = qp // 2; d = euler_phi(qp) // 2
    Mq = Mfull(qp, s)
    edM = [Integer(x) for x in Mq.elementary_divisors() if x != 0]
    ndiv = len([x for x in edM if x % p == 0])
    rkM = Mq.change_ring(GF(p)).rank()
    rkH = Mfull(qp, hs).change_ring(GF(p)).rank()
    for m in [B * r, B * r - 1]:
        if m < 2: continue
        v = 1
        while 2 * m >= p**v * qp or euler_phi(p**v) < 2:
            v += 1
        q = p**v * qp
        t0 = walltime()
        I = Img(q, p, m, 2); W = I.W()
        rkS, _ = moment_rank(qp, m, p)
        Sv = [0] * qp
        for j in range(1, m + 1):
            Sv[j % qp] += j; Sv[(-j) % qp] -= j
        Bo = ceil(m / r)
        sform = all(Sv[c] == Bo * hs(qp, c) for c in range(qp)) and Bo % 2 == 1
        ok = sform and (W[1] == rkM == rkH == rkS == d - ndiv) and ((W[1] == d) == (hm[qp] % p != 0))
        tot += 1; bad += (not ok)
        print("q'=%2d p=%2d B=%d m=%4d v=%d q=%6d | S=B*hat_s %s | W=%s | rank_p M=%d Mhat=%d [S]=%d | d=%2d #ed(p)=%d h^-=%s p|h^- %s | %.1fs %s"
              % (qp, p, B, m, v, q, sform, W, rkM, rkH, rkS, d, ndiv, hm[qp], hm[qp] % p == 0, walltime(t0), "OK" if ok else "FAIL"))
        sys.stdout.flush()
print("END-TO-END TOTAL=%d fail=%d" % (tot, bad))

print("")
print("=== (C) controls: p | B (the hypothesis p does not divide ceil(m/r)) ===")
for (qp, p, B) in [(8, 3, 3), (20, 3, 3), (28, 3, 3), (52, 3, 3), (16, 5, 5), (24, 5, 5)]:
    r = qp // 2; d = euler_phi(qp) // 2
    rkM = Mfull(qp, s).change_ring(GF(p)).rank()
    for m in [B * r, B * r - 1]:
        v = 1
        while 2 * m >= p**v * qp or euler_phi(p**v) < 2:
            v += 1
        q = p**v * qp
        I = Img(q, p, m, 2); W = I.W()
        rkS, _ = moment_rank(qp, m, p)
        print("q'=%2d p=%d B=%d m=%4d v=%d q=%6d | W=%s rank[S]=%d rank_p M=%d d=%2d | prediction %s"
              % (qp, p, B, m, v, q, W, rkS, rkM, d, "holds" if W[1] == rkM else "FAILS (as expected)"))
        sys.stdout.flush()
