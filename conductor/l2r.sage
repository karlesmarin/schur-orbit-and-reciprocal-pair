# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# l2r.sage -- item 3: column lattice of M_{2r} equals 2 * column lattice of M_r (odd r), compatible odd row representatives;
#             subfamily function h_{2r} spans L_r; end-to-end W_1 for q'=2r families and subfamilies.
# What: ZZ column modules (HNF), Smith forms, Delta; layers by inst.sage (K=2).   Run: sh sg.sh l2r.sage
load("inst.sage")
hm = {}
for ln in open("stickelberger_smith_OUT.txt"):
    t = ln.split()
    if t and t[0].isdigit():
        hm[int(t[0])] = Integer(t[4])

def s(n, c):
    c %= n
    return 0 if c == 0 else 2 * c - n

def reps(r):
    return [u for u in range(1, r) if gcd(u, r) == 1 and u < r - u]

def Mr(r):
    return matrix(ZZ, [[s(r, inverse_mod(u, r) * c) for c in range(1, r)] for u in reps(r)])

def odd_rep(u, r):
    return u if u % 2 else u + r

def M2r(r):
    n = 2 * r
    return matrix(ZZ, [[s(n, inverse_mod(odd_rep(u, r), n) * c) for c in range(1, n)] for u in reps(r)])

def h2r(r, c):
    c %= 2 * r
    if c == 0 or c == r: return 0
    return c if c < r else c - 2 * r

def H2r(r):
    n = 2 * r
    return matrix(ZZ, [[h2r(r, inverse_mod(odd_rep(u, r), n) * c) for c in range(1, n)] for u in reps(r)])

def colmod(M):
    return M.transpose().row_module()

bad = 0; tot = 0
for r in range(3, 46, 2):
    d = euler_phi(r) // 2
    A = Mr(r); B = M2r(r); C = H2r(r)
    Lr = colmod(A); L2r = colmod(B); Lh = colmod(C)
    L2 = colmod(2 * A)
    eq = (L2r == L2); eqh = (Lh == Lr)
    edA = [x for x in A.elementary_divisors() if x != 0]; edB = [x for x in B.elementary_divisors() if x != 0]
    ed_ok = (edB == [2 * x for x in edA])
    DA = prod(edA) if len(edA) == d else 0; DB = prod(edB) if len(edB) == d else 0
    del_ok = (DB == 2**d * DA)
    tot += 1; bad += not (eq and eqh and ed_ok and del_ok)
    print("r=%2d d=%2d rank(M_r)=%d rank(M_2r)=%d | L_2r==2L_r %s | L(h_2r)==L_r %s | ed(M_2r)==2ed(M_r) %s | Delta_2r=2^d Delta_r %s (Delta_r=%s, h^-(r)=%s)" % (
        r, d, A.rank(), B.rank(), eq, eqh, ed_ok, del_ok, factor(DA) if DA else 0, hm.get(r)))
    sys.stdout.flush()
print("LATTICE TOTAL r=%d fail=%d" % (tot, bad))

# end-to-end
tot = bad = 0
for r in [5, 7, 11, 13, 23, 29, 31]:
    qp = 2 * r; d = euler_phi(r) // 2
    for p in [3, 5, 7]:
        if qp % p == 0: continue
        rkMr = Mr(r).change_ring(GF(p)).rank()
        cases = []
        for N in [1, 2, 3, 4, 5, 7]:
            if N % p: cases += [("m=Nq'", N, N * qp), ("m=Nq'-1", N, N * qp - 1)]
        for Bo in [1, 3, 5, 7]:
            if Bo % p: cases += [("m=Br", Bo, Bo * r), ("m=Br-1", Bo, Bo * r - 1)]
        for (fam, N, m) in cases:
            v = 1
            while 2 * m >= p**v * qp or euler_phi(p**v) < 2:
                v += 1
            q = p**v * qp
            I = Img(q, p, m, 2); W = I.W()
            rkS, _ = moment_rank(qp, m, p)
            ok = (W[1] == rkMr == rkS)
            if fam in ("m=Nq'", "m=Nq'-1"):
                ok = ok and ((W[1] == d) == (hm[r] % p != 0))
            tot += 1; bad += (not ok)
            print("r=%2d q'=%d p=%d %-8s N/B=%d m=%d v=%d q=%d | W=%s rank_Fp(M_r)=%d rank[S]=%d d=%d p|h^-=%s %s" % (
                r, qp, p, fam, N, m, v, q, W, rkMr, rkS, d, hm[r] % p == 0, "OK" if ok else "FAIL"))
            sys.stdout.flush()
print("END-TO-END TOTAL=%d fail=%d" % (tot, bad))
