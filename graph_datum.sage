# Authors: Carles Marin + Claude (AI assistant).
# WHICH DATUM IS MISSING AT EVEN t?
# Candidate: my signature used d_3 = |a1+a2-b1-b2|, the ABSOLUTE value.  The sign
# formula (eq. sign) carries sgn(a1+a2-b1-b2) as a separate factor, and Question
# 3.7 states that eps_lambda * sgn(sigma) * sgn(a1+a2-b1-b2) is constant on each
# cell (t, r_A, r_B).  So the missing datum should be the ORIENTATION of d_3.
Sym = SymmetricFunctions(QQ); p, s = Sym.p(), Sym.s()

def sgn_sigma(beta, t):
    order = sorted(range(len(beta)), key=lambda j: (beta[j] % t, -beta[j]))
    sg = 1
    for a in range(len(order)):
        for b in range(a+1, len(order)):
            if order[a] > order[b]: sg = -sg
    return sg

def sig(lam, t, mode):
    N = t+2
    lam = list(lam) + [0]*(N-len(lam))
    beta = [lam[i]+N-1-i for i in range(N)]
    prof = [0]*t; cls = {}
    for b in beta:
        prof[b % t] += 1; cls.setdefault(b % t, []).append(b)
    if 0 in prof: return None
    big = sorted([r for r in range(t) if prof[r] >= 2])
    if len(big) == 2:
        rA, rB = big
        A = sorted(cls[rA], reverse=True); B = sorted(cls[rB], reverse=True)
    else:
        rA = rB = big[0]; P = sorted(cls[rA], reverse=True)
        A = [P[0], P[1]]; B = [P[1], P[2]]
    d1 = A[0]-A[1]; d2 = B[0]-B[1]; e = A[0]+A[1]-B[0]-B[1]
    base = (tuple(prof), rA, rB, d1, d2, abs(e), sgn_sigma(beta, t))
    if mode == "plain":    return base
    if mode == "oriented": return base + (sign(e) if e != 0 else 0,)
    if mode == "a1b1":     return base + (sign(A[0]-B[0]),)
    if mode == "both":     return base + (sign(e) if e != 0 else 0, sign(A[0]-B[0]))

for t in (2, 3, 4, 5, 6):
    N = t+2
    K = CyclotomicField(t) if t > 2 else QQ
    zt = [K.zeta(t)^j for j in range(t)] if t > 2 else [K(1), K(-1)]
    R = LaurentPolynomialRing(K, ['z']); z = R.gen()
    def phi(lam):
        out = R(0)
        for rho, c in p(s[list(lam)] if lam else s.one()):
            term = R(c)
            for k in rho: term *= R(sum(K(w)^k for w in zt) + z^k + z^(-k))
            out += term
        return out
    MAX = 13 if t <= 3 else (11 if t <= 4 else 9)
    shapes = [(list(l), phi(list(l))) for size in range(0, MAX+1)
              for l in Partitions(size, max_length=N)]
    line = "  t=%d :" % t
    for mode in ("plain", "oriented", "a1b1", "both"):
        bk = {}
        for lam, v in shapes:
            d = sig(lam, t, mode)
            if d is None: continue
            bk.setdefault(d, set()).add(v)
        conf = sum(1 for vs in bk.values() if len(vs) > 1)
        line += "  %-9s %2d" % (mode, conf)
    print(line)
print("")
print("  columns: conflicts remaining under each signature.")
print("  'oriented' = the graph datum with sgn(a1+a2-b1-b2) added, i.e. d_3 as an")
print("  ORIENTED distance rather than a length.")
print("DONE")
