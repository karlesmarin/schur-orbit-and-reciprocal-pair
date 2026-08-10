# Authors: Carles Marin + Claude (AI assistant).
# Check EVERY formula asserted in the rewritten Sec. 8, literally as written,
# so no claim in the paper rests on a plausible-looking manipulation.

Sym = SymmetricFunctions(QQ)
p, h, O, SP = Sym.p(), Sym.h(), Sym.o(), Sym.sp()

def addletter(f, c):
    out = p.zero()
    for rho, coeff in p(f):
        term = p.one()
        for k in rho:
            term = term * (p[k] + c^k)
        out += coeff * term
    return out

print("="*74)
print("F1 -- the recurrence  h_m(W,c) - c h_{m-1}(W,c) = h_m(W)")
print("="*74)
bad = 0
for c in (1, -1, 2, -3):
    for m in range(0, 9):
        lhs = addletter(h[m] if m > 0 else h.one(), c)
        prev = addletter(h[m-1] if m-1 > 0 else h.one(), c) if m >= 1 else p.zero()
        if p(lhs - c*prev) != p(h[m] if m > 0 else h.one()):
            bad += 1
print("  c in {1,-1,2,-3}, m <= 8 : %d failures" % bad)
assert bad == 0
print("  and at c = -1 it reads  h_m(W,-1) + h_{m-1}(W,-1) = h_m(W),")
print("  which is [AK25, Proposition 2.1] verbatim.")

def hk(k):
    return h.one() if k == 0 else (h[k] if k > 0 else h.zero())

def detperm(E):
    from itertools import permutations
    n = len(E)
    if n == 0: return h.one()
    tot = h.zero()
    for perm in permutations(range(n)):
        s = 1
        for a in range(n):
            for b in range(a+1, n):
                if perm[a] > perm[b]: s = -s
        t = h.one()
        for i in range(n): t = t * E[i][perm[i]]
        tot += s*t
    return tot

print("")
print("="*74)
print("F2 -- the DEFINING determinants, as the paper cites them from [AK25]")
print("="*74)
def o_det(lam):
    n = len(lam)
    if n == 0: return h.one()
    return detperm([[hk(lam[i]-(i+1)+(j+1)) - hk(lam[i]-(i+1)-(j+1))
                     for j in range(n)] for i in range(n)])
def so_det(lam):
    n = len(lam)
    if n == 0: return h.one()
    return detperm([[hk(lam[i]-(i+1)+(j+1)) + hk(lam[i]-(i+1)-(j+1)+1)
                     for j in range(n)] for i in range(n)])
def sp_det(lam):
    n = len(lam)
    if n == 0: return h.one()
    return detperm([[hk(lam[i]-(i+1)+(j+1)) + hk(lam[i]-(i+1)-(j+1)+2)
                     for j in range(n)] for i in range(n)]) / 2
b1 = b2 = 0
for size in range(0, 7):
    for lam in Partitions(size):
        lam = list(lam)
        if p(o_det(lam)) != p(O[lam] if lam else O.one()): b1 += 1
        if p(sp_det(lam)) != p(SP[lam] if lam else SP.one()): b2 += 1
print("  o_lambda  determinant vs Sage's o basis  : %d failures" % b1)
print("  sp_lambda determinant vs Sage's sp basis : %d failures" % b2)
assert b1 == 0 and b2 == 0
print("  (so_lambda has no Sage basis; it is pinned by F3 below.)")

print("")
print("="*74)
print("F3 -- THE COLUMN OPERATION, exactly as the proof states it")
print("="*74)
print("  Claim to test: the entries of the so-determinant are obtained from those")
print("  of the o-determinant evaluated at (W,1) by  C_j -> C_j - C_{j-1}.")
print("  A column operation of that shape leaves the determinant unchanged, so if")
print("  the entrywise claim holds, so does  iota_1 o_nu = so_nu.")
print("")
bad3 = badsign = 0
for size in range(0, 7):
    for lam in Partitions(size):
        lam = list(lam)
        n = len(lam)
        if n == 0: continue
        # E = entries of o_lambda evaluated at (W, 1)
        E = [[addletter(hk(lam[i]-(i+1)+(j+1)) - hk(lam[i]-(i+1)-(j+1)), 1)
              for j in range(n)] for i in range(n)]
        # S = entries of so_lambda over W
        S = [[p(hk(lam[i]-(i+1)+(j+1)) + hk(lam[i]-(i+1)-(j+1)+1))
              for j in range(n)] for i in range(n)]
        for i in range(n):
            if p(E[i][0]) != S[i][0]:
                bad3 += 1
            for j in range(1, n):
                if p(E[i][j]) - p(E[i][j-1]) != S[i][j]:
                    bad3 += 1
                # the WRONG sign, which the first draft of the proof wrote
                if p(E[i][j]) + p(E[i][j-1]) == S[i][j]:
                    badsign += 1
print("  C_j -> C_j - C_{j-1} reproduces every so entry : %d failures" % bad3)
print("  cases where C_j + C_{j-1} would also work      : %d" % badsign)
assert bad3 == 0
print("  >> the sign is MINUS at c = 1.  [AK25, Lemma 3.2] writes PLUS because its")
print("     step is c = -1; the general operation is C_j -> C_j - c C_{j-1}.")

print("")
print("="*74)
print("F4 -- and therefore the two steps, then the composition")
print("="*74)
s1 = s2 = s3 = 0
for size in range(0, 8):
    for lam in Partitions(size):
        lam = list(lam)
        f = O[lam] if lam else O.one()
        if p(addletter(f, 1)) != p(so_det(lam)): s1 += 1
        if p(addletter(so_det(lam), -1)) != p(SP[lam] if lam else SP.one()): s2 += 1
        if p(addletter(addletter(f, 1), -1)) != p(SP[lam] if lam else SP.one()): s3 += 1
print("  iota_1 o_nu = so_nu           , |nu| <= 7 : %d failures" % s1)
print("  iota_{-1} so_nu = sp_nu       , |nu| <= 7 : %d failures" % s2)
print("  iota_{-1} iota_1 o_nu = sp_nu , |nu| <= 7 : %d failures" % s3)
assert s1 == s2 == s3 == 0

print("")
print("="*74)
print("F5 -- nu'_1 = l(nu), so the trichotomy is on LENGTH against the rank")
print("="*74)
bad5 = 0
for size in range(0, 9):
    for nu in Partitions(size):
        if len(nu) and Partition(list(nu)).conjugate()[0] != len(nu):
            bad5 += 1
print("  nu'_1 = l(nu) for every partition, |nu| <= 8 : %d failures" % bad5)
assert bad5 == 0
print("  N/2 = (2r+2)/2 = r+1, so  nu'_1 = N/2  <=>  l(nu) = r+1  and")
print("                            nu'_1 > N/2  <=>  l(nu) > r+1.")
print("DONE")
