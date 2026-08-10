# Authors: Carles Marin + Claude (AI assistant).
# Problem 10.1, the case the paper singles out: sp_lambda at the alphabet (1,-1,z,zbar).
#
# Section 8 rests on   o_nu(W,1,-1) = sp_nu(W)   [AK25 (2.13) + Lemma 3.2 composed],
# so at r = 1 the orthogonal universal character on A = (1,-1,z,zbar) is a RANK-ONE
# symplectic character, and its law follows.  The symplectic column of Problem 10.1
# has no law recorded.  The companion question is therefore
#
#         what is   sp_nu(W, 1, -1)  ?
#
# Adding the letter c is the ring map p_k -> p_k + c^k, so this is an identity in
# Lambda and can be tested for every partition with no length hypothesis.

Sym = SymmetricFunctions(QQ)
p, h, e, O, SP, SO = Sym.p(), Sym.h(), Sym.e(), Sym.o(), Sym.sp(), Sym.s()

def addletter(f, c):
    out = p.zero()
    for rho, coeff in p(f):
        term = p.one()
        for k in rho:
            term = term * (p[k] + c^k)
        out += coeff * term
    return out

def both(f):
    return p(addletter(addletter(f, 1), -1))

print("="*78)
print("CONTROL -- the identity Section 8 already rests on must come out")
print("="*78)
ok = bad = 0
for size in range(0, 8):
    for lam in Partitions(size):
        lam = tuple(lam)
        f = O[list(lam)] if lam else O.one()
        g = SP[list(lam)] if lam else SP.one()
        if both(f) == p(g): ok += 1
        else: bad += 1
print("  o_nu(W,1,-1) = sp_nu(W):  %d hold, %d fail   <-- must be 0 fail" % (ok, bad))
assert bad == 0

print("")
print("="*78)
print("THE QUESTION -- expand sp_nu(W,1,-1) in the three universal bases")
print("="*78)
for size in range(0, 7):
    for lam in Partitions(size):
        lam = tuple(lam)
        g = SP[list(lam)] if lam else SP.one()
        val = both(g)
        insp = SP(val)
        ino = O(val)
        print("  nu=%-12s  in sp: %-34s  in o: %s"
              % (str(lam), str(insp)[:34], str(ino)[:40]))
