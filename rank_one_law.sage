# Authors: Carles Marin + Claude (AI assistant).
# Does the paper's own Lemma 8.8 answer part of its own open Problem 10.1
# ("the other classical types")?  That problem reports, as a measurement with no
# law, that o_lambda vanishes on 83.8% of the l(lambda)>=2 shapes at t=2.
# At t=2 the alphabet is {1,-1,z,1/z} = our A at r=1, and Lemma 8.8 says
# o_nu(A) = sp_nu(z), a rank-ONE symplectic character.  So o_nu(A) must vanish
# exactly when sp_nu dies in rank 1 -- which is a LAW, not a percentage.
Sym = SymmetricFunctions(QQ)
p, O, SP = Sym.p(), Sym.o(), Sym.sp()
R = LaurentPolynomialRing(QQ, ['z']); z = R.gen()

def o_at_A(nu):
    out = R(0)
    for rho, c in p(O[list(nu)] if nu else O.one()):
        t = R(c)
        for k in rho: t *= R(1 + (-1)^k + z^k + z^(-k))
        out += t
    return out
def sp_rank1(nu):
    out = R(0)
    for rho, c in p(SP[list(nu)] if nu else SP.one()):
        t = R(c)
        for k in rho: t *= R(z^k + z^(-k))
        out += t
    return out

print("="*76)
print("t = 2, i.e. r = 1: is o_nu(1,-1,z,1/z) = sp_nu(z)?  and when is it 0?")
print("="*76)
tot = zero = 0
bad = 0
by_len = {}
for size in range(0, 13):
    for nu in Partitions(size):
        nu = tuple(nu)
        v = o_at_A(nu)
        if v != sp_rank1(nu): bad += 1
        tot += 1
        L = len(nu)
        by_len.setdefault(L, [0,0])
        by_len[L][0] += 1
        if v == 0:
            zero += 1
            by_len[L][1] += 1
print("  identity o_nu(A) = sp_nu(z) : %d shapes, %d mismatches" % (tot, bad))
assert bad == 0
print("")
print("  %-10s %-10s %-10s %s" % ("l(nu)", "shapes", "vanish", "all vanish?"))
for L in sorted(by_len):
    a, b = by_len[L]
    print("  %-10d %-10d %-10d %s" % (L, a, b, "YES" if a == b else ("no" if b == 0 else "partly")))
ge2 = sum(a for L,(a,b) in by_len.items() if L >= 2)
ge2z = sum(b for L,(a,b) in by_len.items() if L >= 2)
print("")
print("  l(nu) >= 2 : %d shapes, %d vanish = %.1f%%" % (ge2, ge2z, 100.0*ge2z/ge2))
print("")
print("  >> at rank one the symplectic basis is {sp_(m)}, one row only.  So")
print("     sp_nu(z) = 0 for l(nu) = 2 and folds for l(nu) > 2; the vanishing of")
print("     o_nu on this alphabet is therefore governed by a LAW -- Lemma 8.8 plus")
print("     the rank-one modification -- not by a percentage.")
print("DONE")
