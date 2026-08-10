# Authors: Carles Marin + Claude (AI assistant).
# CORRECTIONS.md speculates: if no non-standard label has o_nu(A) != 0, then in the
# unstable range only standard labels survive, Littlewood's rule applies verbatim,
# and conj:crit becomes a THEOREM for r >= 3.  Lemma lem:AtoSp makes this decidable.
Sym = SymmetricFunctions(QQ)
p, O = Sym.p(), Sym.o()

def o_at_A(nu, r, R, zs):
    out = R(0)
    for rho, c in p(O[list(nu)] if nu else O.one()):
        t = R(c)
        for k in rho:
            t *= R(1 + (-1)^k + sum(zs[j]^k + zs[j]^(-k) for j in range(r)))
        out += t
    return out

print("Does EVERY non-standard label vanish on A?  (nu'_1 + nu'_2 > N = 2r+2)")
print("%-4s %-8s %-10s %-10s %-10s %s" % ("r","N","|nu|<=","non-std","vanish","NONZERO"))
for r in (1, 2, 3):
    N = 2*r+2
    R = LaurentPolynomialRing(QQ, ['z%d'%(j+1) for j in range(r)])
    zs = list(R.gens())
    MAX = 14 if r == 1 else (12 if r == 2 else 11)
    ns = van = nz = 0
    examples = []
    for size in range(0, MAX+1):
        for nu in Partitions(size):
            nu = tuple(nu)
            c = list(Partition(list(nu)).conjugate()) if nu else [0]
            c1 = c[0]; c2 = c[1] if len(c) > 1 else 0
            if c1 + c2 <= N:
                continue
            ns += 1
            v = o_at_A(nu, r, R, zs)
            if v == 0: van += 1
            else:
                nz += 1
                if len(examples) < 3: examples.append(nu)
    print("%-4d %-8d %-10d %-10d %-10d %d  %s" % (r, N, MAX, ns, van, nz, examples))
print("")
print("If the NONZERO column is 0 for r>=3, CORRECTIONS.md's target is real.")
print("If it is not, the non-standard labels do NOT all die; they fold, which is")
print("what Lemma lem:AtoSp says they do.")
