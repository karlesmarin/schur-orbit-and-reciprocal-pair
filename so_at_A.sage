# Problem 10.1 has three columns -- o, sp, so -- and the paper answers only two of them.
#
#   o :  Lemma 8.8,  o_nu(W,1,-1) = sp_nu(W).  The two letters COLLAPSE the character.
#   sp:  the branching expansion added in Section 10, sp_nu(W,1,-1) = sum c_{nu/gamma} sp_gamma(W).
#   so:  never examined.
#
# Same one-line question, same machinery: adjoining 1 and -1 is the ring map p_k -> p_k + c^k, so
# what follows is an identity in Lambda with no length hypothesis.
#
# Authors: Carles Marin + Claude (AI assistant).

Sym = SymmetricFunctions(QQ)
p, O, SP, SO = Sym.p(), Sym.o(), Sym.sp(), Sym.so()

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

print("=" * 78)
print("CONTROL -- the identity Section 8 rests on")
print("=" * 78)
ok = bad = 0
for size in range(0, 8):
    for lam in Partitions(size):
        lam = tuple(lam)
        f = O[list(lam)] if lam else O.one()
        g = SP[list(lam)] if lam else SP.one()
        if both(f) == p(g): ok += 1
        else: bad += 1
print("  o_nu(W,1,-1) = sp_nu(W): %d hold, %d fail   <-- must be 0" % (ok, bad))
assert bad == 0

print("")
print("=" * 78)
print("THE UNEXAMINED COLUMN -- so_nu(W,1,-1), expanded in each universal basis")
print("=" * 78)
single_so = single_sp = single_o = 0
tested = 0
for size in range(0, 7):
    for lam in Partitions(size):
        lam = tuple(lam)
        g = SO[list(lam)] if lam else SO.one()
        val = both(g)
        tested += 1
        a, b, c = SO(val), SP(val), O(val)
        if len(list(a)) == 1: single_so += 1
        if len(list(b)) == 1: single_sp += 1
        if len(list(c)) == 1: single_o += 1
        if size <= 4:
            print("  nu=%-11s so: %-26s sp: %-24s o: %s"
                  % (str(lam), str(a)[:26], str(b)[:24], str(c)[:26]))
print("")
print("  shapes tested: %d" % tested)
print("  single term in the so basis : %d" % single_so)
print("  single term in the sp basis : %d" % single_sp)
print("  single term in the o  basis : %d   <-- a collapse like Lemma 8.8 would give all of them"
      % single_o)
