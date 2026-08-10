# Authors: Carles Marin + Claude (AI assistant).
# Section 8 ("Outside Littlewood's range"), derived rather than verified.
#
# The three facts that section recorded as computed-over-a-range are consequences
# of one identity, which is published:
#
#   [AK25, (2.13)]   so_lambda(X) = o_lambda(X, Xbar, 1, 0, 0, ...)   (Koike-Terada)
#   [AK25, Lemma 3.2] so_lambda(X, Xbar, -1) = sp_lambda(X)
#   ---------------------------------------------------------------
#   composing:       o_nu(A) = sp_nu(z_1, ..., z_r),  A as in (psir)
#
# and of the SYMPLECTIC modification rules, which are one of the two classical
# cases Fauser-Jarvis-King-Wybourne exempt.  Everything below is exact.
#
# Adding the letter c to an alphabet is the ring map p_k -> p_k + c^k on Lambda,
# so the identity can be tested for EVERY partition, with no length hypothesis.

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

print("="*76)
print("STEP 1 -- the identity, as an identity in Lambda (no length hypothesis)")
print("="*76)
ok = bad = 0
for size in range(0, 9):
    for lam in Partitions(size):
        lam = tuple(lam)
        f = O[list(lam)] if lam else O.one()
        lhs = p(addletter(addletter(f, 1), -1))
        rhs = p(SP[list(lam)] if lam else SP.one())
        if lhs == rhs: ok += 1
        else: bad += 1
print("  o_nu(W,1,-1) = sp_nu(W)  for |nu| <= 8 :  %d hold, %d fail" % (ok, bad))
assert bad == 0
print("  CONTROL: each letter alone")
c1 = sum(1 for size in range(0,6) for lam in Partitions(size)
         if p(addletter(O[list(lam)] if lam else O.one(), 1))
            == p(SP[list(lam)] if lam else SP.one()))
c2 = sum(1 for size in range(0,6) for lam in Partitions(size)
         if p(addletter(O[list(lam)] if lam else O.one(), -1))
            == p(SP[list(lam)] if lam else SP.one()))
tot = sum(1 for size in range(0,6) for lam in Partitions(size))
print("    o(W,1)=sp(W) in %d of %d;  o(W,-1)=sp(W) in %d of %d" % (c1, tot, c2, tot))
assert c1 == 1 and c2 == 1
print("    >> only the empty shape.  BOTH letters are needed -- the reciprocal pair.")

def sp_at_torus(nu, r, R, zs):
    """sp_nu evaluated on the Sp(2r) torus {z_j, 1/z_j}."""
    out = R(0)
    for rho, c in p(SP[list(nu)] if nu else SP.one()):
        t = R(c)
        for k in rho:
            t *= R(sum(zs[j]^k + zs[j]^(-k) for j in range(r)))
        out += t
    return out

def o_at_A(nu, r, R, zs):
    """o_nu evaluated at A = {1,-1,z_1,1/z_1,...}."""
    out = R(0)
    for rho, c in p(O[list(nu)] if nu else O.one()):
        t = R(c)
        for k in rho:
            t *= R(1 + (-1)^k + sum(zs[j]^k + zs[j]^(-k) for j in range(r)))
        out += t
    return out

print("")
print("="*76)
print("STEP 2 -- the three facts of Sec. 8, now as symplectic statements")
print("="*76)
print("  nu'_1 = l(nu), so the trichotomy on nu'_1 vs N/2 = r+1 is a trichotomy")
print("  on the LENGTH of nu against the rank r of Sp(2r):")
print("")
for r in (1, 2, 3):
    N = 2*r + 2
    R = LaurentPolynomialRing(QQ, ['z%d' % (j+1) for j in range(r)])
    zs = list(R.gens())
    MAX = 10 if r == 1 else (9 if r == 2 else 8)
    n_short = n_crit = n_long = 0
    crit_nonzero = long_bad = 0
    basis = {}
    for size in range(0, MAX+1):
        for nu in Partitions(size):
            nu = tuple(nu)
            v = o_at_A(nu, r, R, zs)
            w = sp_at_torus(nu, r, R, zs)
            assert v == w, ("identity broke", nu)
            L = len(nu)
            if L <= r:
                n_short += 1
                basis[nu] = v
            elif L == r + 1:
                n_crit += 1
                if v != 0:
                    crit_nonzero += 1
            else:
                n_long += 1
                if v != 0 and not any(v == u or v == -u for u in basis.values()):
                    long_bad += 1
    print("  r=%d (N=%d), |nu| <= %d" % (r, N, MAX))
    print("     l(nu) <= r      : %4d  -- irreducible Sp(2r) characters, the basis of (Cmu)"
          % n_short)
    print("     l(nu) = r+1     : %4d  -- %d nonzero   [FACT 1: self-associate labels vanish]"
          % (n_crit, crit_nonzero))
    print("     l(nu) >  r+1    : %4d  -- %d not +- a basis element  [FACTS 2 and 3]"
          % (n_long, long_bad))
    assert crit_nonzero == 0 and long_bad == 0

print("")
print("="*76)
print("STEP 3 -- and the linear independence of (Cmu) is now a theorem, not an")
print("          observation: distinct irreducible characters of Sp(2r).")
print("="*76)
for r in (1, 2, 3):
    R = LaurentPolynomialRing(QQ, ['z%d' % (j+1) for j in range(r)])
    zs = list(R.gens())
    vals = {}
    dup = 0
    for size in range(0, 9):
        for nu in Partitions(size):
            nu = tuple(nu)
            if len(nu) > r: continue
            v = sp_at_torus(nu, r, R, zs)
            if v in vals.values(): dup += 1
            vals[nu] = v
    M = Matrix(QQ, [[QQ(c) for c in
                     [vals[k].dict().get(m, 0) for m in
                      sorted(set().union(*[set(u.dict()) for u in vals.values()]))]]
                    for k in sorted(vals)])
    print("  r=%d : %d characters, rank of their coefficient matrix = %d, repeats = %d"
          % (r, len(vals), M.rank(), dup))
    assert M.rank() == len(vals) and dup == 0

print("")
print("="*76)
print("WHAT IS NOW PROVED, AND WHAT IS NOT")
print("="*76)
print("""
  PROVED, on published results:  every value o_nu(A) is a symplectic character;
  the labels of length r+1 vanish; longer labels fold onto the basis with a sign;
  and the basis is linearly independent.  Equation (Cmu) therefore holds with no
  appeal to a computed range, and with no modification rule outside the two
  classical cases.

  STILL OPEN:  Conjecture (iso) -- the existence of an isolating mu -- which is a
  statement about the Littlewood coefficients c_nu, not about the values o_nu(A).
  Nothing here touches it.
""")
print("DONE")
