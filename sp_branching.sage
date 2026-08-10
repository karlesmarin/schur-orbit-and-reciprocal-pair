# Authors: Carles Marin + Claude (AI assistant).
# Problem 10.1, symplectic column: WHY there is no law of the orthogonal kind, and what replaces it.
#
# Section 8 rests on  o_nu(W,1,-1) = sp_nu(W): adjoining the two letters collapses the orthogonal
# universal character to a SINGLE symplectic one.  The symplectic column has no such law, and the
# reason is the general branching rule for universal symplectic functions,
#
#     sp_nu(X u Y) = sum_{gamma subset nu} sp_gamma(X) * sp_{nu/gamma}(Y),
#
# which at Y = (1,-1) makes the coefficient of sp_gamma(W) in sp_nu(W,1,-1) the value of a SKEW
# universal symplectic character at (1,-1).  For the orthogonal character those values happen to be
# concentrated on gamma = nu; for the symplectic one they are not.
#
# The rule predicts something testable without computing any skew character: the coefficient must
# depend on nu and gamma ONLY through the skew shape nu/gamma.  That is what this checks.

Sym = SymmetricFunctions(QQ)
p, O, SP = Sym.p(), Sym.o(), Sym.sp()

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

def skewkey(nu, gamma):
    """the skew shape nu/gamma, normalised: the tuple of (row, first col, last col) cells."""
    g = list(gamma) + [0]*len(nu)
    cells = []
    for i, n in enumerate(nu):
        if n > g[i]:
            cells.append((i, g[i], n))
    # translate rows so the shape is read up to vertical shift
    if not cells:
        return ()
    r0 = cells[0][0]
    return tuple((r - r0, a, b) for (r, a, b) in cells)

print("="*78)
print("CONTROL -- the orthogonal identity of Section 8")
print("="*78)
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
print("="*78)
print("PREDICTION -- the coefficient depends only on the skew shape nu/gamma")
print("="*78)
seen = {}
clash = 0
tested = 0
for size in range(0, 9):
    for nu in Partitions(size):
        nu = tuple(nu)
        g = SP[list(nu)] if nu else SP.one()
        exp = SP(both(g))
        coeffs = {}
        for gamma, c in exp:
            coeffs[tuple(gamma)] = c
        for gamma in Partitions(size):
            gamma = tuple(gamma)
            if len(gamma) > len(nu):
                continue
            gg = list(gamma) + [0]*len(nu)
            if any(gg[i] > nu[i] for i in range(len(nu))):
                continue
        for gamma, c in coeffs.items():
            k = skewkey(nu, gamma)
            tested += 1
            if k in seen and seen[k][0] != c:
                clash += 1
                if clash <= 5:
                    print("  CLASH %s: nu=%s gives %s, nu=%s gave %s"
                          % (str(k), str(nu), c, str(seen[k][1]), seen[k][0]))
            else:
                seen.setdefault(k, (c, nu))
print("  %d coefficients over %d distinct skew shapes: %d clashes   <-- 0 confirms the rule"
      % (tested, len(seen), clash))

print("")
print("="*78)
print("THE VALUES -- sp_{nu/gamma}(1,-1) read off, by skew shape size")
print("="*78)
bysize = {}
for k, (c, nu) in seen.items():
    n = sum(b - a for (_, a, b) in k)
    bysize.setdefault(n, set()).add(c)
for n in sorted(bysize):
    print("  |nu/gamma| = %-2d : values %s" % (n, sorted(bysize[n])))
