# THE VANISHING CRITERION, DERIVED -- self-complementarity of lambda, for every r.
#
# Derivation (not a fit).  N = 2r+2 letters (1, -1, z_1,1/z_1, ..., z_r,1/z_r), beta = lam + rho.
# Laplace the bialternant on the two mu_2 rows -- the all-ones row 1^{beta_j} and the parity row
# (-1)^{beta_j}.  Their 2x2 minor on columns {j1,j2} is (-1)^{beta_j2} - (-1)^{beta_j1}: zero on
# equal parity, +-2 on opposite. So
#
#     Psi  =  2 * sum_{e in E, o in O}  sigma(e,o) * a_{beta \ {e,o}}(Z, Zbar)
#
# where a_S is the 2r x 2r alternant with rows z_i^{+-b}, b in S. The alphabet (Z,Zbar) has
# PRODUCT 1, hence
#     a_{S+c} = (prod of the letters)^c * a_S = a_S          (the det-twist is trivial)
#     a_{-S}  = (-1)^r a_S                                   (r row swaps)
# and therefore  a_S = (-1)^r a_{c-S}  for every c.  The terms pair up under the reflection
# S -> c-S. A pairing of the whole sum with opposite Laplace signs kills it, and
# beta \ {e,o} -> c - (beta \ {e,o}) lands back in the term set exactly when beta = c - beta.
#
# beta_i + beta_{N+1-i} = c  <=>  lam_i + lam_{N+1-i} = c - N + 1 =: w  constant
#                          <=>  LAMBDA IS SELF-COMPLEMENTARY in a w x N rectangle.
# N is even, so c even <=> w odd.
#
#   CLAIM.  Psi = 0  <=>  (a) all beta_j have the same parity   [the -1 row duplicates the +1 row]
#                    OR   (b) lam is self-complementary in a w x N rectangle with w ODD.
#
# At r=1 this is equivalent to the published branch (b): with |E|=2 and beta sorted, the only
# parity pairing compatible with sum(E)=sum(O) is {b1,b4} vs {b2,b3}, which IS the symmetry.
# The ratio law of ratio_law.sage is the moment consequence: if c is even the reflection preserves
# parity classes, so E and O are each symmetric about c/2 and (r+1)sum(E) = k*T follows.
# The stable-range theorem (odd rectangles (k^{r+1}), §14) is the sub-case lam_N = 0, w = lam_1.
#
# Verified below TWO-SIDED against exact Laurent determinants, r = 1,2,3.
# Carles Marin + Claude (AI assistant).

import random
random.seed(int(20260730))

def selfcomp_odd(lam, N):
    """lam padded to N parts: is lam_i + lam_{N+1-i} constant and odd?"""
    L = list(lam) + [0] * (N - len(lam))
    if len(L) != N:
        return False
    w = L[0] + L[N - 1]
    if w % 2 == 0:
        return False
    return all(L[i] + L[N - 1 - i] == w for i in range(N))

def run(r_pairs, NMAX):
    N = 2 * r_pairs + 2
    rho = list(range(N - 1, -1, -1))

    def hgen(alpha, D, ring):
        PK = PowerSeriesRing(ring, 's', default_prec=D + 2); s = PK.gen()
        G = prod(1 / (1 - ring(a) * s) for a in alpha)
        return [G[k] for k in range(D + 1)]

    def jt(lam, alpha, ring):
        L = list(lam) + [0] * (N - len(lam))
        D = L[0] + N + 1
        h = hgen(alpha, D, ring)
        H = lambda k: ring(0) if k < 0 else h[k]
        return matrix(ring, N, N, lambda i, j: H(L[i] - i + j)).det()

    A_ONE = [QQ(1), QQ(-1)] + [QQ(1)] * (2 * r_pairs)
    PTS = []
    for _ in range(4):
        A = [QQ(1), QQ(-1)]
        for _ in range(r_pairs):
            z = QQ(int(random.randint(3, 120))) / QQ(int(random.randint(3, 120)))
            A += [z, 1 / z]
        PTS.append(A)
    Rz = LaurentPolynomialRing(QQ, ['z%d' % i for i in range(1, r_pairs + 1)])
    A_EXACT = [Rz(1), Rz(-1)]
    for z in Rz.gens():
        A_EXACT += [z, z**-1]

    def branch_a(lam):
        L = list(lam) + [0] * (N - len(lam))
        b = [L[i] + rho[i] for i in range(N)]
        p = set(x % 2 for x in b)
        return len(p) == 1

    def predict(lam):
        return branch_a(lam) or selfcomp_odd(lam, N)

    print("=" * 74)
    print("r = %d,  N = %d,  |lam| <= %d" % (r_pairs, N, NMAX))
    fn, fp, nvan, ntot = [], [], 0, 0
    for S in range(0, NMAX + 1):
        for lam in Partitions(S, max_length=N):
            ntot += 1
            p = predict(lam)
            # exact status, with the cheap necessary sieve first
            if jt(lam, A_ONE, QQ) != 0 or any(jt(lam, A, QQ) != 0 for A in PTS):
                z = False
            else:
                z = (jt(lam, A_EXACT, Rz) == 0)
            if z:
                nvan += 1
            if z and not p:
                fn.append(tuple(lam))
            if p and not z:
                fp.append(tuple(lam))
    print("  scanned %d, exact vanishers %d" % (ntot, nvan))
    print("  FALSE NEGATIVES (vanishes, criterion silent) : %d" % len(fn))
    for l in fn[:10]:
        print("     ", l)
    print("  FALSE POSITIVES (criterion says zero, is not): %d" % len(fp))
    for l in fp[:10]:
        print("     ", l)
    if not fn and not fp:
        print("  ==> CRITERION EXACT on this range, both directions.")
    # how the two branches share the work, and that (b) is not vacuous
    nb = len([1 for S in range(NMAX + 1) for lam in Partitions(S, max_length=N)
              if selfcomp_odd(lam, N) and not branch_a(lam)])
    na = len([1 for S in range(NMAX + 1) for lam in Partitions(S, max_length=N) if branch_a(lam)])
    print("  branch (a) only: %d    branch (b) self-complementary-odd, not (a): %d" % (na, nb))

    # is the parity of w doing real work?  count self-complementary with w EVEN that do NOT vanish
    we = [tuple(lam) for S in range(NMAX + 1) for lam in Partitions(S, max_length=N)
          if (lambda L: len(L) == N and (L[0] + L[N - 1]) % 2 == 0
              and all(L[i] + L[N - 1 - i] == L[0] + L[N - 1] for i in range(N)))
             (list(lam) + [0] * (N - len(lam)))]
    print("  CONTROL: self-complementary with w EVEN: %d in range (criterion predicts NONZERO)" % len(we))
    bad = [l for l in we if not branch_a(l) and jt(l, A_ONE, QQ) == 0
           and not any(jt(l, A, QQ) != 0 for A in PTS) and jt(l, A_EXACT, Rz) == 0]
    print("           of those, actually zero: %d   (should be 0 -- w parity is load-bearing)" % len(bad))
    for l in bad[:6]:
        print("     ", l)

run(1, 28)
run(2, 26)
run(3, 22)
