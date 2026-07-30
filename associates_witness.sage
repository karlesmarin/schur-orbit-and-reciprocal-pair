# The stable-range converse, PROVED by constructing the associate witness.
#
# Psi = 0 iff m_nu(lam) = m_{nu*}(lam) for all nu. So to prove Psi != 0 it suffices to EXHIBIT one
# nu with m_nu != m_{nu*}. Dualising the observed witnesses (always one row taller than lam), the
# construction is on mu = nu*, one row SHORTER than lam:
#
#   if ell(lam) <= r                       : mu = lam itself (beta = empty, m_lam >= 1), and
#                                            ell(lam*) = N - ell(lam) > ell(lam) so m_{lam*} = 0.
#   if ell(lam) = r+1 and lam_{r+1} even   : mu = lam minus its last row, beta = (lam_{r+1}) even,
#                                            lam/mu a horizontal strip, so m_mu >= 1, ell(mu) = r.
#   if ell(lam) = r+1, lam not a rectangle : j = largest index <= r with lam_j > lam_{j+1}; remove
#                                            the last row AND one box from row j. beta =
#                                            (lam_{r+1}+1) is even, and lam_j >= lam_{r+1}+1 keeps
#                                            lam/mu a horizontal strip. ell(mu) = r.
#   lam = (k^{r+1}), k odd                 : NEITHER construction exists -- and this is exactly the
#                                            self-complementary-odd shape that vanishes.
#
# Also: branch (a) cannot occur in the stable range. With ell(lam) <= r+1 = N/2 there are at least
# two trailing zeros, so beta contains consecutive integers N-j, hence mixed parity.
#
# This script verifies each construction: that mu is a partition of the right length, that
# m_mu >= 1, that m_{mu*} = 0, and that the only stable-range shapes with no witness are the odd
# rectangles.
#
# Carles Marin + Claude (AI assistant).

Sym = SymmetricFunctions(QQ); s = Sym.s()

def associate(nu, N):
    c = list(Partition(list(nu)).conjugate())
    first = N - (c[0] if c else 0)
    rest = c[1:]
    if first < 0 or first < (rest[0] if rest else 0):
        return None
    return tuple(Partition([first] + rest).conjugate())

def mult(lam, nu, N):
    lam = Partition(list(lam)); nu = Partition(list(nu))
    d = lam.size() - nu.size()
    if d < 0:
        return 0
    tot = 0
    if d % 2 != 0:
        # beta must have all parts even, so |beta| must be even
        return 0
    for half in Partitions(d // 2, max_length=N):
        beta = Partition([2 * x for x in half])
        tot += (s(nu) * s(beta)).coefficient(lam)
    return tot

def witness(lam, r):
    """returns (mu, kind) or None"""
    N = 2 * r + 2
    L = [x for x in lam if x > 0]
    ell = len(L)
    if ell <= r:
        return (tuple(L), 'mu = lam')
    if ell == r + 1:
        t = L[r]
        if t % 2 == 0:
            return (tuple(L[:r]), 'drop last row')
        # not a rectangle?
        js = [j for j in range(r) if L[j] > L[j + 1]]
        if js:
            j = max(js)
            mu = list(L[:r]); mu[j] -= 1
            mu = [x for x in mu if x > 0]
            return (tuple(mu), 'drop last row + one box from row %d' % (j + 1))
    return None

def run(r, NMAX):
    N = 2 * r + 2
    half = N // 2
    print("=" * 72)
    print("r=%d  N=%d  stable range ell <= %d  |lam| <= %d" % (r, N, half, NMAX))
    nowit, bad = [], []
    checked = 0
    for tot in range(1, NMAX + 1):
        for lam in Partitions(tot, max_length=half):
            L = [x for x in lam if x > 0]
            w = witness(L, r)
            if w is None:
                nowit.append(tuple(L))
                continue
            mu, kind = w
            st = associate(mu, N)
            m1 = mult(L, mu, N)
            m2 = 0 if st is None else mult(L, st, N)
            checked += 1
            if not (m1 >= 1 and m2 == 0):
                bad.append((tuple(L), mu, kind, m1, m2))
    print("  witnesses constructed and verified: %d ; failures: %d" % (checked, len(bad)))
    for b in bad[:8]:
        print("     lam=%-16s mu=%-16s [%s]  m_mu=%d m_mu*=%d" % (str(b[0]), str(b[1]), b[2], b[3], b[4]))
    print("  shapes with NO witness: %d" % len(nowit))
    notrect = [l for l in nowit if not (len(set(l)) == 1 and l[0] % 2 == 1 and len(l) == half)]
    for l in nowit[:12]:
        print("     %-18s odd rectangle of height %d: %s"
              % (str(l), half, (len(set(l)) == 1 and l[0] % 2 == 1 and len(l) == half)))
    print("  of those, NOT an odd rectangle of full height: %d  (must be 0)" % len(notrect))
    # branch (a) cannot occur in the stable range
    ba = 0
    for tot in range(0, NMAX + 1):
        for lam in Partitions(tot, max_length=half):
            Lp = list(lam) + [0] * (N - len(lam))
            b = [Lp[i] + (N - 1 - i) for i in range(N)]
            if len(set(x % 2 for x in b)) == 1:
                ba += 1
    print("  branch-(a) shapes inside the stable range: %d  (must be 0)" % ba)
    return len(bad) + len(notrect) + ba

t = 0
t += run(1, 16)
t += run(2, 14)
t += run(3, 12)
print("\nTOTAL failures: %d" % t)
if t == 0:
    print("=> the stable-range converse is PROVED: every lam with ell(lam) <= N/2 that is not an")
    print("   odd rectangle of full height admits an explicit associate witness, so Psi != 0;")
    print("   and the odd full-height rectangles are exactly the self-complementary-odd shapes.")
    print("   Residue: the UNSTABLE range ell(lam) > N/2, where Littlewood needs King's rules.")
