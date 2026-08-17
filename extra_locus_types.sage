# Problem 10.6, the computable half: does Theorem 5.2 (the extra locus) have an analogue for the
# universal symplectic and orthogonal characters?
#
# THE QUESTION.  [AK25, Thm 5.3] says f_lambda(x, Y, zeta Y, ..., zeta^{t-1} Y) = f_lambda(x) iff
# lambda is a t-core, for f of type s, sp or o.  Theorem 5.2 of this paper says that at OUR
# specialization -- the free part being one reciprocal pair, and equality allowed up to SIGN -- the
# criterion acquires exactly one extra family, and only for t even.  Problem 10.6 asks whether sp
# and o behave the same way, and whether the extra family is again indexed by an order-two element.
#
# HOW IT IS TESTED.  Everything happens in Lambda, where adjoining letters is a specialization of
# power sums:
#
#     the full orbit mu_t   ->   p_k  gets  + t   when t | k,  + 0 otherwise
#     one reciprocal pair   ->   p_k  gets  + z^k + z^-k
#
# so the two sides of the criterion are
#
#     LHS = f_lambda  at  p_k = t*[t|k] + z^k + z^-k          (the orbit and the pair)
#     RHS = f_lambda  at  p_k = z^k + z^-k                    (the pair alone)
#
# and the locus is { lambda : LHS = +- RHS }.  For type s this is exactly the test behind Figure 4,
# which is the control: it must return the t-cores plus, at even t, the family of Theorem 5.2 and
# nothing else.
#
# TWO EVALUATION POINTS.  An equality of Laurent polynomials is tested at two unrelated rational z,
# and a lambda is admitted only if both agree; one point alone would let a coincidence through.
#
# Authors: Carles Marin + Claude (AI assistant).

Sym = SymmetricFunctions(QQ)
p = Sym.p()
BASES = [("s", Sym.s()), ("sp", Sym.sp()), ("o", Sym.o())]

# THE RANGE IS PART OF THE QUESTION, and getting it wrong makes the test measure something else.
# The right-hand side f_lambda(free part alone) VANISHES once ell(lambda) exceeds the rank of the
# free part: s_lambda(z,zbar) = 0 for ell(lambda) >= 3, and sp_lambda(z) = 0 for ell(lambda) >= 2.
# Beyond that range "LHS = +- RHS" is not the independence criterion at all -- it degenerates into
# "LHS = 0", which is the zero locus and a different question. A first run of this script with three
# rows and a rank-one free part reported hundreds of "extras" that were nothing but vanishing shapes.
# So each type is run with a free part of rank at least ell(lambda), and never past it.
#
#     s   : free part (z, 1/z) is rank 2 for Schur      -> ell(lambda) <= 2
#     sp  : free part (z_1, ..., z_R) is rank R         -> ell(lambda) <= R
#     o   : same
#
ZS = [QQ(3) / 2, QQ(5) / 3]          # the two evaluation points
TS = [3, 4, 5, 6]
MAXSIZE = 14
MAXROWS = 2                          # the range where the criterion is a criterion
#
# A CONFOUND, and the control for it.  A Schur alphabet of p reciprocal pairs is 2p variables and
# admits ell(lambda) <= 2p; a symplectic or orthogonal one of p pairs has rank p and admits only
# ell(lambda) <= p.  So comparing type s at p=1 with sp and o at p=2 compares two-letter free parts
# against four-letter ones, and any difference could be the SIZE rather than the TYPE.  Type s is
# therefore run at p=2 as well, on the same two-row range: if its extra family survives there, the
# difference is the type; if it does not, the extra family is a fact about the minimal alphabet.
CONFIGS = [("s", 1), ("s", 2), ("sp", 2), ("o", 2)]
FREE = None


def specialize(f, t, z, with_orbit, pairs):
    """f in Lambda at the alphabet (mu_t if with_orbit) u {z_i, 1/z_i : i <= pairs}, the z_i
    being z, z^2, ... so that they are multiplicatively independent enough for the test."""
    out = QQ(0)
    for rho, c in p(f):
        term = QQ(1)
        for k in rho:
            v = QQ(0)
            for i in range(1, pairs + 1):
                w = z**i
                v += w**k + w**(-k)
            if with_orbit and k % t == 0:
                v += t
            term *= v
        out += c * term
    return out


def tcore(lam, t):
    """the t-core, by repeatedly removing rim hooks via the beta set"""
    n = max(len(lam), 1) + t
    beta = [lam[i] + n - 1 - i if i < len(lam) else n - 1 - i for i in range(n)]
    changed = True
    while changed:
        changed = False
        s = set(beta)
        for b in sorted(beta, reverse=True):
            if b - t >= 0 and (b - t) not in s:
                beta.remove(b); beta.append(b - t); changed = True
                break
    beta = sorted(beta, reverse=True)
    out = [beta[i] - (n - 1 - i) for i in range(n)]
    return tuple(v for v in out if v > 0)


print("=" * 92)
print("PROBLEM 10.6 -- the extra locus for the three types")
print("=" * 92)

RESULT = {}
BASIS = dict(BASES)
for tag, pairs in CONFIGS:
    B = BASIS[tag]
    for t in TS:
        core_hits, extra_hits = [], []
        for size in range(0, MAXSIZE + 1):
            for lam in Partitions(size, max_length=MAXROWS):
                lam = tuple(lam)
                f = B[list(lam)] if lam else B.one()
                ok = True
                for z in ZS:
                    L = specialize(f, t, z, True, pairs)
                    R = specialize(f, t, z, False, pairs)
                    if R == 0:          # out of range: the criterion has degenerated
                        ok = False
                        break
                    if L != R and L != -R:
                        ok = False
                        break
                if not ok:
                    continue
                (core_hits if tcore(lam, t) == lam else extra_hits).append(lam)
        RESULT[(tag, pairs, t)] = (core_hits, extra_hits)
        print("  %-3s %d pair(s), %d free letters, t=%d :  locus %3d   t-cores %3d   EXTRA %3d   %s"
              % (tag, pairs, 2 * pairs, t, len(core_hits) + len(extra_hits), len(core_hits),
                 len(extra_hits), ("  extras: " + str(extra_hits[:6])) if extra_hits else ""))

print("")
print("=" * 92)
print("CONTROL -- type s must reproduce Theorem 5.2: extras only at even t, and of the stated form")
print("=" * 92)
for t in TS:
    core_hits, extra_hits = RESULT[("s", 1, t)]
    two_row = [l for l in extra_hits if len(l) == 2]
    predicted = []
    if t % 2 == 0:
        for l2 in range(t // 2, t):
            l1 = l2 + 3 * t // 2 - 1
            if l1 + l2 <= MAXSIZE:
                predicted.append((l1, l2))
    match = sorted(two_row) == sorted(predicted)
    print("  t=%d : extras %-28s predicted %-28s  %s"
          % (t, str(sorted(extra_hits)), str(sorted(predicted)),
             "MATCH" if match and len(two_row) == len(extra_hits) else "*** MISMATCH ***"))

print("")
print("=" * 92)
print("THE ANSWER, per type: does the criterion acquire extras, and only at even t?")
print("=" * 92)
for tag, pairs in CONFIGS:
    odd = sum(len(RESULT[(tag, pairs, t)][1]) for t in TS if t % 2)
    even = sum(len(RESULT[(tag, pairs, t)][1]) for t in TS if t % 2 == 0)
    print("  %-3s with %d pair(s) : extras at odd t = %d ; at even t = %d  ->  %s"
          % (tag, pairs, odd, even,
             "an extra family, only at even t" if odd == 0 and even > 0
             else ("NO extra family at all" if odd == 0 and even == 0
                   else "extras at odd t too")))

print("")
print("Read the two type-s rows against each other. If s keeps its extra family at two pairs, the")
print("difference from sp and o is the TYPE. If it loses it, the extra family belongs to the")
print("minimal alphabet, and Problem 10.6 is asking about minimality rather than about type.")
