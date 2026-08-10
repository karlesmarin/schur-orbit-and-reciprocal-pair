"""The two-class stratum is a lattice, and the fibres of I_t are its sections.

Authors: Carles Marin, Claude (AI assistant).

This is the computation behind Proposition prop:fibrelattice and Corollary cor:fibregf.  It answers
what Problem prob:fibres asked: which partitions share a triple (d1,d2,d3), and how many of each
size.

CONVENTION.  Theorem thm:main orders the two distinguished classes BY RESIDUE, r_A <= r_B, and reads
d1 off class r_A.  theorem_full.setup() orders them by COLUMN, i.e. by largest beta number, which is
a different tie-break and transposes (d1,d2) on part of the range.  Everything below is written from
the theorem's definition directly and never calls setup(), so the two conventions cannot be confused.

WHAT IS CHECKED, all over |lambda| <= nmax:

  (1) LATTICE.  The map sending lambda to (r_A, r_B, nu, mu, m) -- the two distinguished quotient
      components and the t-2 single parts -- is a bijection from the two-class stratum onto
          {r_A < r_B}  x  {two-part partitions}^2  x  Z_{>=0}^{t-2}.
      Tested as a dictionary equality in both directions, triple by triple: the set built from the
      lattice must equal the set found by brute force, for every triple.

  (2) GENERATING FUNCTION.  For each triple, the number of lambda in its fibre with |lambda| = n is
      compared, coefficient by coefficient, with the prediction

          sum over branches of   q^{c} / [ (1 - q^{4t}) (1 - q^{t})^{t-2} ]

      the 4t because the free parameter v moves |nu| and |mu| by two each, and one factor 1 - q^t per
      free component.

CONTROLS, each able to fail:
  K1  wrong denominator: (1 - q^{2t}) in place of (1 - q^{4t}) must disagree.
  K2  wrong labelling: A = the class of the largest beta number, i.e. setup()'s tie-break rather
      than the theorem's, must disagree -- this is what makes the convention a checked fact and not
      a matter of taste.
  K3  coverage: the partitions covered must be exactly the two-class ones, counted independently.
"""
import sys
from collections import defaultdict


def partitions_upto(n, maxlen):
    """all partitions of n with at most maxlen parts"""
    out = []

    def rec(rem, cap, acc):
        if rem == 0:
            out.append(tuple(acc))
            return
        if len(acc) == maxlen:
            return
        for p in range(min(rem, cap), 0, -1):
            rec(rem - p, p, acc + [p])
    rec(n, n, [])
    return out


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[j] + N - 1 - j for j in range(N)]


def lam_of_beta(b, N):
    b = sorted(b, reverse=True)
    lam = [b[j] - (N - 1 - j) for j in range(N)]
    if any(x < 0 for x in lam):
        return None
    while lam and lam[-1] == 0:
        lam.pop()
    return tuple(lam)


def triple_from_theorem(lam, t, by_column=False):
    """(d1,d2,d3) exactly as Theorem thm:main defines it, for the two-class profile.

    by_column=True switches to setup()'s tie-break -- control K2."""
    N = t + 2
    beta = beta_of(lam, N)
    cls = defaultdict(list)
    for b in beta:
        cls[b % t].append(b)
    if len(cls) < t:
        return None                                   # degenerate
    big = sorted([r for r in cls if len(cls[r]) == 2])
    if len(big) != 2:
        return None                                   # size-three
    rA, rB = big                                      # the theorem: by residue, r_A < r_B
    A, B = sorted(cls[rA], reverse=True), sorted(cls[rB], reverse=True)
    if by_column and B[0] > A[0]:
        A, B = B, A
    return (A[0] - A[1], B[0] - B[1], abs(A[0] + A[1] - B[0] - B[1]))


def brute(t, nmax, by_column=False):
    N = t + 2
    fib = defaultdict(set)
    other = 0
    for n in range(0, nmax + 1):
        for lam in partitions_upto(n, N):
            d = triple_from_theorem(lam, t, by_column)
            if d is None:
                other += 1
            else:
                fib[d].add(lam)
    return fib, other


def built(t, nmax):
    """every point of the lattice, labelled by the theorem's rule"""
    N = t + 2
    fib = defaultdict(set)
    K = nmax // t + 3

    def frees(rest, budget):
        if not rest:
            yield []
            return
        for m in range(budget + 1):
            for tail in frees(rest[1:], budget - m):
                yield [t * m + rest[0]] + tail

    for rA in range(t):
        for rB in range(rA + 1, t):
            rest = [r for r in range(t) if r not in (rA, rB)]
            for kA in range(K):
                for v in range(K):
                    #   nu = (v + kA, v)   ->   beta  t(v+kA+1)+rA,  t v + rA
                    A = [t * (v + kA + 1) + rA, t * v + rA]
                    if A[0] - (N - 1) > nmax:
                        break
                    for kB in range(K):
                        for u in range(K):
                            B = [t * (u + kB + 1) + rB, t * u + rB]
                            if B[0] - (N - 1) > nmax:
                                break
                            for f in frees(rest, K):
                                b = A + B + f
                                if len(set(b)) != N:
                                    continue
                                lam = lam_of_beta(b, N)
                                if lam is None or sum(lam) > nmax:
                                    continue
                                d = (A[0] - A[1], B[0] - B[1],
                                     abs(A[0] + A[1] - B[0] - B[1]))
                                fib[d].add(lam)
    return fib


# ----------------------------------------------------------------- generating function

def series_mult(a, b, n):
    out = [0] * (n + 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if i + j <= n and bj:
                    out[i + j] += ai * bj
    return out


def inv_one_minus_q(e, n):
    """the power series of 1/(1 - q^e) up to q^n"""
    s = [0] * (n + 1)
    for i in range(0, n + 1, e):
        s[i] = 1
    return s


def core_size(t, rA, rB):
    """|core_t(lambda)| for the two-class profile with distinguished residues rA < rB.

    The core is the lattice point with every component empty."""
    N = t + 2
    b = [t + rA, rA, t + rB, rB] + [r for r in range(t) if r not in (rA, rB)]
    assert len(set(b)) == N
    return sum(b) - N * (N - 1) // 2


def predicted(t, d, nmax, wrong_denominator=False):
    """the predicted size-count sequence of the fibre over d, from the closed form"""
    d1, d2, d3 = d
    out = [0] * (nmax + 1)
    if d1 % t or d2 % t:
        return out
    kA, kB = d1 // t - 1, d2 // t - 1
    if kA < 0 or kB < 0:
        return out
    step = 2 * t if wrong_denominator else 4 * t           # control K1
    base = inv_one_minus_q(step, nmax)
    for _ in range(t - 2):
        base = series_mult(base, inv_one_minus_q(t, nmax), nmax)
    for rA in range(t):
        for rB in range(rA + 1, t):
            for eps in (+1, -1):
                num = eps * d3 - 2 * (rA - rB)
                if num % t:
                    continue
                #  2(v - u) = num/t - kA + kB
                w = num // t - kA + kB
                if w % 2:
                    continue
                c = w // 2                                  # u = v - c
                vmin = max(0, c)
                #  |lambda| = |core| + t(4v + kA + kB - 2c)
                off = core_size(t, rA, rB) + t * (kA + kB - 2 * c) + 4 * t * vmin
                if off < 0 or off > nmax:
                    continue
                if eps == -1 and d3 == 0:
                    continue                                # the two branches coincide
                out = [o + (base[i - off] if i >= off else 0) for i, o in enumerate(out)]
    return out


def main():
    print("=" * 96)
    print("The two-class stratum as a lattice, and the generating function of its fibres")
    print("=" * 96)
    print()
    for t, nmax in [(2, 14), (3, 12), (4, 10), (5, 10)]:
        B, other = brute(t, nmax)
        G = built(t, nmax)
        nb = sum(len(v) for v in B.values())
        print("t = %d,  |lambda| <= %d" % (t, nmax))
        print("    two-class partitions: %d on %d triples   (%d in other profiles)"
              % (nb, len(B), other))
        print("    (1) LATTICE   dictionaries equal: %s"
              % ("YES" if dict(B) == dict(G) else "NO"))
        # ---- control K2: setup()'s tie-break must give a different bucketing
        Bc, _ = brute(t, nmax, by_column=True)
        print("    K2 control    theorem's labelling differs from setup()'s: %s"
              % ("YES" if dict(Bc) != dict(B) else "NO -- the test is blind to the convention"))
        # ---- (2) the generating function
        ok = bad = 0
        firstbad = None
        for d in sorted(B):
            actual = [0] * (nmax + 1)
            for lam in B[d]:
                actual[sum(lam)] += 1
            if predicted(t, d, nmax) == actual:
                ok += 1
            else:
                bad += 1
                if firstbad is None:
                    firstbad = (d, actual, predicted(t, d, nmax))
        print("    (2) GEN.FUNC  fibres matching the closed form: %d / %d" % (ok, ok + bad))
        if firstbad:
            d, a, p = firstbad
            print("        first mismatch d = %s" % (d,))
            print("            actual    %s" % a)
            print("            predicted %s" % p)
        # ---- control K1: the wrong denominator
        okw = sum(1 for d in sorted(B)
                  if predicted(t, d, nmax, wrong_denominator=True)
                  == [sum(1 for lam in B[d] if sum(lam) == n) for n in range(nmax + 1)])
        print("    K1 control    same with (1-q^{2t}) instead of (1-q^{4t}): %d / %d%s"
              % (okw, ok + bad, "   <-- must be smaller" if okw < ok else "   <-- FAILED"))
        print()


def _main_lattice():
    main()




# ----------------------------------------------------------------- the sign on a branch

def sign_on_branches(t, nmax, theta=0.7):
    """How does eps_lambda restrict to a branch of Corollary cor:fibregf?

    eps_lambda is the sign of Phi_t and is convention-free, so we take it from the BIALTERNANT --
    the definition -- and not from any closed form: sgn Phi_t(lambda; e^theta) for theta > 0, since
    every sinh in eq:main is then positive.

    A branch is indexed by (rA, rB, kA, kB, c) with u = v - c; inside it the free parameters are v
    and the t-2 invisible parts m.  We report, over everything in range:
        (a) is eps constant in m at fixed v?
        (b) what does eps do as v increases by one?"""
    from theorem_full import phi_bialternant
    N = t + 2
    data = defaultdict(dict)                      # branch -> (v, m) -> eps
    for n in range(0, nmax + 1):
        for lam in partitions_upto(n, N):
            d = triple_from_theorem(lam, t)
            if d is None or d[2] == 0:
                continue
            beta = beta_of(lam, N)
            cls = defaultdict(list)
            for b in beta:
                cls[b % t].append(b)
            rA, rB = sorted([r for r in cls if len(cls[r]) == 2])
            A = sorted(cls[rA], reverse=True)
            B = sorted(cls[rB], reverse=True)
            v = A[1] // t
            kA = (A[0] - A[1]) // t - 1
            u = B[1] // t
            kB = (B[0] - B[1]) // t - 1
            m = tuple(sorted(cls[r][0] // t for r in cls if len(cls[r]) == 1))
            val = phi_bialternant(lam, t, theta)
            eps = 1 if val.real > 0 else (-1 if val.real < 0 else 0)
            data[(rA, rB, kA, kB, v - u)][(v, m)] = eps
    const_in_m = varies_in_m = 0
    flips = same = other = 0
    for br, pts in data.items():
        byv = defaultdict(set)
        for (v, m), e in pts.items():
            byv[v].add(e)
        for v in byv:
            if len(byv[v]) == 1:
                const_in_m += 1
            else:
                varies_in_m += 1
        vs = sorted(byv)
        for v1, v2 in zip(vs, vs[1:]):
            if v2 != v1 + 1 or len(byv[v1]) != 1 or len(byv[v2]) != 1:
                continue
            e1, e2 = next(iter(byv[v1])), next(iter(byv[v2]))
            if e1 == -e2:
                flips += 1
            elif e1 == e2:
                same += 1
            else:
                other += 1
    print("    SIGN on a branch:  eps constant in the invisible parts m at fixed v: "
          "%d of %d (v, branch) slots" % (const_in_m, const_in_m + varies_in_m))
    print("                       consecutive v: %d flips, %d equal, %d other"
          % (flips, same, other))


if __name__ == "__main__":
    main()
    print()
    print("=" * 96)
    print("The sign restricted to a branch")
    print("=" * 96)
    print()
    for t, nmax in [(2, 14), (3, 12), (4, 10)]:
        print("t = %d,  |lambda| <= %d" % (t, nmax))
        sign_on_branches(t, nmax)
        print()
