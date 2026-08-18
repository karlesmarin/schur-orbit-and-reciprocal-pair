# -*- coding: utf-8 -*-
# t = 4, r = 2: what separates the vanishing concentric shapes from the rest?
#
# Concentric is necessary and not sufficient in the interior.  This battery scores candidate extra
# conditions on the same two sets, and a candidate only counts if it fires on ALL the vanishing ones
# and on NONE of the others.  Printing the two columns side by side is the point: a feature shared
# by both decides nothing, however suggestive it looks on the five.
#
# The candidates are not arbitrary.  They are the shapes the two settled lines suggest:
#   - the t = 2 branch is a self-complementarity, so: is beta symmetric about the common centre,
#     b -> C - b, as a whole set?  (It maps classes to other classes, not each to itself.)
#   - the r = 1 branch is one concentric pair, so: does the interior need TWO?
#   - and the structural readings: equal-size concentric classes, every class symmetric, the
#     quotient components of the concentric pair conjugate or equal, the complement of the pair
#     concentric in its turn.
#
# Authors: Carles Marin, Claude (AI assistant).

t, r = 4, 2
N = t + 2 * r
MAX = 17
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
zs = list(R.gens())
Sym = SymmetricFunctions(QQ)
p, s = Sym.p(), Sym.s()


def psi(lam):
    out = R(0)
    for rho, c in p(s[list(lam)] if lam else s.one()):
        term = R(c)
        for k in rho:
            term *= R((t if k % t == 0 else 0) + sum(zz ** k + zz ** (-k) for zz in zs))
        out += term
    return out


def info(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = [sorted([b for b in beta if b % t == i], reverse=True) for i in range(t)]
    return beta, cls


def conc_pairs(cls):
    out = []
    big = [(i, c) for i, c in enumerate(cls) if len(c) >= 2]
    for x in range(len(big)):
        for y in range(x + 1, len(big)):
            (i, a), (j, b) = big[x], big[y]
            if a[0] + a[-1] == b[0] + b[-1]:
                out.append((i, j, a[0] + a[-1]))
    return out


def quot(cls):
    return [tuple((b - i) // t - (len(cls[i]) - 1 - j) for j, b in enumerate(cls[i]))
            for i in range(t)]


# ---- the candidate conditions -------------------------------------------------------------
def f_beta_symmetric(beta, cls, cp):
    return any(set(C - b for b in beta) == set(beta) for _, _, C in cp)


def f_two_pairs(beta, cls, cp):
    return len(cp) >= 2


def f_equal_sizes(beta, cls, cp):
    return any(len(cls[i]) == len(cls[j]) for i, j, _ in cp)


def f_each_class_symmetric(beta, cls, cp):
    for _, _, C in cp:
        if all(set(C - b for b in cls[k]) == set(cls[k]) for k in range(t) if cls[k]):
            return True
    return False


def f_pair_classes_symmetric(beta, cls, cp):
    return any(set(C - b for b in cls[i]) == set(cls[i]) and
               set(C - b for b in cls[j]) == set(cls[j]) for i, j, C in cp)


def f_quot_equal(beta, cls, cp):
    q = quot(cls)
    return any(q[i] == q[j] for i, j, _ in cp)


def f_rest_concentric(beta, cls, cp):
    for i, j, C in cp:
        rest = [k for k in range(t) if k not in (i, j) and len(cls[k]) >= 2]
        for x in range(len(rest)):
            for y in range(x + 1, len(rest)):
                a, b = cls[rest[x]], cls[rest[y]]
                if a[0] + a[-1] == b[0] + b[-1]:
                    return True
    return False


def f_centre_equals_max(beta, cls, cp):
    return any(C == max(beta) for _, _, C in cp)


def f_all_centres_equal(beta, cls, cp):
    cs = [c[0] + c[-1] for c in cls if len(c) >= 2]
    return len(set(cs)) == 1 and len(cs) >= 2


def f_empty_quot_in_pair(beta, cls, cp):
    """One of the two concentric classes is an arithmetic progression of step t -- equivalently its
    quotient component is empty.  Read off the eight vanishing shapes by hand before being tested
    here, which is why it needs the second column more than any other candidate on the list."""
    q = quot(cls)
    return any(sum(q[i]) == 0 or sum(q[j]) == 0 for i, j, _ in cp)


def f_both_empty_quot(beta, cls, cp):
    q = quot(cls)
    return any(sum(q[i]) == 0 and sum(q[j]) == 0 for i, j, _ in cp)


CAND = [("ONE of the pair has empty quotient", f_empty_quot_in_pair),
        ("BOTH of the pair have empty quotient", f_both_empty_quot),
        ("beta symmetric about the centre", f_beta_symmetric),
        ("two or more concentric pairs", f_two_pairs),
        ("the concentric classes have equal size", f_equal_sizes),
        ("every class symmetric about the centre", f_each_class_symmetric),
        ("both classes of the pair symmetric", f_pair_classes_symmetric),
        ("the pair's quotient components equal", f_quot_equal),
        ("the complementary pair is concentric too", f_rest_concentric),
        ("the centre equals max(beta)", f_centre_equals_max),
        ("all excess classes share one centre", f_all_centres_equal)]

zeros, others = [], []
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        lam = list(l)
        beta, cls = info(lam)
        if any(len(c) == 0 for c in cls):
            continue
        cp = conc_pairs(cls)
        if not cp:
            continue
        rec = (lam, beta, cls, cp)
        (zeros if psi(lam) == 0 else others).append(rec)

print("=" * 78)
print("t = 4, r = 2, |lambda| <= %d : %d concentric shapes, %d of them vanishing"
      % (MAX, len(zeros) + len(others), len(zeros)))
print("=" * 78)
print("")
print("  %-42s %8s %10s" % ("candidate extra condition", "on the 5", "on the rest"))
print("  " + "-" * 64)
for name, f in CAND:
    a = sum(1 for rec in zeros if f(rec[1], rec[2], rec[3]))
    b = sum(1 for rec in others if f(rec[1], rec[2], rec[3]))
    flag = ""
    if a == len(zeros) and b == 0:
        flag = "   <-- SEPARATES"
    elif a == len(zeros):
        flag = "   (necessary, not sufficient)"
    elif b == 0:
        flag = "   (sufficient, not necessary)"
    print("  %-42s %4d/%-4d %6d/%-5d%s" % (name, a, len(zeros), b, len(others), flag))

print("")
print("  the vanishing ones:")
for lam, beta, cls, cp in zeros:
    print("    lam=%-24s beta=%-30s conc=%s"
          % (str(lam), str(beta), str(cp)))

print("")
print("DONE")
