# -*- coding: utf-8 -*-
# *** SU CONTEO ES FALSO. NO CITAR EL 30/30 NI EL 0/635. ***
#
# closed_under_reflection() reporta que el conjunto de terminos es cerrado bajo T -> C-T con todas
# las parejas de signo opuesto en las 30 formas nulas y en ninguna de las 635 restantes, y lo lei
# como un criterio.  Es falso.  parity2.sage calcula el mismo signo por dos rutas independientes,
# pasa un test de aceptacion sobre el caso impreso en quadrant_debug_OUT.txt, y encuentra que la
# cancelacion POR PAREJAS ocurre en 2 formas, no en 30.  Las otras 28 se anulan por una
# combinacion de los terminos que no es un emparejamiento.
#
# El error de metodo: di por bueno un conteo agregado sin contrastarlo con un solo caso, y luego
# el caso que mire a mano resulto ser una de las dos que si encajaban.  Se conserva el fichero
# porque el fallo es la parte util.
#
# The 65 zeros that are not self-complementary: what cancels?
#
# Structure first, features second.  In the bialternant the t columns carrying the roots of unity
# have entries zeta^{k beta_i}, which depend on beta_i ONLY THROUGH ITS RESIDUE mod t.  Two rows in
# the same class are therefore equal on those columns, so a Laplace expansion along them keeps only
# the selections that take one row per residue class -- a transversal -- and the t x t minor is the
# same Vandermonde constant for every one of them.  Hence
#
#     Psi  ~  sum_T  sgn(T) A(T),
#
# T running over the 2r-subsets whose complement is a transversal, and A(T) the 2r x 2r alternant
# on the columns z_1, 1/z_1, ..., z_r, 1/z_r.  Vanishing is a cancellation among those terms, and
# the natural candidate for what pairs them is the reflection about the concentric centre,
# T -> C - T, since concentric is exactly what the sweep found to be necessary.
#
# Tested here:
#   CONTROL   the decomposition must reproduce Psi up to the constant: sum_T sgn A(T) = 0 exactly
#             when Psi = 0, on both the zeros and a sample of non-zeros.
#   R1        is the term set closed under T -> C - T?
#   R2        and when it is, do the paired terms carry OPPOSITE signs -- which is what makes the
#             sum collapse -- or the same, which would make it double instead?
#
# Authors: Carles Marin, Claude (AI assistant).

t, r = 4, 2
N = t + 2 * r
MAX = 24
R = LaurentPolynomialRing(QQ, ['z0', 'z1'])
zs = list(R.gens())
COLS = []
for zz in zs:
    COLS.append(('p', zz))
    COLS.append(('m', zz))


def alt(T):
    """The 2r x 2r alternant on rows T, columns z_1, 1/z_1, ..., z_r, 1/z_r."""
    M = matrix(R, len(T), len(T),
               lambda a, b: (zs[b // 2] ** T[a]) if b % 2 == 0 else (zs[b // 2] ** (-T[a])))
    return M.det()


def info(lam):
    lam = list(lam) + [0] * (N - len(lam))
    beta = [lam[i] + N - 1 - i for i in range(N)]
    cls = {}
    for k, b in enumerate(beta):
        cls.setdefault(b % t, []).append(k)
    return beta, cls


def decompose(beta, cls):
    """[(sorted index tuple of the 2r kept rows, sign)] over all transversals."""
    out = []
    keys = sorted(cls)
    if len(keys) < t:
        return out
    import itertools
    for pick in itertools.product(*[cls[k] for k in keys]):
        rest = tuple(sorted(set(range(N)) - set(pick)))
        perm = list(pick) + list(rest)
        sgn = 1
        for a in range(N):
            for b in range(a + 1, N):
                if perm[a] > perm[b]:
                    sgn = -sgn
        out.append((rest, sgn))
    return out


def conc_centres(beta, cls):
    out = []
    big = [(k, sorted([beta[i] for i in v], reverse=True)) for k, v in cls.items() if len(v) >= 2]
    for x in range(len(big)):
        for y in range(x + 1, len(big)):
            if big[x][1][0] + big[x][1][-1] == big[y][1][0] + big[y][1][-1]:
                out.append(big[x][1][0] + big[x][1][-1])
    return sorted(set(out))


print("=" * 78)
print("The Laplace decomposition of the vanishing, t = 4, r = 2, |lambda| <= %d" % MAX)
print("=" * 78)

rows = []
for size in range(MAX + 1):
    for l in Partitions(size, max_length=N):
        lam = list(l)
        beta, cls = info(lam)
        if len(cls) < t:
            continue
        C = conc_centres(beta, cls)
        if not C:
            continue
        dec = decompose(beta, cls)
        total = sum(sg * alt([beta[i] for i in rest]) for rest, sg in dec)
        rows.append((lam, beta, cls, C, dec, total == 0))

zeros = [x for x in rows if x[5]]
nonz = [x for x in rows if not x[5]]
print("")
print("  concentric shapes: %d, of which the signed sum vanishes: %d" % (len(rows), len(zeros)))
print("  CONTROL: the decomposition is the value up to a constant, so these counts must match")
print("           the earlier sweep at the same range.")


def selfcomp(beta, C):
    return any(set(c - b for b in beta) == set(beta) for c in C)


def closed_under_reflection(beta, dec, C):
    """Is the multiset of kept-row VALUE sets closed under T -> c - T, and with which signs?"""
    for c in C:
        table = {}
        for rest, sg in dec:
            T = tuple(sorted(beta[i] for i in rest))
            table[T] = table.get(T, 0) + sg
        ok, opposite, same = True, 0, 0
        for T, w in table.items():
            Tr = tuple(sorted(c - x for x in T))
            if Tr not in table:
                ok = False
                break
            if w == -table[Tr]:
                opposite += 1
            elif w == table[Tr]:
                same += 1
        if ok:
            return (True, opposite, same)
    return (False, 0, 0)


for label, group in (("VANISHING", zeros), ("NOT vanishing", nonz)):
    n_sc = n_closed = n_opp = 0
    for lam, beta, cls, C, dec, z in group:
        if selfcomp(beta, C):
            n_sc += 1
        cl, opp, sam = closed_under_reflection(beta, dec, C)
        if cl:
            n_closed += 1
            if opp and not sam:
                n_opp += 1
    print("")
    print("  %s: %d shapes" % (label, len(group)))
    print("    self-complementary                      : %d" % n_sc)
    print("    term set closed under T -> C - T        : %d" % n_closed)
    print("    ... and every matched pair opposite sign: %d" % n_opp)

print("")
print("  the non-self-complementary zeros, with their term count:")
shown = 0
for lam, beta, cls, C, dec, z in zeros:
    if selfcomp(beta, C):
        continue
    cl, opp, sam = closed_under_reflection(beta, dec, C)
    print("    lam=%-24s beta=%-30s terms=%-4d centres=%s closed=%s"
          % (str(lam), str(beta), len(dec), str(C), "yes" if cl else "no"))
    shown += 1
    if shown >= 20:
        break

print("")
print("DONE")
