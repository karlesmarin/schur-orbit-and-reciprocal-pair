# -*- coding: utf-8 -*-
# WHY THE TOP-DEGREE PART MISSES EXACTLY 20 SHAPES: the stratum is a smaller copy of the problem.
#
# *** TWO CORRECTIONS OF MINE, BOTH FROM THIS FILE'S OWN OUTPUT ***
# (1) The T3 equality below compared a STRING against a LIST, so it printed False on a case that is
#     true and, worse, was not in the pass condition.  A guard that cannot say yes is not a guard.
#     Fixed; T3 is now fatal on it.
# (2) The reading in the next paragraph is TOO STRONG and the contingency table refutes it.
#     Closure of the maximiser set under T -> C - T is NECESSARY for [Phi]_top = 0 -- 0 shapes in
#     12613 vanish without it -- but it is very far from sufficient: 4438 shapes are closed and do
#     NOT vanish.  Closure is the analogue of concentricity; what is still missing at stratum level
#     is the analogue of condition (ii), i.e. the SIGN.  See stratum_sign.sage.
#
# middle_block.sage proved that a unique degree-maximiser forces Phi_t != 0, and measured that the
# whole top-degree part decides 12404 of the 12424 criterion-failing shapes, missing 20 -- every one
# of them with (i) FAILING and (ii) HOLDING, and every one with exactly TWO maximisers.
#
# THE READING TO TEST.  The proved vanishing direction is: when S is concentric about C and both
# fixed classes are excess, the terms pair under T -> C - T with ratio -1 and the sum collapses.
# The 20 misses look like the same mechanism acting on a PROPER SUBSET: S is not concentric, but the
# handful of terms that attain the top degree is closed under T -> C - T all the same.  If that is
# right, the top-degree part is not a weaker tool that happens to fail -- it is the SAME criterion
# applied to a sub-configuration, and the converse needs a filtration, one stratum at a time.
#
# Worked by hand, beta = [10,9,7,4,3,2,1,0], t = 4, r = 2, so S = beta, n = 8, e = 4, C = 0+10 = 10.
# The classes are 0:{4,0}, 1:{9,1}, 2:{10,2}, 3:{7,3}; the middle block S[2:6] = {7,4,3,2} has class
# 3 twice, so it is not a transversal and the a priori bound 19-1 = 18 is not attained.  The maximum
# is 16 and exactly two transversals attain it:
#       g = {4,1,2,7}  ->  T = (10,9,3,0),  H = {10,9}, L = {3,0}
#       g = {4,9,2,3}  ->  T = (10,7,1,0),  H = {10,7}, L = {1,0}
# and 10 - T1 = (10,7,1,0) = T2.  Their high-minus-low differences agree as multisets --
#       {10-3, 9-0} = {7,9}  and  {10-0, 9-3} = {10,6}   against   {10-1, 7-0} = {9,7} and
#       {10-0, 7-1} = {10,6}
# -- so the two products a_H(z) a_L(1/z) are the same up to sign, and with opposite w they cancel.
# S itself is NOT concentric: 10 - 4 = 6 is not in S.
#
# TESTED, each able to fail:
#   T1  is the maximiser set closed under T -> C - T, with C = min S + max S?  Correlate with
#       [Phi]_top = 0 over ALL shapes.  The contingency table is the test: closure and vanishing
#       must agree, or the reading is wrong.
#   T2  a DECOY centre: closure under T -> C' - T with C' = min(union of maximisers) + max(same).
#       If the decoy scores as well as C, then C is not doing the work.
#   T3  the hand case above must come out of the code exactly as written, printed in full.
#   T4  how much of "Phi = 0 => (ii)" is already PROVED: among the shapes where (ii) fails, how
#       many have a unique degree maximiser (settled by the corollary) and how many only survive
#       through the full top-degree part (still measured).
#
# CONTROLS:
#   K1  non-vacuity: closure must occur AND fail among criterion-failing shapes.
#   K2  forced: on criterion-holding shapes closure must always hold; a miss kills the reading.
#
# Authors: Carles Marin, Claude (AI assistant).

import itertools

CONF = [(4, 1, 30), (4, 2, 24), (4, 3, 18), (6, 2, 18), (6, 3, 14), (8, 2, 16)]


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


def perm_sign(q):
    n = len(q)
    seen = [False] * n
    s = 1
    for i in range(n):
        if seen[i]:
            continue
        j, ln = i, 0
        while not seen[j]:
            seen[j] = True
            j = q[j]
            ln += 1
        if ln % 2 == 0:
            s = -s
    return s


def topdeg_dict(T, r):
    D = {}
    n = 2 * r
    for a in itertools.permutations(range(r)):
        for b in itertools.permutations(range(r)):
            q = [0] * n
            e = [0] * r
            for i in range(r):
                q[i] = 2 * a[i]
                e[a[i]] += T[i]
            for i in range(r):
                q[r + i] = 2 * b[i] + 1
                e[b[i]] -= T[r + i]
            k = tuple(e)
            D[k] = D.get(k, 0) + perm_sign(q)
    return {k: v for k, v in D.items() if v != 0}


def analyse(beta, t, r):
    N = len(beta)
    cl = {}
    for i, b in enumerate(beta):
        cl.setdefault(b % t, []).append(i)
    if len(cl) < t:
        return None
    E = sorted(k for k in cl if len(cl[k]) >= 2)
    S = sorted((beta[i] for k in E for i in cl[k]), reverse=True)
    C = S[0] + S[-1]
    conc = set(C - b for b in S) == set(S)
    fixed = [k for k in range(t) if (2 * k - C) % t == 0]
    cond_ii = (len(fixed) == 2 and all(k in E for k in fixed))
    tm = []
    for pick in itertools.product(*[cl[k] for k in sorted(cl)]):
        P = sorted(pick)
        w = perm_sign([beta[i] % t for i in P])
        if sum(P) % 2:
            w = -w
        Ps = set(P)
        tm.append((w, tuple(beta[i] for i in range(N) if i not in Ps),
                   tuple(sorted(beta[i] for i in P))))
    degs = [sum(T[:r]) - sum(T[r:]) for _, T, _ in tm]
    Dmax = max(degs)
    G = [tm[i] for i, d in enumerate(degs) if d == Dmax]
    top = {}
    for w, T, _ in G:
        for k, v in topdeg_dict(T, r).items():
            top[k] = top.get(k, 0) + w * v
    top = {k: v for k, v in top.items() if v != 0}
    Ts = set(T for _, T, _ in G)
    closedC = all(tuple(sorted((C - x for x in T), reverse=True)) in Ts for T in Ts)
    vals = [x for T in Ts for x in T]
    Cp = min(vals) + max(vals)
    closedCp = all(tuple(sorted((Cp - x for x in T), reverse=True)) in Ts for T in Ts)
    return dict(crit=conc and cond_ii, conc=conc, cond_ii=cond_ii, S=S, C=C, Cp=Cp,
                Dmax=Dmax, G=G, nG=len(G), top=top, closedC=closedC, closedCp=closedCp, e=len(E))


# ---------------------------------------------------------------- T3, the hand case --------------
print("=" * 104)
print("T3  the hand case, printed by the code: beta = [10,9,7,4,3,2,1,0], t = 4, r = 2")
print("=" * 104)
hb = [10, 9, 7, 4, 3, 2, 1, 0]
ha = analyse(hb, 4, 2)
print("")
print("   S = %s   C = %d   concentric: %s   (ii): %s"
      % (ha['S'], ha['C'], ha['conc'], ha['cond_ii']))
print("   Dmax = %d, attained by %d transversals:" % (ha['Dmax'], ha['nG']))
for w, T, g in ha['G']:
    print("      kept g = %-18s T = %-18s H = %-12s L = %-12s w = %+d"
          % (str(list(g)), str(list(T)), str(list(T[:2])), str(list(T[2:])), w))
    print("         a_H(z) a_L(1/z) monomials: %s"
          % (sorted(topdeg_dict(list(T), 2).items(), reverse=True),))
T1 = ha['G'][0][1]
refl = tuple(sorted((ha['C'] - x for x in T1), reverse=True))
same = (refl == ha['G'][1][1])
print("   C - T1 = %s   equals T2: %s" % (str(list(refl)), same))
print("   the two w are opposite: %s" % (ha['G'][0][0] == -ha['G'][1][0]))
print("   [Phi]_top = %s  ->  %s" % (ha['top'], "VANISHES" if not ha['top'] else "survives"))
if ha['top'] or ha['conc'] or not ha['cond_ii'] or ha['nG'] != 2 or not same \
        or ha['G'][0][0] != -ha['G'][1][0]:
    print("")
    print("   T3 FAILED -- the hand case is not what the header says.")
    raise SystemExit(1)
print("")
print("   T3 PASS: the case comes out exactly as written by hand.")

# ---------------------------------------------------------------- T1, T2, K1, K2 -----------------
print("")
print("=" * 104)
print("T1/T2  closure of the maximiser set under T -> C - T, against the vanishing of [Phi]_top")
print("=" * 104)
print("")
print("     t   r   shapes | crit: closed  top=0 | non-crit: closed & top=0   closed only"
      "   top=0 only  neither")
print("  " + "-" * 100)

K2bad = 0
CT = {}
CTp = {}
for t, r, MAX in CONF:
    N = t + 2 * r
    nsh = ccl = ctop = a11 = a10 = a01 = a00 = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None:
                continue
            nsh += 1
            key = (a['closedC'], not a['top'])
            CT[key] = CT.get(key, 0) + 1
            CTp[(a['closedCp'], not a['top'])] = CTp.get((a['closedCp'], not a['top']), 0) + 1
            if a['crit']:
                ccl += 1 if a['closedC'] else 0
                ctop += 1 if not a['top'] else 0
                if not a['closedC']:
                    K2bad += 1
            else:
                if a['closedC'] and not a['top']:
                    a11 += 1
                elif a['closedC']:
                    a10 += 1
                elif not a['top']:
                    a01 += 1
                else:
                    a00 += 1
    print("  %4d %3d %8d | %11d %6d | %22d %13d %12d %8d"
          % (t, r, nsh, ccl, ctop, a11, a10, a01, a00))

print("")
print("  contingency over ALL shapes, centre C = min S + max S:")
print("     closed=yes top=0: %6d      closed=yes top!=0: %6d"
      % (CT.get((True, True), 0), CT.get((True, False), 0)))
print("     closed=no  top=0: %6d      closed=no  top!=0: %6d"
      % (CT.get((False, True), 0), CT.get((False, False), 0)))
print("")
print("  T2 the DECOY centre C' = min + max of the maximiser values only:")
print("     closed=yes top=0: %6d      closed=yes top!=0: %6d"
      % (CTp.get((True, True), 0), CTp.get((True, False), 0)))
print("     closed=no  top=0: %6d      closed=no  top!=0: %6d"
      % (CTp.get((False, True), 0), CTp.get((False, False), 0)))
print("")
print("  K2 criterion-holding shapes where closure FAILS: %d (must be 0)" % K2bad)
print("  K1 non-vacuity: closure both occurs and fails among non-criterion shapes: %s"
      % ("PASS" if CT.get((True, False), 0) + CT.get((True, True), 0) > 0
         and CT.get((False, False), 0) > 0 else "VACUOUS"))

# ---------------------------------------------------------------- T4 -----------------------------
print("")
print("=" * 104)
print("T4  how much of  Phi_t = 0  =>  (ii)  is PROVED, and how much is still only measured")
print("=" * 104)
print("")
print("     t   r  (ii) fails | unique maximiser (PROVED)   only via [Phi]_top (measured)   left")
print("  " + "-" * 100)
tp = tm_ = tl = 0
for t, r, MAX in CONF:
    N = t + 2 * r
    n = pr = me = lf = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            a = analyse(beta, t, r)
            if a is None or a['cond_ii']:
                continue
            n += 1
            if a['nG'] == 1:
                pr += 1
            elif a['top']:
                me += 1
            else:
                lf += 1
    print("  %4d %3d %11d %26d %30d %6d" % (t, r, n, pr, me, lf))
    tp += pr
    tm_ += me
    tl += lf
print("")
print("  totals: (ii) fails on %d shapes.  PROVED nonzero: %d (%.1f%%).  Measured-only: %d.  Left: %d."
      % (tp + tm_ + tl, tp, 100.0 * tp / (tp + tm_ + tl), tm_, tl))
print("")
print("DONE")
