# -*- coding: utf-8 -*-
# PROOF that |E| is even, with every step machine-checked, and the control that would make it vacuous.
#
# NOTATION.  All t residue classes occupied (branch (a) is excluded by hypothesis).  For a partition
# lambda with beta set beta, let cls_k = {i : beta_i = k mod t}, let E = {k : |cls_k| >= 2} be the
# EXCESS CLASSES, and let S = {beta_i : i in a class of E} be the excess VALUES.  The criterion is
#     (i)   sigma_C(S) = S,  where sigma_C(b) = C - b  and  C = min S + max S;
#     (ii)  both fixed points of sigma_C on Z/t lie in E.
#
# ------------------------------------------------------------------------------------------------
# THEOREM.  (i) and (ii) imply |E| is even.
#
# Step 1.  E is sigma_C-stable.
#   sigma_C acts on residues by k -> C - k, and the residue of C - b is C - k whenever b = k mod t.
#   Take k in E and b in S with b = k mod t.  By (i), C - b lies in S, so C - b is a beta value
#   sitting in an excess class, and its class is C - k.  Hence C - k is in E.  So sigma_C(E) is
#   contained in E, and sigma_C is an involution, so sigma_C(E) = E.
#
# Step 2.  sigma_C has exactly two fixed points on Z/t.
#   k is fixed iff 2k = C mod t, a linear congruence with gcd(2,t) solutions when gcd(2,t) divides C
#   and none otherwise.  Hypothesis (ii) asserts that fixed points exist and lie in E, so t is even
#   and C is even, and then there are exactly two: k = C/2 and k = C/2 + t/2.
#
# Step 3.  Conclusion.
#   E is a union of sigma_C-orbits by Step 1.  An orbit has size 1 (a fixed point) or 2.  By Step 2
#   there are exactly two fixed points in Z/t, and by (ii) both are in E.  Therefore
#        |E| = 2 + 2 * (number of two-element orbits inside E),
#   which is even.  QED.
#
# COROLLARY (unconditional, no criterion needed).  With every class occupied,
#        sum_{k in E} (|cls_k| - 1) = 2r,
#   because the t - |E| non-excess classes have exactly one element each, so
#   sum_{k in E} |cls_k| = N - (t - |E|) = (t + 2r) - t + |E| = 2r + |E|.  Hence 1 <= |E| <= 2r, and
#   with the theorem, 2 <= |E| <= 2r and even.
#
# ------------------------------------------------------------------------------------------------
# WHAT THIS DOES *NOT* PROVE.  It proves (i)+(ii) => |E| even.  Since (i)+(ii) => vanishing is also
# proved (the terms pair with ratio -1), every shape the criterion calls a zero has |E| even.  That
# EVERY zero satisfies (i)+(ii) is the open converse, so "every zero has |E| even" remains MEASURED
# (121 zeros, 0 exceptions), not proved.  The theorem is about the criterion, not about the locus.
#
# MACHINE CHECKS, each able to fail:
#   V1  Step 1 on every concentric shape: sigma_C(E) = E.
#   V2  Step 2: the number of fixed points of sigma_C is gcd(2,t) when it divides C, else 0.
#   V3  the corollary sum_{k in E}(|cls_k|-1) = 2r, on every shape with all classes occupied.
#   V4  no shape satisfies (i)+(ii) with |E| odd.
#   V5  THE CONTROL THAT WOULD MAKE IT VACUOUS: |E| odd must actually OCCUR among concentric
#       shapes, and among all shapes.  If |E| were always even the theorem would say nothing.
#
# Authors: Carles Marin, Claude (AI assistant).

TS = [3, 4, 5, 6, 8]
L = lcm(TS)
p = next_prime(10 ** 9)
while (p - 1) % L != 0:
    p = next_prime(p)
F = GF(p)
G = F.multiplicative_generator()
for t in TS:
    assert (p - 1) % t == 0
print("field GF(%d); guard on t in %s -> PASS" % (p, TS))


def zeta(t):
    z = G ** ((p - 1) // t)
    assert z ** t == 1 and all(z ** k != 1 for k in range(1, t))
    return z


def beta_of(lam, N):
    lam = list(lam) + [0] * (N - len(lam))
    return [lam[i] + N - 1 - i for i in range(N)]


print("")
print("=" * 96)
print("V1-V5  the three steps, the corollary, and the control")
print("=" * 96)
print("")
print("   t   r   N |lam|<=  shapes  concentric  V1 bad  V2 bad  V3 bad | crit  V4 bad | |E| odd")
print("  " + "-" * 92)

TOT = dict(v1=0, v2=0, v3=0, v4=0, odd_any=0, odd_conc=0, crit=0)
for t, r, MAX in ((3, 2, 22), (4, 1, 26), (4, 2, 22), (4, 3, 18),
                  (5, 2, 20), (6, 2, 18), (6, 3, 16), (8, 2, 16)):
    N = t + 2 * r
    nsh = ncon = v1 = v2 = v3 = v4 = ncrit = oddany = oddconc = 0
    for size in range(MAX + 1):
        for l in Partitions(size, max_length=N):
            beta = beta_of(list(l), N)
            cl = {}
            for k, b in enumerate(beta):
                cl.setdefault(b % t, []).append(k)
            if len(cl) < t:
                continue
            nsh += 1
            E = set(k for k in cl if len(cl[k]) >= 2)
            # --- V3, the unconditional corollary
            if sum(len(cl[k]) - 1 for k in E) != 2 * r:
                v3 += 1
            if len(E) % 2:
                oddany += 1
            S = set(beta[i] for k in E for i in cl[k])
            if not S:
                continue
            C = min(S) + max(S)
            if set(C - b for b in S) != S:
                continue                                     # not concentric
            ncon += 1
            if len(E) % 2:
                oddconc += 1
            # --- V1, Step 1
            if set((C - k) % t for k in E) != E:
                v1 += 1
            # --- V2, Step 2
            fixed = [k for k in range(t) if (2 * k - C) % t == 0]
            want = gcd(2, t) if (C % gcd(2, t) == 0) else 0
            if len(fixed) != want:
                v2 += 1
            # --- V4, the theorem
            if all(k in E for k in fixed) and len(fixed) == 2:
                ncrit += 1
                if len(E) % 2:
                    v4 += 1
    print("  %2d %3d %3d %6d %7d %11d %7d %7d %7d | %5d %7d | %7d"
          % (t, r, N, MAX, nsh, ncon, v1, v2, v3, ncrit, v4, oddconc))
    for a, b in (('v1', v1), ('v2', v2), ('v3', v3), ('v4', v4),
                 ('crit', ncrit), ('odd_any', oddany), ('odd_conc', oddconc)):
        TOT[a] = TOT.get(a, 0) + b

print("")
print("  V1/V2/V3/V4 must all be 0.  The last column is V5: shapes that ARE concentric and have")
print("  an ODD number of excess classes.  It must be NONZERO, or the theorem forbids nothing.")
print("")
print("  totals: V1=%d  V2=%d  V3=%d  V4=%d   |   shapes meeting (i)+(ii): %d"
      % (TOT['v1'], TOT['v2'], TOT['v3'], TOT['v4'], TOT['crit']))
print("          |E| odd among ALL shapes: %d      among CONCENTRIC shapes: %d"
      % (TOT['odd_any'], TOT['odd_conc']))
print("")
if TOT['v1'] == TOT['v2'] == TOT['v3'] == TOT['v4'] == 0 and TOT['odd_conc'] > 0:
    print("  PROVED AND NON-VACUOUS: every step holds, no shape meets (i)+(ii) with |E| odd, and")
    print("  |E| odd does occur among concentric shapes -- so the theorem excludes real cases.")
else:
    print("  SOMETHING FAILED -- read the columns.")
print("")
print("DONE")
