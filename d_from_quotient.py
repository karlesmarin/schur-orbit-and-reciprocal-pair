"""Are the three d's readable off the t-QUOTIENT? (Ayyer-Kumari's own language.)
Conjecture, from the beta-set arithmetic:
    d1 = t(k_A + 1),  d2 = t(k_B + 1),  d3 = | t(|lam^(rA)| - |lam^(rB)|) + 2(rA - rB) |
with k_r = lam^(r)_1 - lam^(r)_2 the sl2 content of the r-th quotient component, and rA, rB
the two residues carrying two beta-numbers.
Authors: Carles Marin, Claude (AI assistant)."""
from law_control import partitions
from theorem_full import setup
from extra_structure import quot

def part(q, i):
    return q[i] if i < len(q) else 0

ok2 = bad2 = ok3 = bad3 = 0
firstbad = []
for t in range(2, 8):
    N = t + 2
    for n in range(0, 20):
        for lam in partitions(n, N):
            st = setup(lam, t)
            if st is None: continue
            beta, Ac, Bc = st
            a1, a2 = beta[Ac[0]], beta[Ac[1]]
            b1, b2 = beta[Bc[0]], beta[Bc[1]]
            d = (a1-a2, b1-b2, abs(a1+a2-b1-b2))
            rA, rB = a1 % t, b1 % t
            q = quot(lam, t, N)
            kA = part(q[rA], 0) - part(q[rA], 1)
            kB = part(q[rB], 0) - part(q[rB], 1)
            sA, sB = sum(q[rA]), sum(q[rB])
            pred = (t*(kA+1), t*(kB+1), abs(t*(sA-sB) + 2*(rA-rB)))
            overlap = len(set(Ac) & set(Bc)) > 0
            if pred == d:
                if overlap: ok3 += 1
                else: ok2 += 1
            else:
                if overlap: bad3 += 1
                else:
                    bad2 += 1
                    if len(firstbad) < 5: firstbad.append((t, lam, d, pred, q))
print(f"two-class profile (rA != rB):   match {ok2}   MISMATCH {bad2}")
print(f"size-3 profile  (rA == rB):     match {ok3}   MISMATCH {bad3}")
for f in firstbad: print("   ", f)
