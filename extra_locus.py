"""The extra independence locus created by the reciprocal specialisation.
Authors: Carles Marin, Claude (AI assistant)."""
from mpmath import mp, mpf, sinh
from theorem_full import closed_form
from ak53_consistency import core
mp.dps = 30
theta = mpf("0.41")
extra = []
for t in range(2, 13):
    for l1 in range(0, 34):
        for l2 in range(0, l1 + 1):
            lam = tuple(x for x in (l1, l2) if x > 0)
            if core(lam, t, t + 2) == lam: continue          # covered by their theorem
            ours = closed_form(lam, t, theta)
            if abs(ours - sinh((l1 - l2 + 1) * theta) / sinh(theta)) < mpf("1e-18"):
                extra.append((t, l1, l2))
pred = [(t, l2 + 3*t//2 - 1, l2) for t in range(2, 13, 2) for l2 in range(t//2, t)
        if l2 + 3*t//2 - 1 < 34]
print(f"extra cases found: {len(extra)}")
print(f"predicted family  (t even, lam = (l2 + 3t/2 - 1, l2),  t/2 <= l2 <= t-1): {len(pred)}")
print("sets identical:", sorted(extra) == sorted(pred))
print("first few:", sorted(extra)[:8])
