"""What IS the extra independence family? core, quotient, and the d-triple.
Authors: Carles Marin, Claude (AI assistant)."""
from mpmath import mp, mpf
from theorem_full import setup
from ak53_consistency import core
mp.dps = 25

def quot(lam, t, n):
    """t-quotient from the beta-set (Macdonald / AK section 2.1)."""
    lam = list(lam) + [0]*(n - len(lam))
    beta = [lam[i] + n - 1 - i for i in range(n)]
    out = []
    for r in range(t):
        v = sorted([b for b in beta if b % t == r], reverse=True)
        tv = [(b - r)//t for b in v]
        out.append(tuple(p for p in [tv[j] - len(v) + 1 + j for j in range(len(v))] if p > 0))
    return tuple(out)

def dtriple(lam, t):
    st = setup(lam, t)
    if st is None: return None
    beta, Ac, Bc = st
    a1, a2 = beta[Ac[0]], beta[Ac[1]]
    b1, b2 = beta[Bc[0]], beta[Bc[1]]
    return (a1-a2, b1-b2, abs(a1+a2-b1-b2))

print("THE EXTRA FAMILY  lam = (l2 + 3t/2 - 1, l2),  t/2 <= l2 <= t-1,  t even")
print(f"{'t':>3} {'lambda':>10} {'core_t':>10} {'quot_t':>26} {'(d1,d2,d3)':>14}  two d = t?")
for t in (2, 4, 6, 8):
    for l2 in range(t//2, t):
        lam = (l2 + 3*t//2 - 1, l2)
        d = dtriple(lam, t)
        q = quot(lam, t, t+2)
        qs = "(" + ",".join(str(list(x)) for x in q) + ")"
        two = sorted(d)[:2] == [t, t]
        print(f"{t:>3} {str(lam):>10} {str(core(lam,t,t+2)):>10} {qs:>26} {str(d):>14}  {two}")

print("\nGENERAL: for which lambda (ANY length <= t+2) do two of the three d's equal t?")
from law_control import partitions
for t in (2, 3, 4, 5, 6):
    hit = []
    for n in range(0, 26):
        for lam in partitions(n, t+2):
            d = dtriple(lam, t)
            if d and d[2] != 0 and sorted(d)[:2] == [t, t]:
                hit.append((lam, d))
    lens = sorted({len(h[0]) for h in hit})
    print(f"  t={t}: {len(hit)} such lambda, lengths present {lens};  first: {hit[:4]}")
