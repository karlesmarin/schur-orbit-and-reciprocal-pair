"""Does the theorem give a product formula for a ROOT-OF-UNITY-weighted tableau count
with a free parameter still alive?  (The enumerative payoff, Ayyer-Behrend style.)

s_lam(mu_t, z, 1/z) = SUM over SSYT(lam, t+2) of  w^{sum_k k*m_k} * z^{m_{t+1} - m_{t+2}}
For t=2 the weight is (-1)^{#2}: a (-1)-enumeration with z FREE.
For lam rectangular, SSYT(lam, N) <-> plane partitions in a box <-> lozenge tilings.

Authors: Carles Marin, Claude (AI assistant)."""
from mpmath import mp, mpf, mpc, exp, pi
from theorem_full import closed_form
mp.dps = 30

def interlace(top):
    """all partitions inner interlacing with top:  top_1 >= in_1 >= top_2 >= in_2 >= ..."""
    if not top: 
        yield (); return
    def rec(i, prev, acc):
        if i == len(top) - 1:
            for v in range(min(prev, top[i]), -1, -1):
                yield tuple(acc + [v])
            return
        lo, hi = top[i+1], min(prev, top[i])
        for v in range(hi, lo - 1, -1):
            yield from rec(i + 1, v, acc + [v])
    yield from rec(0, top[0], [])

def tableau_sum(lam, t, z):
    """the weighted SSYT sum, built from Gelfand-Tsetlin chains (no tableaux materialised)."""
    N = t + 2
    w = exp(2j * pi / t)
    weights = [w ** k for k in range(t)] + [z, 1 / z]
    total = mpc(0)
    def rec(level, shape, acc):
        nonlocal total
        if level == 0:
            if sum(shape) == 0: total += acc
            return
        for inner in interlace(shape):
            m = sum(shape) - sum(inner)
            rec(level - 1, inner, acc * weights[level - 1] ** m)
    rec(N, tuple(list(lam) + [0] * (N - len(lam))), mpc(1))
    return total

z = exp(mpf("0.37"))
print(f"{'t':>2} {'lambda':>12} {'SSYT weighted sum':>24} {'closed form':>24}  match")
ok = bad = 0
cases = [(2,(2,2)), (2,(3,3)), (2,(2,2,2)), (2,(3,2,1)), (2,(4,4)),
         (3,(2,2)), (3,(3,1)), (3,(2,2,2)), (3,(3,3)),
         (4,(2,2)), (4,(3,2)), (4,(2,2,1)), (5,(3,1)), (5,(2,2))]
for t, lam in cases:
    s = tableau_sum(lam, t, z)
    c = closed_form(lam, t, mpf("0.37"))
    m = abs(s - c) < mpf("1e-16")
    ok += m; bad += (not m)
    print(f"{t:>2} {str(lam):>12} {mp.nstr(s.real,14):>24} {mp.nstr(c.real,14):>24}  {m}")
print(f"\nmatched {ok}, failed {bad}")

print("\n(-1)-ENUMERATION with the free parameter, t=2, rectangular lambda = (c^r):")
print("   SSYT(lam,4) <-> plane partitions in a box;  weight (-1)^{#2} z^{#3-#4}")
for (r, c) in [(2,2),(2,3),(3,2),(2,4),(3,3)]:
    lam = tuple([c]*r)
    at1 = tableau_sum(lam, 2, mpf(1))
    print(f"   box {r}x{c}:  signed count at z=1  =  {mp.nstr(at1.real, 10)}")
