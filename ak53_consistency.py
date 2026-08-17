"""Consistency of our closed form with Ayyer-Kumari [2501.00275, Thm 5.3].
Their theorem (t, m=1, Y=(1), free vars (z,1/z), n=t+2) reads:
    s_lambda(z, 1/z, mu_t) = s_lambda(z, 1/z)   <=>   lambda = core_t(lambda),   for l(lambda) <= 2.
We check BOTH directions against our Theorem.
Authors: Carles Marin, Claude (AI assistant)."""
from mpmath import mp, mpf, sinh
from theorem_full import closed_form
mp.dps = 30

def core(lam, t, n):
    """t-core of lam via beta-set packing (n >= l(lam) slots)."""
    lam = list(lam) + [0]*(n - len(lam))
    beta = [lam[i] + n - 1 - i for i in range(n)]
    cls = {}
    for b in beta: cls.setdefault(b % t, []).append(b)
    packed = []
    for r, v in cls.items():
        packed += [r + t*k for k in range(len(v))]
    packed.sort(reverse=True)
    parts = [packed[i] - (n - 1 - i) for i in range(n)]
    return tuple(p for p in parts if p > 0)

theta = mpf("0.41")
agree_indep = agree_dep = viol = 0
for t in range(2, 9):
    for l1 in range(0, 16):
        for l2 in range(0, l1 + 1):
            lam = tuple(x for x in (l1, l2) if x > 0)
            is_core = (core(lam, t, t + 2) == lam)
            ours = closed_form(lam, t, theta)
            theirs = sinh((l1 - l2 + 1) * theta) / sinh(theta)      # s_lambda(z,1/z) = chi_{l1-l2}
            equal = abs(ours - theirs) < mpf("1e-18")
            if is_core and equal: agree_indep += 1
            elif (not is_core) and (not equal): agree_dep += 1
            else:
                viol += 1
                if viol <= 6:
                    print(f"  VIOLATION t={t} lam={lam} core={core(lam,t,t+2)} "
                          f"is_core={is_core} equal={equal}")
print(f"\nlambda IS a t-core and our value = s_lambda(z,1/z):        {agree_indep}")
print(f"lambda NOT a t-core and our value != s_lambda(z,1/z):      {agree_dep}")
print(f"violations of [AK25, Thm 5.3] (either direction):    {viol}")
