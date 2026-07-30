"""Is  Phi = a single sl2 character chi_k  <=>  two of the d_i equal t and the third is 2(k+1)?
The cyclotomic argument, checked. Authors: Carles Marin, Claude (AI assistant)."""
from mpmath import mp, mpf, sinh
from law_control import partitions
from theorem_full import setup, closed_form
mp.dps = 30
theta = mpf("0.37")

ok = bad = 0
for t in range(2, 8):
    N = t + 2
    for n in range(0, 17):
        for lam in partitions(n, N):
            st = setup(lam, t)
            if st is None: continue
            beta, Ac, Bc = st
            a1, a2 = beta[Ac[0]], beta[Ac[1]]
            b1, b2 = beta[Bc[0]], beta[Bc[1]]
            d = (a1-a2, b1-b2, abs(a1+a2-b1-b2))
            if d[2] == 0: continue
            v = closed_form(lam, t, theta)
            # is |v| a single sl2 character chi_k = sinh((k+1)theta)/sinh(theta) for some k >= 0?
            is_single = False; kfound = None
            for k in range(0, 4*sum(d) + 4):
                if abs(abs(v) - sinh((k+1)*theta)/sinh(theta)) < mpf("1e-18"):
                    is_single = True; kfound = k; break
            ds = list(d)
            pred = ds.count(t) >= 2
            if pred and is_single:
                rest = ds[:]; rest.remove(t); rest.remove(t)
                third = rest[0]
                ok += (third == 2*(kfound+1))
                bad += (third != 2*(kfound+1))
            elif pred == is_single:
                ok += 1
            else:
                bad += 1
                if bad <= 5: print("  MISMATCH", t, lam, d, "single:", is_single, "k:", kfound)
print(f"criterion 'two d_i = t' <=> 'Phi is a single chi_k, third d = 2(k+1)':  ok {ok}   bad {bad}")
