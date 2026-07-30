"""Why are the rectangular signed counts squares? The interval configuration degenerates.
Authors: Carles Marin, Claude (AI assistant)."""
from law_control import partitions
from theorem_full import setup

def dtrip(lam, t):
    st = setup(lam, t)
    if st is None: return None
    b, Ac, Bc = st
    a1,a2 = b[Ac[0]], b[Ac[1]]; b1,b2 = b[Bc[0]], b[Bc[1]]
    return (a1-a2, b1-b2, abs(a1+a2-b1-b2))

print("t=2, lambda=(c^r): the three d's, and which pair coincides")
print(f"{'r':>2} {'c':>3} {'(d1,d2,d3)':>16}  repeat?   product/8")
rep_even = rep_odd = 0; tot_even = tot_odd = 0
for r in (1,2,3):
    for c in range(1, 13):
        d = dtrip(tuple([c]*r), 2)
        if d is None: continue
        rep = len(set(d)) < 3
        if c % 2 == 0: tot_even += 1; rep_even += rep
        else:          tot_odd  += 1; rep_odd  += rep
        if c <= 6:
            print(f"{r:>2} {c:>3} {str(d):>16}  {'yes' if rep else 'no ':>5}   "
                  f"{d[0]*d[1]*d[2]//8 if d[2] else 0}")
print(f"\n c EVEN: {rep_even}/{tot_even} have two of the three d's equal")
print(f" c ODD : {rep_odd}/{tot_odd} have two equal")

print("\nGeneral t: does a rectangle always degenerate?  (two d's equal, or d3=0)")
for t in range(2, 8):
    deg = tot = 0
    for r in range(1, t+3):
        for c in range(1, 15):
            d = dtrip(tuple([c]*r), t)
            if d is None: continue
            tot += 1
            deg += (len(set(d)) < 3) or d[2] == 0
    print(f"  t={t}: {deg}/{tot} rectangles degenerate")
