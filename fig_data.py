"""Data for the image-of-the-map figure. Authors: Carles Marin, Claude (AI assistant)."""
from collections import defaultdict
from law_control import partitions
from theorem_full import setup, lambda11

def image(t, maxsize):
    N = t + 2
    pts = defaultdict(lambda: {"n": 0, "signs": set()})
    nz = nlam = 0
    for n in range(0, maxsize + 1):
        for lam in partitions(n, N):
            nlam += 1
            st = setup(lam, t)
            if st is None: continue
            beta, Ac, Bc = st
            a1, a2 = beta[Ac[0]], beta[Ac[1]]
            b1, b2 = beta[Bc[0]], beta[Bc[1]]
            d = (a1 - a2, b1 - b2, abs(a1 + a2 - b1 - b2))
            if d[2] == 0: continue
            l11 = lambda11(beta, Ac, Bc, N)
            eps = (-1)**(t + (t+2)*(t+3)//2) * l11 * (1 if a1+a2-b1-b2 > 0 else -1)
            key = tuple(sorted(d[:2])) + (d[2],)
            pts[key]["n"] += 1; pts[key]["signs"].add(eps)
            nz += 1
    return pts, nlam, nz

if __name__ == "__main__":
  for t in (3, 5):
      for M in (18,):
          pts, nlam, nz = image(t, M)
          mixed = sum(1 for v in pts.values() if len(v["signs"]) > 1)
          d1s = sorted({k[0] for k in pts}); d3s = sorted({k[2] for k in pts})
          print(f"t={t} |lam|<={M}: {nlam} partitions -> {nz} nonzero -> {len(pts)} lattice points"
                f"   mixed-sign points {mixed}")
          print(f"   d1 values: {d1s[:8]}...   d3 values: {d3s[:10]}...")
          print(f"   max multiplicity {max(v['n'] for v in pts.values())}")
