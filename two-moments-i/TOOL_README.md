# 🧮 Part I — the pair {1,3}, and exactly when the count is elementary

**One program that recomputes every count the note quotes, and ten scripts behind the statements
marked *Measured*, each with its archived output. Python 3 with `numpy` and `sympy`; two scripts
need SageMath.**

The note counts the $n$-subsets $T\subseteq\mathbb F_p$ with

$$\sum_{t\in T}t\;=\;\sum_{t\in T}t^{3}\;=\;0,$$

that is, the totally split squarefree polynomials with a non-consecutive pair of prescribed
coefficients, and shows the count is elementary — a polynomial in $p$ with coefficients in
quadratic characters of fixed conductor — exactly for $n=3,4,6$. Every reference `gates/…` in the
note is a file in this package, at the same path.

## ⏱️ Thirty seconds

```
python P3c.py
```

It computes $C(n,p)$ by dynamic programming over (size, sum, sum of cubes), assembles the sieve
from the strata, and checks every closed form, every row of both tables, both propositions of the
closed-form section and the classification against what the note prints. It exits non-zero if
any of them fails. Its run is archived as `P3c_selftest.txt`.

Every script has its run saved beside it as `*_OUT.txt`, so any number the note calls *measured*
can be located without running anything.

## 📐 What is measured, by section of the note

| section | scripts |
|---|---|
| **The singular locus, and the count it produces** | `gates/P3b_singularidades.py` — the singular locus of $V_\lambda$ is the set of sign vectors splitting the weights into two equal-sum blocks, checked case by case |
| **When the count is elementary** | `gates/P3c_superficie_22111.py` — the surface $V_{(2,2,1,1,1)}$, its count through $\chi_p(5)$, and why "the cubic splits" is not "the conics split into rational lines" ($p=71,73,83,149$); `gates/P3c_haz_conicas_22111.py` — its conic bundle, derived step by step; `gates/P3c_familia_a11111.py` — the family $V_{(a,1,1,1,1)}$ and the hypotheses on $p$ it needs; `gates/P3b_haz_conicas_Sa.sage` — the conic bundle of the cubic surface $S_a$, symbolically over $\mathbb Q$ and by fibre counts over $\mathbb F_p$ |
| **The four closed forms** | `gates/P3b_B5_cerrada.py`, `gates/P3b_B6_cerrada.py` — $C(5,p)$ and $C(6,p)$ against the exact census in the $33$ primes $7\le p\le151$, with the traces of the curve of conductor $20$ against the eta product; `gates/P3b_B5_weierstrass.sage` — the plane cubic of $\lambda=(2,1,1,1)$ carried to its minimal model over $\mathbb Q$, with conductor and discriminant |
| **Computations** | `P3c.py` — everything above, in one run |
| **The objects in the tables** (appendix) | `gates/P3c_tres_cubicas.py` — the three cubics of the tables are the three classical ones; `gates/P3c_cayley_segre_sobre_Q.py` — the Segre and Cayley models over the base field, not only over the closure |

The two Sage scripts run as `sage <name>.sage` from `gates/`; their verdict is their last line.

## 🧾 What this package is, and is not

It is the material behind the numbers, kept so that nothing rests on a session that no longer
exists. It is **not** a library: each script is standalone, prints its own verdict, and — the
Python ones — exits non-zero when it fails. Where a check could pass vacuously, the script carries
a negative control, a deliberately wrong input that has to fail. Nothing needs network access or
an external database: the conductor is obtained by Tate's algorithm, and the eta product is an
independent check on it.

Part II, which takes every exponent $k\ge3$, is
[10.5281/zenodo.22886542](https://doi.org/10.5281/zenodo.22886542).

## 📄 Licence

CC0 1.0 — public domain. See `LICENSE`.

Carles Marín · <karlesmarin@gmail.com> · [ORCID 0009-0007-5637-9688](https://orcid.org/0009-0007-5637-9688)
