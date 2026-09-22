# 🔢 Subsets of F_p with two prescribed moments — code and data (Part II)

Companion material for *Subsets of $\mathbb F_p$ with two prescribed moments, Part II: every
exponent* (Carles Marín), deposited on Zenodo:
[doi:10.5281/zenodo.22886542](https://doi.org/10.5281/zenodo.22886542). That is the **concept
DOI** — it always resolves to the current version — and it carries the note in both editions,
English and Spanish.

Both PDFs are here too — `P3d.pdf` and `P3d_es.pdf`, byte for byte the files deposited on Zenodo —
together with the thirteen scripts behind the statements the note calls *measured*, each with its
archived output, the data behind every figure, and the table of certificates.

This is **Part II**; Part I, the pair $\{1,3\}$, is [`two-moments-i/`](../two-moments-i), and this
note cites it for $k=3$. Moment subset sums of this kind first appear in
[`s-matrix-column-ii/`](../s-matrix-column-ii): the case $k=3$ is written there as
$C_{\{1,3\}}(n,p)$, with its main term and an error bound, and its closed form is left open.

## 🎯 The question

Fix $k\ge3$ and an odd prime $p$, and let $C_k(n,p)$ count the $n$-subsets
$T\subseteq\mathbb F_p$ with

$$\sum_{t\in T}t\;=\;\sum_{t\in T}t^{k}\;=\;0.$$

Part I settled $k=3$: the count is elementary — a polynomial in $p$ with coefficients in quadratic
characters of fixed integers — exactly for $n=3,4,6$. This note settles every $k$: which counts
are Artin–Tate, the column $n=3$ completely, and where the count stops being polynomial on
residue classes.

The deciding stratum is again of length four, and whether it is singular depends on $k$ only
through $\gcd(k-1,6)$, by Mann's theorem on relations between roots of unity. For plane cubics
singularity forces a genus-zero normalisation, which is why it decides when $k=3$; in general the
right condition is that the normalisation have genus zero.

## 📐 What the note proves

1. **The classification.** $C_k(n,p)$ is Artin–Tate precisely for $n\in\{3,4,6\}$ when $k=3$,
   for $n\in\{3,4\}$ when $k=5$, and for $n=3$ alone for every other $k$.
2. **The column $n=3$.** With $u=e_2^3/e_3^2$, the count is governed by the roots of an explicit
   polynomial $Q_k(u)$ — of Cauchy–Liouville–Mirimanoff type — together with the cubic
   $Y^3+uY-u$ attached to each root. It is expressible in quadratic characters if and only if
   $Q_k$ is constant, which happens exactly for $k\in\{3,4,5,7\}$.
3. **The list is a dimension.** $\deg Q_k+1=\dim M_{2k}(\mathrm{SL}_2(\mathbb Z))$, so the count
   is elementary exactly at the weights $6,8,10,14$; weight $14$ returns because $M_2=0$ after
   the first cusp form appears in weight $12$.
4. **Not PORC.** Where it is not elementary, some fibre of $Q_k$ has Galois group $S_3$, and a
   permutation-character identity shows the nonabelian part cannot cancel. For $6\le k\le46$ this
   rests on forty printed irreducibility certificates, one prime each; beyond, on a proved bound.
5. **Two different failures.** At $n=3$, $k=6$ the obstruction is a nonabelian Artin-type class
   function and the count is polynomial on Frobenius sets; for $n\ge4$ it is a surviving
   weight-one Frobenius trace, and there the count is not even POFS.

## 🧾 What is measured and what is not

Every computational claim marked *measured* names its script, and every script runs on its own
from this directory. What is proved is proved; what is measured says its range. The
irreducibility of the polynomials $Q_k$ in general is Mirimanoff's question, still open, and the
note does not assume it anywhere.

## 🧰 What runs what

`cd gates && python P3d_familia_Qk.py` is the thirty-second check: the family $Q_k$ in one run.
All thirteen scripts live in `gates/`, at the same paths the note prints, and run from there with
Python 3, `numpy` and `sympy`. `datos_P3d.json` is what every figure draws — the figures read it,
they do not compute — and `python datos_P3d.py` regenerates it. The table of certificates is in
`gates/`, exactly as the note inputs it. See `TOOL_README.md` for the map from each section of the
note to the scripts behind it.

## 📎 Licence and citation

The scripts and data are CC0; the note is CC BY 4.0. Cite the concept DOI
[10.5281/zenodo.22886542](https://doi.org/10.5281/zenodo.22886542).

Carles Marín, independent researcher — `karlesmarin@gmail.com` —
[ORCID 0009-0007-5637-9688](https://orcid.org/0009-0007-5637-9688).
