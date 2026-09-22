# 🔢 Subsets of F_p with two prescribed moments — code and data (Part I)

Companion material for *Subsets of $\mathbb F_p$ with two prescribed moments, Part I: the pair
$\{1,3\}$* (Carles Marín), deposited on Zenodo:
[doi:10.5281/zenodo.22886540](https://doi.org/10.5281/zenodo.22886540). That is the **concept
DOI** — it always resolves to the current version — and it carries the note in both editions,
English and Spanish.

Both PDFs are here too — `P3c.pdf` and `P3c_es.pdf`, byte for byte the files deposited on Zenodo —
together with the program `P3c.py`, which recomputes every count the note quotes, and the ten
scripts behind the statements the note calls *measured*, each with its archived output.

This is **Part I**; Part II, which takes every exponent, is
[`two-moments-ii/`](../two-moments-ii). The count studied here first appears in
[`s-matrix-column-ii/`](../s-matrix-column-ii), as $C_{\{1,3\}}(n,p)$ with its main term and an
error bound, and there its closed form is left open. This note takes it up.

## 🎯 The question

Fix an odd prime $p$ and let $C(n,p)$ count the $n$-subsets $T\subseteq\mathbb F_p$ with

$$\sum_{t\in T}t\;=\;\sum_{t\in T}t^{3}\;=\;0,$$

equivalently, for $p>3$, the totally split squarefree polynomials of degree $n$ whose coefficients
of $X^{n-1}$ and $X^{n-3}$ both vanish. For which $n$ does this count have an **elementary**
expression — a polynomial in $p$ with coefficients in quadratic characters of fixed integers?

The partition sieve attaches to each $\lambda\vdash n$ a weighted hyperplane section of a
diagonal cubic. Its singular points are the sign vectors that split the weights into two blocks of
equal sum — the *biquanimous* splittings of $\lambda$ — and the stratum of length exactly four
decides everything: there the variety is a plane cubic, and a smooth one is an elliptic curve,
whose Frobenius trace no elementary expression can carry.

## 📐 What the note proves

1. **The classification.** For $n\ge3$, $C(n,p)$ is elementary precisely for $n=3,4,6$. The
   uniform witness $\lambda=(n-3,1,1,1)$ settles every $n>6$ at once, since $n-3>3$ leaves no
   equal-sum splitting; $n=5$ falls because its total is odd.
2. **Four closed forms.** $C(3,p)$, $C(4,p)$, $C(5,p)$ and $C(6,p)$ exactly; the case $n=5$
   involves the Frobenius trace of the elliptic curve of conductor $20$, so it is exact without
   being elementary — the two are different things, and it is the second that is classified.
3. **Two strata worked out by pencils of conics.** $V_{(a,1,1,1,1)}$ has line field
   $\mathbb Q(\sqrt{D_a},\sqrt{E_a})$, abelian of exponent two for every $a$, while
   $V_{(2,2,1,1,1)}$ — also with exactly two distinct weights — does not: its line field contains
   the splitting field of a cubic of discriminant $-140$ and Galois group $S_3$.
4. **The extremal cubics, over the base field.** The four rational nodes of $V_{(2,1^4)}$ and the
   matrix carrying it to Cayley's standard form; $V_{(1^6)}$ as the Segre cubic on the nose.

## 🧾 What is measured and what is not

The theorems are proved. The blocks marked *measured* check them against exhaustive enumeration
over $\mathbb F_p$, each naming its range and its script, and each script carries a negative
control where a check could otherwise pass vacuously. `P3c.py` recomputes every count the note
quotes in one run and exits non-zero if any printed value fails; its run is archived as
`P3c_selftest.txt`.

## 🧰 What runs what

`python P3c.py` is the thirty-second check. The ten scripts live under `gates/`, at the same paths
the note prints; eight are plain Python with `numpy` and `sympy`, and two need SageMath. See
`TOOL_README.md` for the map from each section of the note to the scripts behind it.

## 📎 Licence and citation

The scripts are CC0; the note is CC BY 4.0. Cite the concept DOI
[10.5281/zenodo.22886540](https://doi.org/10.5281/zenodo.22886540).

Carles Marín, independent researcher — `karlesmarin@gmail.com` —
[ORCID 0009-0007-5637-9688](https://orcid.org/0009-0007-5637-9688).
