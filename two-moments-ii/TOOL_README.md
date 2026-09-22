# 🧮 Part II — every exponent: the classification, the column $n=3$, and where it stops being PORC

**Thirteen scripts behind the statements marked *Measured*, each with its archived output, and
the data behind every figure and the table of certificates. Python 3 with `numpy` and `sympy`,
exact integers wherever a count is exact.**

The note counts the $n$-subsets $T\subseteq\mathbb F_p$ with

$$\sum_{t\in T}t\;=\;\sum_{t\in T}t^{k}\;=\;0$$

for every $k\ge3$: it classifies when that count is elementary, decides the column $n=3$ by the
Cauchy–Liouville–Mirimanoff polynomials $Q_k$, and shows where the count stops being PORC. Every
reference `gates/…` in the note is a file in this package, at the same path.

## ⏱️ Thirty seconds

```
cd gates
python P3d_familia_Qk.py
```

The family $Q_k$ in one run: Waring's formula against the recursion, $\dim M_{2k}=\deg Q_k+1$
checked two ways, the product expansion of $\Delta$, the cubic of $k=11$, and the genus-two curve
$D$ of the genus-five Jacobian.

Every script has its run saved beside it as `*_OUT.txt`, so any number the note calls *measured*
can be located without running anything.

## 📐 What is measured, by section of the note

| section | scripts |
|---|---|
| **The criterion, and why only $\gcd(k-1,6)$** | `P3d_tricotomia_Sm.py` — the trichotomy of $S_m$ by exhaustive enumeration over $\mu_m^4$, its collapse to $\gcd(m,6)$, the counterexample $(3,1,1,1)$; `P3d_mann_minimalidad.py` — the minimality hypothesis of Mann's theorem and the $2+2$ step it requires |
| **The counts that survive** | `literatura_super_vandermonde.py` — the case $k=4$, $n=3$ is a super-Vandermonde set of size three |
| **Two different failures** | `P3d_columna_n3.py` — $C_k(3,p)$ orbit by orbit; `P3d_caracter_regular_S3.py` — $6N_{\mathcal O}(p)=\chi(\mathrm{Frob}_p)$ at $k=6$, the left side from a direct count of $C_6(3,p)$ and the right from the factorisation of the cubic, so neither is computed from the other; `P3d_familia_Qk.py`, `P3d_grado_Qk.py` — the family $Q_k$ and its degree; `P3d_grupo_del_caracter.py` — the group of the character for $3\le k\le14$; `P3b_pofs_vs_porc.py` — the two failures are of different kinds |
| **What is left open** | `P3d_weil_D.py` — the Weil polynomials of $D$ from $\#D(\mathbb F_p)$ **and** $\#D(\mathbb F_{p^2})$, ordinarity, the four Howe–Zhu comparisons; `P3d_jacobiana_genero5.py` — the irreducibility certificates and the two negative controls, the Weil restriction $T^4-T^2+121$ and the CM curve $v^2=t^5-1$ |
| **Appendix A: proof of Proposition 13** | `P3d_s3_certificados.py` — the forty certificates for $6\le k\le46$, one printed prime each, and the table that prints them; `P3d_lc_Qk.py` — $\mathrm{lc}(Q_k)\mid k$ up to $k=400$, the finite exceptional set, and when $Q_k$ is monic |

All scripts live in `gates/` and run from there.

## 🗂️ Data

- `datos_P3d.json` — what every figure draws: $\deg Q_k$, the coefficients of $Q_k$, and the
  values of $6C_k(3,p)/(p-1)$ with the primes that give them. The figures read it; they do not
  compute. `python datos_P3d.py` regenerates it.
- `gates/tabla_certificados.tex` and `gates/tabla_certificados_pie.tex` — the table of certificates
  and its caption, exactly as the note inputs them, written by `P3d_s3_certificados.py`. Each row
  is one irreducibility check modulo one printed prime, verifiable in any computer algebra system.

## 🧾 What this package is, and is not

It is the material behind the numbers, kept so that nothing rests on a session that no longer
exists. It is **not** a library: each script is standalone, prints its own verdict, and exits
non-zero when it fails. Where a check could pass vacuously, the script carries a negative
control, a deliberately wrong input that has to fail.

Part I, the pair $\{1,3\}$, is [10.5281/zenodo.22886540](https://doi.org/10.5281/zenodo.22886540).

## 📄 Licence

CC0 1.0 — public domain. See `LICENSE`.

Carles Marín · <karlesmarin@gmail.com> · [ORCID 0009-0007-5637-9688](https://orcid.org/0009-0007-5637-9688)
