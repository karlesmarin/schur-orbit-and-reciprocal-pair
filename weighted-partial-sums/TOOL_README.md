# 🧮 `segmento_II.py` — zeros of weighted partial sums of completely multiplicative functions, computed

**One file. Python 3 standard library. No Sage, no numpy, nothing to install. Exact integers
throughout.**

A completely multiplicative function $\varepsilon$ with values $\pm1$ can make
$S_\varepsilon(m)=\sum_{n\le m}n\,\varepsilon(n)$ vanish only when $m\equiv0,3\pmod 4$, and the note
proves that some $\varepsilon$ does for every such $m$. This tool exhibits those zeros, re-checks the
certificates of the proof, and computes the zeros that Legendre symbols carry through a digit
recursion — each routine named after the statement it implements.

## ⏱️ Thirty seconds

```
python segmento_II.py autotest
```

It recomputes, from nothing, the worked examples of the note — explicit zeros and the counts
$a(m)$, a sample of the 104 block certificates, the digit sets $R_q$ and the zeros they propagate,
the loops at $q=59,131,347$, the weighted half sum on all 1510 odd characters of odd conductor
$f\le121$, the null cells of Conrey's function at $q=103$ — runs the decoys that must fail,
exercises every guard, and ends with

```
SELF-TEST: all correct
```

If it says anything else, the note and the code disagree, and the note is the one to distrust. The
archived run is in `autotest_OUT.txt` (about 12 seconds).

## 🧮 What it does

| command | what it answers | statement of the note |
|---|---|---|
| `cero m [full]` | a completely multiplicative $\varepsilon$ with $S_\varepsilon(m)=0$, the sum recomputed $n$ by $n$; for $m\equiv1,2\pmod4$, why none exists | Theorem *every admissible length* |
| `certificado [j]` | re-verifies the 104 interval certificates of $[5000,10^8)$, or block $j$, against `II_certificados_A_bloques.csv` | its proof, Lemmas *second moment* and *intervals of subset sums* |
| `digitos q [s] [X]` | $R_q$, $h(-q)$, and the zeros of $S_\chi$ up to $X$, each by direct summation **and** by the recursion, plus every part of the corollary that applies | Theorem *a digit recursion*, Corollary *zeros from digits* |
| `semisuma f` | the weighted half sum of **every** odd character mod $f$ against its closed form, exactly in $\mathbb Z[\zeta_N]$ | Theorem *the weighted half sum*, Corollaries *how many components vanish* and *Dirichlet's family* |
| `conrey q` | the cells on which $F_q$, hence Conrey's sine series of $(\cdot\mid q)$, vanishes identically, compared with $R_q$ | §*the same two sums in Conrey's sine series*, Proposition *positivity empties the digit set* |

As a library:

```python
from segmento_II import cero, certificado, digitos, bloque, semisuma, conrey

cero(23)[0]                 # {2: -1, 3: -1, 5: 1, ..., 13: -1, ...}   S = 0, recomputed
certificado(1)["ratio"]     # Fraction(229805468750, 516488006241)    B*(b)^2/D^2 <= 0.4449
digitos(23, 1, 10**6)["ceros"]   # [11, 264, 6083, 139920]            the zeros (23^j - 1)/2
bloque(59, (1, 22, 21, 0))  # (True, 283200, [(16708800, 0), ...])    a loop at h(-59) = 3
conrey(103)["nulas"]        # [47, 51, 55]                           = R_103
```

## 🚧 What it refuses to do

- `cero` builds explicit zeros up to $m=10^6$ (about a second). Above that it **declines and computes
  nothing**: a zero exists there by the theorem — the certificates up to $10^8$ and the analytic
  tail beyond — but the tool does not claim a witness it has not produced. Every witness it does
  return has had its sum recomputed from the signs.
- `digitos` declines when $q$ is not an odd prime, and sums directly only up to $10^7$; beyond, the
  recursion `S_recursion` evaluates single values.
- `semisuma` declines on even $f$ and above $f=1000$. The count $d/t$ of vanishing components is
  given only for **prime** $f$: for composite $f$ imprimitive characters add zeros (at $f=21$ it
  predicts one and there are two), and `cuantos` declines there.
- `conrey` declines unless $q$ is a prime $\equiv3\pmod4$, $q>3$, the setting in which the note uses
  Conrey's identity. It computes $F_q$ exactly and does not sum the series. For $q\equiv3\pmod8$ it
  says the proposition applies and reports $R_q$; it does not test the positivity the proposition
  assumes.
- For $q\equiv3\pmod8$, `digitos` reports whether $R_q$ is empty and marks it as measured, not
  proved, as the note does.

## 🗂️ Files

| file | what it is |
|---|---|
| `segmento_II.py` | the tool — one file, standard library only |
| `II_certificados_A_bloques.csv` | the 104 block certificates, as written by `II_certificados_A.py`; read by `certificado` |
| `autotest_OUT.txt` | the archived run of `python segmento_II.py autotest` |
| `LICENSE` | CC0 1.0 Universal |

The scripts listed below this section are the gates of the note's *measured* statements and travel
with the note; the tool does not import any of them.

## 📜 Licence

**CC0 1.0 Universal** — public domain. Use it, cut it up, rewrite it or ship it inside something
else, without asking and without attribution. A citation is welcome and is not required.

## 📎 Citing

This tool travels with the note it implements, *Zeros of weighted partial sums of completely
multiplicative functions: Dirichlet characters, the Euler factor at 2 and Conrey's sine series*,
concept DOI [10.5281/zenodo.22847374](https://doi.org/10.5281/zenodo.22847374) — it always resolves to the current version.

## ✍️ Author

Carles Marín, independent researcher — `karlesmarin@gmail.com` —
[ORCID 0009-0007-5637-9688](https://orcid.org/0009-0007-5637-9688).
Claude (Anthropic) was used as a research assistant; every statement was checked by the author, who
takes full responsibility for the content.

---

# 🧮 Tools for *Zeros of weighted partial sums of completely multiplicative functions*

Scripts that generate or verify every number labelled *measured* in the note, and every figure, with
their archived output (`*_OUT.txt`). Python 3 (standard library; some use NumPy or Matplotlib);
the `.sage` files run in SageMath. Public domain (CC0).

**Every admissible length (Theorem A)**

| script | what it does |
|---|---|
| `II_certificados_A.py` | builds and checks, in integer arithmetic, the 104 block certificates for 5000 ≤ m < 10⁸; writes `II_certificados_A_bloques.csv` (~15 s) |
| `II_cola_analitica.py` | the inequalities of the proof for m ≥ 10⁸ (Dusart, Sárközy–Lev), with π(10⁸) exact |
| `II_ruta_A_exacta.py` | an explicit sign pattern with zero weighted sum for every admissible m in a range, sum recomputed (`python II_ruta_A_exacta.py 20000 5000`) |
| `II_segundo_momento.py` | the second-moment lemma, computed exactly from squarefree kernels |
| `II_lema_pares_directo.py` | the pairs lemma against the true subset sums of the first blocks |
| `II_lema_signos.py`, `II_ruta_A.py` | the uncovered width K(m) and the measured remainder (Figure: margin) |

**The diagonal and the rows**

| script | what it does |
|---|---|
| `diagonal_chi2.sage`, `primitividad_ceros.py` | the diagonal formula and the count d/t, with controls that must fail |
| `caracteres_orden_mayor.py`, `patron_cuadratico.sage`, `patron2.sage`, `patron3.sage` | vanishing components of higher order and the patterns off the diagonal |
| `filas_hipotesis.py`, `filas_sucesion.py`, `frontera_filas*.sage` | the rows: hypotheses, the counts a(m), the converse with a control |
| `II_exceso_mod4.py`, `II_mezcla.py` | the exact counts a(m) up to m = 120 (`II_exceso_mod4.py`, whose output `II_mezcla.py` reads) against the Gaussian and the mixture predictions |
| `II_momento_orden_r.py` | the second moment for values in μ_r and the constant C_4 |

**Along one function**

| script | what it does |
|---|---|
| `II_ceros_digitos.py` | the digit recursion at every m ≤ 10⁷, and the rules for zeros |
| `II_ceros_R.py`, `II_R_vacio.py` | the digit sets R_q for all primes up to 2000 / 10⁵ |
| `II_ceros_orden_r.py` | the zeros at (q^j − 1)/2 for odd characters of every order, exactly in ℤ[ζ_r] |
| `II_ceros_heegner.py` | the class-number-one chains, by direct summation |
| `II_bloques_ceros.py` | the block lemma at q = 59 and 347, R_103, and (q−5)/6 ∈ R_q |
| `II_s_menos_y_131.py` | S(q²m) = q²S(m) for s = −1, and the loop family at q = 131 |
| `II_automata_y_conrey.py` | the automaton of zeros (first zeros at 59 and 347; none at 179, 227 up to 40 digits) and the sine-series identity |
| `II_celdas_nulas_conrey.py` | the cells where the sine series vanishes identically, with controls |
| `SegmentoII.lean` | Lean 4 / Mathlib: the digit recursion, the propagation of zeros, and the pairs lemma; no `sorry`, standard axioms only (`lake env lean SegmentoII.lean` in a Mathlib project) |
| `figuras_II.py` | draws the five figures, in English or Spanish (`python figuras_II.py es`), from the archived outputs in this folder |
