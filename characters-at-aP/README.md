# 🔢 Integer character values on the rho-line of Sp(2m) — code and data

Companion material for *Integer character values on the $\rho$-line of $\mathrm{Sp}(2m)$*
(Carles Marín Muñoz), **version 2**, deposited on Zenodo. Concept DOI
[10.5281/zenodo.22683560](https://doi.org/10.5281/zenodo.22683560) — it always resolves to the
current version; version 2 is [10.5281/zenodo.24073781](https://doi.org/10.5281/zenodo.24073781).
The record carries the paper in both editions, English (36 pp) and Spanish (37 pp), and an
ancillary archive of 107 files.

🗺️ **Blueprint** — [read the paper statement by statement](https://karlesmarin.github.io/blueprints/characters-at-aP/):
every result with its number, what it uses and where it is used, what was checked besides the
written proof and how far, the dependency graph and a one-page map, two labs in the browser (the
four families of denominators of Theorem 5.1, and the character at a power of $a_P$ for any $m$, $k$
and $\lambda$), and every formalized statement linked to its Lean declarations.

Both PDFs are here too — `characters_at_aP.pdf` and `characters_at_aP_es.pdf`, byte for byte the
files deposited on Zenodo — together with the 107 files of the ancillary archive: the code, the
exact certificates, the gates, the archived output, and the two LaTeX sources with their figures,
so that either edition also rebuilds from this directory alone.

The `blueprint/` folder beside this README is an earlier site, kept so that old links still
resolve; it predates version 2 and is numbered differently. The current site is the one above.

This is **not** a third part of the series in the root of this repository. Those two papers are
about Schur polynomials at a root-of-unity orbit; this one is about characters of $\mathrm{Sp}(2m)$
at the rational points of Kostant's $\rho$-line. It lives here, rather than in a repository of its
own, so that the material of one line of work is found in one place.

## 🎯 What the paper proves

Kostant's element $a_P = \rho(\zeta_t)$, $t = 2m+2$, is one point of the line
$\theta \mapsto \exp(2\pi i\theta\rho^{\sharp})$; the alcove argument that makes every character
$0$ or $\pm 1$ there needs $\rho$ to be regular for the relevant arrangement. At the **powers** of
$a_P$ it is not. Three things are proved:

1. **A classification (Theorem 5.1).** Of all rational points $g_{1/q}$ of the line, exactly four
   families make *every* character of $\mathrm{Sp}(2m)$ an integer: $q \mid 2m$, $q \mid 2m+1$,
   $q \mid 2m+2$, and a sporadic $q = 6$ with $m \equiv 1 \pmod 3$. Integrality forces the
   characteristic polynomial of the natural representation into $\mathbb{Z}[X]$, so the spectrum
   must be Galois-stable; classifying the stable spectra is elementary. The necessity also follows
   from Gannon–Schopieray (2019). The powers of $a_P$ are the family $q \mid 2m+2$ — one case of a
   list, not the boundary of what can be said.
2. **A vanishing criterion (Theorem 3.1).** With $\ell = \lambda+\rho$ and
   $N_q(x) = \#\{\alpha > 0 : q \mid (x,\alpha)\}$, the character is non-zero exactly when
   $N_q(\ell) = N_q(\rho)$. The count can never come out below; equality is the whole condition.
3. **The value, in type C, for every power and with its sign (Theorem 7.1)** —
   $\chi_\lambda(a_P^{\,k}) = (-1)^{E_{p,q}(\lambda)} Z_q(\lambda+\rho)/Z_q(\rho)$ on the surviving
   weights and $0$ elsewhere. That modulus **factors as $\kappa(\lambda)\in\{1,2\}$ times a product
   of dimensions of irreducible representations** of $GL_n$, $Sp_{2n}$ and $SO_{2n+1}$, one per
   folded class (Theorem 10.1). A closed count of the surviving weights follows for every $q$
   (Theorem 11.1), and $G_2$ is settled in closed form at $q = 3$.

## 📚 What is not new, stated plainly

At $a_P$ **itself, nothing here is new.** The criterion and its sign are Corollary 3.2 of
Cellini–Möseneder Frajria–Papi (2007), in every type and on this very line; read on coroots the
count is Corollary 3.7 of Nadimpalli–Pattanayak–Prasad; the principal specialisation the argument
runs on is classical; and at $q = 2$ the factorisation into dimensions has the two-block form of
Theorem 5.1 of Karmakar (2025). The paper carries a table saying line by line which statement is
whose. What is done here is the powers, where $\rho$ is no longer regular — and, in type C, the
exact value with its sign, the factorisation into dimensions at $q \ge 3$, and the closed count.

## 🔒 What is formalized, and what is not

In **Lean 4 with Mathlib** — [karlesmarin/formalizations](https://github.com/karlesmarin/formalizations)
at commit [`6283de8`](https://github.com/karlesmarin/formalizations/commit/6283de8ed4582c42d045e8b9e71b1a8512340c77),
no `sorry`, standard axioms only — are the type-$C_m$ statements on the root products $F_x$ and
their quotient $Q$ (including that $F_\rho$ divides $F_\ell$ over any commutative ring), the
spectral classification behind the first theorem, and, stated for $Q(\zeta)$, the criterion, the
value with its sign and the factorisation into dimensions in both parities.
**Outside Lean:** the identification of $Q(\zeta)$ with the character (Weyl's formula), the
equivalence between spectral and character integrality, that the $\nu^{(c)}$ are partitions, and
the step from a dominant $\lambda$ to $\ell$. The paper's section *Formalization* gives the scope
statement by statement. The exact certificates below are checks on finite ranges, not proofs.

## 📦 What is here

The two PDFs and the 107 files of the ancillary archive, flat, named after what they do — the same
files as `characters_at_aP_ancillary.zip` on Zenodo.

| | |
|---|---|
| `characters_at_aP.pdf`, `characters_at_aP_es.pdf` | the paper, 36 pp in English and 37 pp in Spanish — the same files Zenodo carries |
| `*.py` with a `*_OUT.txt` beside it | the computation behind a framed formula or a quoted count, and the run that produced it |
| `certificar_exacto.py` | the exact certificates C1–C4, which replace the floating-point checks of version 1, each with a control built to fail |
| `check_*.py` | the eleven verification gates that run over the manuscript itself — figures, formulas, numbers, numbering, coverage, cited attributions, bibliography, EN/ES parity, the paragraph on Gannon–Schopieray, the restructuring, and the overlap with the author's own [arXiv:2609.02630](https://arxiv.org/abs/2609.02630) |
| `characters_at_aP.tex`, `characters_at_aP_es.tex` | the two editions; `pdflatex` needs nothing that is not in this directory |
| `fig_*.pdf`, `fig_*.py`, `figs.py` | the six figures and the scripts that draw them |
| `root_data.json`, `biblio_registro.json` | the root data and the bibliographic register the gates read |

Forty-seven of the fifty scripts carry their archived run. The three that do not are
`fig_value.py` and `fig_walls.py`, whose output is the figure itself and it travels, and
`verifica_informe.py`, which recomputes four statements from scratch and prints them.

## 🧪 How to read a claim

Find the script named after it, read its header — each says what it tests and what would falsify
it — and either read the `*_OUT.txt` beside it or run it again. The gates are built with a
**control that must fail**: a gate whose decoy also passes is reported as untested, not as
confirmed.

Everything is Python 3. Most scripts use the standard library only; six need `mpmath`, two (the
certificates and `comprobar_v253.py`) need `sympy`, the three figure scripts need `matplotlib`, and
the three gates that measure the rendered page need `pymupdf` — the text layer of a PDF is not the
page, and those three look at the page. A script `foo.py` writes `foo_OUT.txt`.

## 🚫 What is deliberately not here

A twelfth gate exists and does not travel: its subject is a letter, not a statement of the paper.
It is held out **by location, not by a list of names** — a name list cannot see what is inside a
file, and excluding one by name leaves the references to it dangling in what does ship. A check
refuses to release a bundle in which anything imports or reads a file that is not in it; this
bundle passes it.

## 📄 Licence and citation

Cite the paper — concept DOI [10.5281/zenodo.22683560](https://doi.org/10.5281/zenodo.22683560)
(version 2: [10.5281/zenodo.24073781](https://doi.org/10.5281/zenodo.24073781)). The Zenodo record
is released under the Apache License 2.0. The scripts are research code, written to be read and
doubted rather than to be depended on.

This work was carried out with an AI assistant; the author checked every statement and takes full
responsibility for the content.
