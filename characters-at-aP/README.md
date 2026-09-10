# Integer character values on the rho-line of Sp(2m) — code and data

Companion material for *Integer character values on the $\rho$-line of $\mathrm{Sp}(2m)$*
(Carles Marín), deposited on Zenodo:
[doi:10.5281/zenodo.22683560](https://doi.org/10.5281/zenodo.22683560). That is the **concept
DOI** — it always resolves to the current version — and it carries the paper in both editions,
English (33 pp) and Spanish (33 pp), under CC BY 4.0.

The PDFs are not duplicated here: Zenodo carries them, and one copy is enough. What is here is the
code, the archived output, and the two LaTeX sources with their figures, so that either edition
rebuilds from this directory alone.

This is **not** a third part of the series in the root of this repository. Those two papers are
about Schur polynomials at a root-of-unity orbit; this one is about characters of $\mathrm{Sp}(2m)$
at the rational points of Kostant's $\rho$-line. It lives here, rather than in a repository of its
own, so that the material of one line of work is found in one place.

## What the paper proves

Kostant's element $a_P = \rho(\zeta_t)$, $t = 2m+2$, is one point of the line
$\theta \mapsto \exp(2\pi i\theta\rho^{\sharp})$; the alcove argument that makes every character
$0$ or $\pm 1$ there needs $\rho$ to be regular for the relevant arrangement. At the **powers** of
$a_P$ it is not. Three things are proved:

1. **A classification.** Of all rational points $g_{1/q}$ of the line, exactly four families make
   *every* character of $\mathrm{Sp}(2m)$ an integer: $q \mid 2m$, $q \mid 2m+1$, $q \mid 2m+2$,
   and a sporadic $q = 6$ with $m \equiv 1 \pmod 3$. The powers of $a_P$ are the family
   $q \mid 2m+2$ — one case of a list, not the boundary of what can be said.
2. **A vanishing criterion.** With $\ell = \lambda+\rho$ and
   $N_q(x) = \#\{\alpha > 0 : q \mid (x,\alpha)\}$, the character is non-zero exactly when
   $N_q(\ell) = N_q(\rho)$. The count can never come out below; equality is the whole condition.
3. **The value, in type C, for every power and with its sign** —
   $\chi_\lambda(a_P^{\,k}) = (-1)^{E_{p,q}(\lambda)} Z_q(\lambda+\rho)/Z_q(\rho)$ on the surviving
   weights and $0$ elsewhere — whose modulus **factors as a product of dimensions of irreducible
   representations** of $GL_n$, $Sp_{2n}$ and $SO_{2n+1}$, one per folded class. Two closed counts
   of the surviving weights follow for every $q$, and $G_2$ is settled in closed form at $q = 3$.

## What is not new, stated plainly

At $a_P$ **itself, nothing here is new.** The criterion and its sign are Corollary 3.2 of
Cellini–Möseneder Frajria–Papi (2007), in every type and on this very line; read on coroots the
count is Corollary 3.7 of Nadimpalli–Pattanayak–Prasad; the principal specialisation the argument
runs on is classical; and at $q = 2$ the factorisation into dimensions is Theorem 5.1 of Karmakar
(2025). The paper carries a table saying line by line which statement is whose. What is done here
is the powers, where $\rho$ is no longer regular — and, in type C, the exact value with its sign,
the factorisation into dimensions at $q \ge 3$, and the closed counts.

## What is here

101 files, flat, named after what they do.

| | |
|---|---|
| `*.py` with a `*_OUT.txt` beside it | the computation behind a framed formula or a quoted count, and the run that produced it |
| `check_*.py` | the ten verification gates that run over the manuscript itself — figures, formulas, numbering, cited attributions, bibliography, EN/ES parity, and the overlap with the author's own [arXiv:2609.02630](https://arxiv.org/abs/2609.02630) |
| `characters_at_aP.tex`, `characters_at_aP_es.tex` | the two editions; `pdflatex` needs nothing that is not in this directory |
| `fig_*.pdf`, `fig_*.py`, `figs.py` | the six figures and the scripts that draw them |
| `root_data.json`, `biblio_registro.json` | the root data and the bibliographic register the gates read |

Forty-four of the forty-seven scripts carry their archived run. The three that do not are
`fig_value.py` and `fig_walls.py`, whose output is the figure itself and it travels, and
`verifica_informe.py`, which recomputes four statements from scratch and prints them.

## How to read a claim

Find the script named after it, read its header — each says what it tests and what would falsify
it — and either read the `*_OUT.txt` beside it or run it again. The gates are built with a
**control that must fail**: a gate whose decoy also passes is reported as untested, not as
confirmed.

Everything is plain Python 3 on the standard library, except the three figure scripts, which need
`matplotlib`, and the three gates that measure the rendered page with `pymupdf` — the text layer of
a PDF is not the page, and those three look at the page. A script `foo.py` writes `foo_OUT.txt`.

## What is deliberately not here

An eleventh gate exists and does not travel: its subject is a letter, not a statement of the paper.
It is held out **by location, not by a list of names** — a name list cannot see what is inside a
file, and excluding one by name leaves the references to it dangling in what does ship. A check
refuses to release a bundle in which anything imports or reads a file that is not in it; this
bundle passes it.

## Licence and citation

Cite the paper — concept DOI [10.5281/zenodo.22683560](https://doi.org/10.5281/zenodo.22683560).
The scripts may be used freely; they are research code, written to be read and doubted rather than
to be depended on.
