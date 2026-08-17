# 🔢 Schur polynomials at a root-of-unity orbit and reciprocal pairs — scripts and data

Verification code, **saved output** and the paper itself for

> **Factorization and vanishing of Schur polynomials twisted by roots of unity and reciprocal pairs**
> Carles Marín
> [arXiv:2608.09619](https://arxiv.org/abs/2608.09619) · [doi:10.48550/arXiv.2608.09619](https://doi.org/10.48550/arXiv.2608.09619) · `math.CO` (cross-list `math.RT`) · CC BY 4.0

**This repository tracks v2.** The first version, announced 11 August 2026, was titled *Factorization
of Schur polynomials twisted by roots of unity and a reciprocal pair*: 39 pages, one reciprocal pair,
the factorization and its three consequences. The title changed when the vanishing results for an
arbitrary number of pairs were added — the zero locus at `t = 2` for every `r`, the odd-`t` criterion,
and the reflection of the excess part — which is most of what the paper grew by. The arXiv identifier
and the DOI are unchanged, and the old title stays on v1 in the version history.

The object is the Schur polynomial at a full root-of-unity orbit together with one free reciprocal
pair,

```
Phi_t(lambda; z) = s_lambda( 1, zeta, ..., zeta^(t-1), z, 1/z ),    zeta = exp(2 pi i / t),
```

a polynomial in `N = t + 2` variables, and — for the zero locus — its several-pair generalization

```
Psi_r(lambda) = s_lambda( 1, -1, z_1, 1/z_1, ..., z_r, 1/z_r ),     N = 2r + 2.
```

The published rank-one case at `t = 2` is Part IV of the series,
[doi:10.5281/zenodo.21463000](https://doi.org/10.5281/zenodo.21463000), whose own scripts live in
[`schur-nonidentity-o4`](https://github.com/karlesmarin/schur-nonidentity-o4). This repository is the
general-`t` sequel and is kept separate so that each paper has one artifact.

## 📦 What is here

The **scripts and their archived stdout**, and the paper in both languages with its sources and
figures. Every count the paper's verification table quotes can be located in [`outputs/`](outputs)
without running anything, and either edition rebuilds from what is in this directory —
`pdflatex orbit_pair.tex` needs nothing that is not here.

| | |
|---|---|
| **Group 1** — the evaluation and its sharpness | plain Python, needs only `mpmath` (plus `numpy`/`matplotlib` for figures) |
| **Group 2** — the zero locus for every `r` | Sage: exact multivariate Laurent rings and linear systems over `Q` |
| [`outputs/`](outputs) | full stdout of all **50** runs, 2026-07-30 to 2026-08-14 |
| [`orbit_pair.pdf`](orbit_pair.pdf) | **the paper**, 66 pp — the file arXiv carries as v2 — with its source `orbit_pair.tex` |
| [`orbit_pair_es.pdf`](orbit_pair_es.pdf) | **the Spanish edition**, 68 pp, with its source `orbit_pair_es.tex` |
| `orbit_pair_Z.pdf`, `orbit_pair_Z_es.pdf` | the same two with the **long abstract** |
| `fig_*.pdf` | the figures the two sources include — 19 each; the twentieth is drawn in the source |

```bash
pip install mpmath numpy matplotlib     # group 1
python theorem_full.py

sage selfcomp_law.sage                  # group 2
```

### 🇪🇸 Edición española

`orbit_pair_es.pdf` es el mismo artículo en castellano. arXiv sólo admite una lengua por envío, así
que la edición española vive aquí. No es un resumen: son las mismas secciones, los mismos enunciados
y los mismos números. `check_parity.py` compara los dos fuentes número a número (todo número de dos o
más cifras, con multiplicidad), entorno numerado a entorno numerado y etiqueta a etiqueta, y da **0
divergencias**; encontró dos defectos reales el 2026-07-30, que es por lo que existe. Los dos fuentes
están aquí, así que la comparación se corre tal cual:

```bash
python check_parity.py orbit_pair.tex orbit_pair_es.tex
```

Van dos PDF por lengua y **el cuerpo de los dos es el mismo**: salen del mismo fuente, y lo único que
los separa es un interruptor de compilación que elige entre el resumen corto —el que cabe en el campo
de metadatos de arXiv— y el largo, que explica el resultado con más sitio. Si vas a leer sólo uno, lee
el del resumen largo.

La edición inglesa es la de referencia: si alguna vez discreparan, manda la de arXiv.

## ✅ The load-bearing numbers

| what | script | result |
|---|---|---|
| **every displayed formula, against its definitions** | `AUDIT_FORMULAS.py` | **3297 evaluations over 9 formulas, 0 failures** |
| the evaluation, sign included | `theorem_full.py` | 959 exact, 476 zeros, **0 failures**; L4 724/724 |
| an independent second implementation | `law_control.py` | 959 / 476 / 0 along a different code path |
| Theorem 3.1 again, from the printed statement alone | `thm_main_independent.py` | **749 shapes, 0 failures**, sign included; decoy fails 35/133 |
| Proposition 3.10 the same way, both directions | `invariant_separates.py` | **676 shapes, 0 collisions**; decoy 104 |
| the sharpness controls (D3), (D4) | `falsify.py` | orbit 600/0, coset 383/**217**, free pair 200/**400** |
| the three arguments from the `t`-quotient | `d_from_quotient.py` | 2970/2970 |
| the short form of the sign (statement) | `sign_ayyer_idiom.py` | ordered **1496/1496**; 112 cells, 0 mixed |
| the short form of the sign (proof, step by step) | `sign_proof_check.py` | 1529 / 4000 / 1496 / 198, **0 failures** |
| the invariant is minimal, and the shift law | `invariant_minimal.py` | **75640 invariants, 0 collisions**; shift law 826/826, control fires at every even `t` |
| the sign of the extra family, both halves | `extra_sign.py` | 28 extras + 219 cores, **all +1**; `eps = -1` on 80 of 3038, so the claim is not vacuous |
| the concentric locus, in quotient coordinates | `concentric_locus.py` | 1331/1331; **0** concentric at odd `t`, 43 at even |
| the value is constant on each fibre of the invariant | `fig_fibres.py` | 860 + 883 partitions, worst spread **1e-37**; 121 invariants → **110** values at `t=4` |
| the second half of the involution, measured | `involution_runs.py` | max \|sigma\| = 1; **0 of 1406** runs break the alternation; 333/333 |
| the same, proved step by step | `alternation_proof.py` | 4823 / 4823 / 1341 / 217, **0 failures** |
| the vanishing criterion, both directions | `selfcomp_law.sage` | 9961 shapes, 318 vanishers, **0 false negatives, 0 false positives** |
| the converse at one pair | `close_X_r1.sage` | 14950 beta sets, **0 violations** |
| the converse inside Littlewood's range | `associates_witness.sage` | 372 witnesses, **0 failures** |
| the isolating-`mu` statement at `r = 2, 3` | `prove_W.sage` | 80 witnesses, **0 residue** |
| the core does not decide the vanishing | `core_vs_criterion.py` | 260 profiles, 0 profile→core clashes; 1338 shapes, **5 core classes carrying both behaviours** |
| everything in §8, rebuilt from the definitions | `AUDIT_ALL.sage` | **3414 checks, 0 failures** |

## 🧪 Controls that are meant to fail

A test that cannot fail proves nothing, so several of these scripts print a deliberate negative:

- **`falsify.py`** applies the same identity to the coset alphabet and to a *free* (non-reciprocal)
  pair. It fails in 217 and 400 of 600 cases — the reciprocity hypothesis is doing work.
- **`sign_ayyer_idiom.py`** prints, before the real test, the same comparison with the two blocks in
  arbitrary order: **592 agree, 904 disagree**. `lambda11` and `sgn(a1+a2-b1-b2)` are each
  antisymmetric under exchanging the blocks, so their product is invariant while the candidate
  formula is not; without a rule fixing which block is `A` the right-hand side is not a function of
  `lambda` at all.
- **`sign_proof_check.py`** checks the three steps of the proof rather than its conclusion, so a
  disagreement can be localised: the parity count on random words, where it can fail; the resulting
  formula for `inv(w) - inv(b_S)`; and the arithmetic fact that closes it. Its own control repeats
  the middle step with the ordering hypothesis dropped, and it fails 904 of 1496 — the hypothesis
  enters the proof at exactly one point, and this is that point.
- **`invariant_minimal.py`** grades the same population by a deliberately lossy invariant, the
  SUM `d_1+d_2+d_3` with the sign instead of the multiset: it collides **75284 times of 75640**.
  And it re-runs the shift law with the exponent `(-1)^(t+1)` dropped, which agrees at odd `t`,
  where the law is trivial, and fails on every shape at even `t`, which is where it has content.
- **`extra_sign.py`** would be satisfied for free by a sign routine that always returned `+1`, so
  it reports how many two-row shapes in the same range carry `eps = -1`: **80 of 3038**. It also
  refutes a near-miss closed form for the inversion count, `2m^2 + j`, on all 21 members with
  `j > 0`, and confirms that no shape outside the two families satisfies the equality (0 of 2835).
- **`concentric_locus.py`** claims that the concentric branch is empty at odd `t`, which a script
  can satisfy by finding nothing anywhere. So it prints the even-`t` count on the same range next to
  it: **43**. A zero beside a zero would prove nothing.
- **`selfcomp_law.sage`** runs an even-width control: self-complementary shapes of *even* width, none
  of which vanishes. The parity hypothesis is load-bearing rather than decorative.

## ⚠️ Two cautions worth reading before you believe a number

**`close_X_r1.sage` contains an exploratory sub-check that does not hold**, and it prints roughly
4180 mismatches out of 6655. That sub-claim ("the 4×4 determinant is the signed sum of single `S`'s")
appears nowhere in the paper, and the script now says so on the line above it. The result the paper
rests on is the case analysis printed afterwards, whose ten buckets sum to 14950 with 0 violations.

**`prove_W.sage` uses a different range per `r`** — `run(1,12)`, `run(2,10)`, `run(3,9)` — so its
counts of non-standard labels, 33 / 14 / 0, are **not comparable across `r`**. A non-standard label
needs `|nu| >= 2r+3`, so the usable slack above that threshold is 7 / 3 / 0: the zero at `r = 3` is
the sweep stopping at the threshold, not the obstruction disappearing. Reading a trend off those
three numbers is a mistake, and it is one this repository made before catching it.

## 🧾 Sums, stated

A few table entries are sums of printed per-`r` lines rather than a single printed total. They are:

```
9961 = 2157 + 4542 + 3262     shapes scanned
 318 =  242 +   55 +   21     vanishers
 372 =   76 +  144 +  152     associate witnesses
14950 = the ten parity buckets of close_X_r1
 184 =   94 +   90            criterion vs exact object
 236 =  139 +   97            type-D labels reduced (unstable_closed; NOT the 97+67
                              of fig_reduction, which counts label classes per panel)
  80 =   55 +   25            isolating witnesses
  21 =   18 +    3            residual labels, reduced exactly
```

## 🔄 Regenerating everything

```bash
python _audit_group1.py     # group 1 -> _out/
bash _audit_table.sh        # group 2, inside a Sage container
bash _save_outputs.sh       # consolidate into outputs/
```

Sage runs in Docker if you would rather not install it:

```bash
docker run --rm -v "$PWD:/work" -w /work sagemath/sagemath:latest sage selfcomp_law.sage
```

## 📐 Conventions

`t` is the order of the root of unity and `z` is the free reciprocal variable. Much of the
surrounding literature uses `t` for the free variable; these scripts do not.

In `theorem_full.setup` the two distinguished residue classes are ordered **by column**, which is
what the proof of the main theorem needs. The short form of the sign instead needs them ordered **by
residue**. `lambda11` is antisymmetric under exchanging the two blocks, so it must be recomputed
after any reordering — reusing it is a silent sign error.

Determinants use Gaussian elimination with partial pivoting rather than `mpmath.det`, because the
singular case is the one that has to be scored rather than raised.

## 📄 Licence and authorship

Carles Marín, with Claude (Anthropic) as an AI research assistant. Released for verification and
reuse; if the scripts are useful in your own work, a citation of the preprint is welcome:

```bibtex
@misc{marin2026orbitpair,
  author = {Carles Mar\'in},
  title  = {Factorization and vanishing of {S}chur polynomials twisted by roots of unity and reciprocal pairs},
  year   = {2026},
  eprint = {2608.09619},
  archivePrefix = {arXiv},
  primaryClass  = {math.CO},
  doi    = {10.48550/arXiv.2608.09619}
}
```
