# 🔢 Schur polynomials at a root-of-unity orbit and a reciprocal pair — scripts and data

Verification code and **saved output** for

> **Factorization of Schur polynomials twisted by roots of unity and a reciprocal pair**
> Carles Marín · preprint forthcoming

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

This repository carries the **scripts and their archived stdout**, not the paper. Every count the
paper's verification table quotes can be located in [`outputs/`](outputs) without running anything.

| | |
|---|---|
| **Group 1** — the evaluation and its sharpness | plain Python, needs only `mpmath` (plus `numpy`/`matplotlib` for figures) |
| **Group 2** — the zero locus for every `r` | Sage: exact multivariate Laurent rings and linear systems over `Q` |
| [`outputs/`](outputs) | full stdout of all 22 runs, 2026-07-30 |

```bash
pip install mpmath numpy matplotlib     # group 1
python theorem_full.py

sage selfcomp_law.sage                  # group 2
```

## ✅ The load-bearing numbers

| what | script | result |
|---|---|---|
| the evaluation, sign included | `theorem_full.py` | 959 exact, 476 zeros, **0 failures**; L4 724/724 |
| an independent second implementation | `law_control.py` | 959 / 476 / 0 along a different code path |
| the sharpness controls (D3), (D4) | `falsify.py` | orbit 600/0, coset 383/**217**, free pair 200/**400** |
| the three arguments from the `t`-quotient | `d_from_quotient.py` | 2970/2970 |
| the short form of the sign | `sign_ayyer_idiom.py` | ordered **1496/1496**; 112 cells, 0 mixed |
| the vanishing criterion, both directions | `selfcomp_law.sage` | 9961 shapes, 318 vanishers, **0 false negatives, 0 false positives** |
| the converse at one pair | `close_X_r1.sage` | 14950 beta sets, **0 violations** |
| the converse inside Littlewood's range | `associates_witness.sage` | 372 witnesses, **0 failures** |
| the isolating-`mu` statement at `r = 2, 3` | `prove_W.sage` | 80 witnesses, **0 residue** |
| everything in §8, rebuilt from the definitions | `AUDIT_ALL.sage` | **3414 checks, 0 failures** |

## 🧪 Controls that are meant to fail

A test that cannot fail proves nothing, so three of these scripts print a deliberate negative:

- **`falsify.py`** applies the same identity to the coset alphabet and to a *free* (non-reciprocal)
  pair. It fails in 217 and 400 of 600 cases — the reciprocity hypothesis is doing work.
- **`sign_ayyer_idiom.py`** prints, before the real test, the same comparison with the two blocks in
  arbitrary order: **592 agree, 904 disagree**. `lambda11` and `sgn(a1+a2-b1-b2)` are each
  antisymmetric under exchanging the blocks, so their product is invariant while the candidate
  formula is not; without a rule fixing which block is `A` the right-hand side is not a function of
  `lambda` at all.
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
 206 =  139 +   67            type-D labels
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
reuse; if the scripts are useful in your own work, a citation of the preprint is welcome.
