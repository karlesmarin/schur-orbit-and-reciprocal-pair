# 🗃️ Ancillary files — *The ring of character values depends only on the segment: SU(n), Sp(2m) and GL(n) at a torsion element*

Everything the note calls *measured* is produced here, and the tool that does it is meant to be
used, not just archived.

## 🧰 The tool

**`segmento.py`** — one file, Python 3 standard library only. No Sage, no numpy, no install.

```
python segmento.py autotest          reproduce the numbers printed in the note
python segmento.py familias n q      which primes can divide the index          (Prop. support)
python segmento.py perfil   n q p    the layer profile W_0, W_1, ...            (Thm. local)
python segmento.py indice   n q p    v_p of the index, the slow honest path     (Thm. local (d))
python segmento.py puente   n q p    the same for GL(n), at half the degree     (Thm. bridge)
python segmento.py ceros    f        vanishing components on the diagonal       (Cor. how many)
python segmento.py fila     m        can row m be occupied at all               (Prop. rows)
```

As a library: `from segmento import familias, perfil, indice, puente, ceros, fila`.

`autotest` is the point. It recomputes, from nothing, the numbers the note prints: the profiles
`(1,11,21,22)` at `q=207` and `(1,15,29,30)` at `q=279`, the index `v_3([O:B_23(207)]) = 33`, the
shortcut against the slow path in five cells, the counts `d/t` at `f = 7, 23, 31, 47`, and which
rows can be occupied. If it prints anything but *all correct*, the note and the code disagree and
the note is the one to distrust.

**`medir_atajo.py`** — times the shortcut against the slow path and checks they agree. Nine cells,
the same integer in every cell, about half the wall clock at the largest. The times are a property
of the machine; `medir_atajo_OUT.txt` is the archived run.

## 🚧 What the tool refuses to do

Where the note marks a statement *measured* and not proved, the routine says so instead of handing
back a number it cannot justify:

- `puente` warns when `u <= d/2`, because `R = A` is proved only above that;
- `fila` distinguishes its two halves: `(2|n) = -1` empties the row and is **proved**; `(2|n) = +1`
  allows it and is **measured** on 39 rows and not proved.

That is the same honesty floor the note runs on, moved into the code.

## 🤝 Who might want this

The object is an order in a cyclotomic field and the question is its conductor — the everyday
computation in the algorithmics of complex multiplication, where the conductor of a CM order is
what one computes to identify an endomorphism ring or to place a curve in an isogeny graph. Two
routines replace a computation by a formula: `familias` decides maximality by divisibility alone,
and `puente` computes the non-real index from the real one, in half the degree.

## 📜 Licence

`segmento.py` and `medir_atajo.py` are released into the public domain under **CC0 1.0**: use,
cut up, rewrite or ship them without asking and without attribution. The note itself is not.

They also travel **outside** this note, as a software record of their own with their own DOI —
concept DOI [10.5281/zenodo.22834098](https://doi.org/10.5281/zenodo.22834098), mirrored at
[github.com/karlesmarin/schur-orbit-and-reciprocal-pair/tree/main/segmento](https://github.com/karlesmarin/schur-orbit-and-reciprocal-pair/tree/main/segmento)
— so that someone who never reads the note can find them, use them and cite them. The deposited
copy is byte for byte the one here; the publisher of that record checks it before depositing,
along with running `autotest` on the copy it is about to deposit.

## 🖼️ Sources and figures

The directory is flat: every script reads and writes in its own directory, and nothing outside it
is used.

| file | what it is |
|---|---|
| `segmento.tex`, `segmento_es.tex` | the LaTeX sources of the note, English and Spanish (`pdflatex`, three passes) |
| `fig_*.pdf`, `fig_*_es.pdf` | the figures, in both languages |
| `figuras.py` | draws every figure (`python figuras.py`; needs `numpy` and `matplotlib`) |
| `plano_denso.sage` (loads `inst_tipoA.sage`) | the data of the plane figure: `v_p` from the elementary divisors of the lattice, not from the formula |
| `plano_denso_OUT.txt` | its archived run, which `figuras.py` reads |
| `caracteres_orden_mayor.py` | vanishing components by characters of every order, in exact arithmetic (read by `figuras.py`) |
| `caracteres_orden_mayor_OUT.txt` | its archived run |

The `.sage` files need SageMath 10.x; the tool above does not. Every script here, and its archived
output, is CC0 like the tool; the sources and the figures are the note, under CC BY 4.0.

## 📏 Where each *measured* block comes from

One row per `\medido` block of `segmento.tex` (line numbers of the English source). The `.sage`
files require SageMath 10.x and load `inst_tipoA.sage` where needed; run them from this directory,
e.g. `docker run --rm -v "$PWD:/work" -w /work <sage image> sage anillos_iguales.sage`. The `.py`
files need Python 3 (`qlucas_signo.py` also needs `sympy`) and import only `segmento.py` from here.

| line | the block | script | archived output |
|---|---|---|---|
| 341 | level-rank: the value sets of SU(n) and SU(q-n) coincide, 25 of 25 pairs, up to 462 values | `complemento_levelrank.sage` (requires SageMath) | `complemento_levelrank_OUT.txt` |
| 363 | 287 cells for (a), as lattices; 90 complements and 82 periods, as lattices | not archived here | — |
| 549 | type-C law 18 of 18; type A 44 of 44, 15 of 18 outside, 3 fail with W_1=0 | `capa1_tipoA.sage` (requires SageMath) | `capa1_tipoA_OUT.txt` |
| 549 | statement (b): 193 of 193 stable, 98 of 127 unstable | `recuento_estable.sage` (requires SageMath) | `recuento_estable_OUT.txt` |
| 632 | 807 cells in the families, all stable | `familia_estable_r.py` | `familia_estable_r_OUT.txt` |
| 662 | 26 of 26 with v_p >= d-1; 23 of 23 with p dividing the index | not archived here | — |
| 694 | uI=I only for u=±1 (3.1e7 checks); the union-of-strata intervals at q'=6 | `inter_autocorrelacion.py` | `inter_autocorrelacion_OUT.txt` |
| 694 | the affine form (1.9e6 checks) | not archived here | — |
| 813 | 7714 cells outside the families, 80 stable, 26 maximal | not archived here | — |
| 813 | 195 cells, 39 with index > 1 | not archived here | — |
| 813 | Proposition *tangente*: 3060 of 3060 cells | not archived here | — |
| 813 | 147 926 unstable segments, none with all moments divisible by p | `soporte_por_momentos.py` | `soporte_por_momentos_OUT.txt` |
| 839 | r divides n+1: 543 cells, no defect | not archived here | — |
| 914 | Propositions *par2* (95 of 95) and *par1* (48 of 48) | `npar_y_soporte.sage` (requires SageMath) | `npar_y_soporte_OUT.txt` |
| 914 | 7580 comparisons, 288 cells, the odd columns losing rank | not archived here | — |
| 936 | the shadow rank drops to 1 exactly at n = 0,1 mod q (36 cells) | not archived here | — |
| 950 | twisted Gaussian binomials: 135 checks | `socratico_tipoA.sage`, part P3 (requires SageMath) | `socratico_tipoA_OUT.txt` |
| 976 | e_rs over 44 pairs: 154 coefficients, the unsigned formula failing in 30 | not archived here | — |
| 1007 | Theorem *sombra* 75 of 75; the description of W_1 81 of 81 | `sombra_capas.sage` (requires SageMath) | `sombra_capas_OUT.txt` |
| 1090 | 100 cells by the slow path | not archived here | — |
| 1250 | 78 of 78 lengths, 70 of 70 residues, 68 exact identities | not archived here | — |
| 1322 | the centred representatives at q=135 | not archived here | — |
| 1322 | the profile (1,0,1,2,2,2) at q=45 | not archived here | — |
| 1362 | r=41, p=11, q=451: profile (1,18,20,20), lengths 21 and 19 | not archived here | — |
| 1390 | the displayed identities, 306 of 306 | not archived here | — |
| 1390 | the index law, 13 of 13 by integer lattices | `puente_indice.sage` (requires SageMath; 12 cells, and the cell q=207 on its own) | `puente_indice_OUT.txt`, `puente_indice_207_OUT.txt` |
| 1390 | the seven steps of the proof, 66 of 66 | `puente_prueba.sage` (requires SageMath) | `puente_prueba_OUT.txt` |
| 1527 | the shortcut against the slow path, 9 cells | `medir_atajo.py` | `medir_atajo_OUT.txt` |

`datos_indice.json` is the table of lattice indices published with the previous note of the
series, byte for byte.
