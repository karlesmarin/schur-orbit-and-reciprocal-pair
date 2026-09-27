# Ancillary files: the conductor note

Scripts and archived outputs behind every number labelled *measured* or *observed*, and every figure,
in the note on the conductor of the orders `A_m(q)`. Statement, figure and section numbers refer to the
English edition. The directory is flat: every script reads and writes in its own directory, and nothing
outside it is used.

Code comments are in Spanish or English.

## Requirements

- Python 3 with `numpy`, `sympy`, `matplotlib` (checked with Python 3.14, numpy 2.4, sympy 1.14,
  matplotlib 3.10).
- SageMath 10.x for the `.sage` files (the outputs were produced with SageMath 10.9 in Docker).
  Every `.sage` file is run from this directory, e.g.
  `docker run --rm -v "$PWD:/work" -w /work <sage image> sage item3.sage`.
  The `.sh` drivers are bash and are meant to run inside that container, in this directory.

## Where each number comes from

| in the note | script | output |
|---|---|---|
| orders, conductors, Figure 3 | `orden_generado.sage`, `factorizar_conductor.sage`, `h2_residuo.sage`; gathered by `datos_indice.py` | `celdas_m*_OUT.txt`, `h2_OUT.txt`, `r2_OUT.txt`, `datos_indice.json` |
| §9 *Orders*: closure = Sage order (22 cases) | `orden_generado_control.sage` | `orden_generado_control_OUT.txt` |
| §9 *Layers*: 184 pairs; Theorem 5 | `perfil_W.py`, `ley_conductor.sage` | `perfil_W_OUT.txt`, `ley_conductor_ancho_OUT.txt`, `ley_conductor_repro_OUT.txt` |
| §7 *Observed and certified*: 1183 pairs, 215 observed | `catalogo_e.py` | `catalogo_e_OUT.txt` |
| §7: the 215 certified by Proposition 6(b); §9: Proposition 8 and the symmetry of Proposition 10 on the 168 lattice pairs inside the hypotheses of §5 | `capas_locales.py` | `capas_locales_OUT.txt` |
| §7: profiles against `v_p(N)` | `catalogo_N.py` | `catalogo_N_OUT.txt` |
| Figure 6, Theorem 12 | `bernoulli_W1.py`, `bernoulli_2m1.py` | `bernoulli_W1_OUT.txt`, `bernoulli_2m1_OUT.txt` |
| §9: exact `h^-` for the 22 primes of Figure 6 | `h_menos_exacto.py` | `h_menos_exacto_OUT.txt` |
| Remark 15 (984 pairs) | `rango_compuesto.py` | `rango_compuesto_OUT.txt` |
| §9: `Psi_q = Psi_{q'}^E (mod p)`; ranks by orbits | `psi_y_orbitas.py` | `psi_y_orbitas_OUT.txt` |
| Proposition 13: 1505 pairs; exact `h^-` for the 66 values `n < 90`, `n != 2 mod 4` (used by Figure 7, `item3*.sage`) | `stickelberger_smith.sage` | `stickelberger_smith_OUT.txt` |
| §9: `h^-` against PARI, with proof for `phi(n) <= 20`, under GRH for `20 < phi(n) <= 32` | `h_menos_pari.sage`, `h_menos_pari_grh.sage` | `h_menos_pari_OUT.txt`, `h_menos_pari_grh_OUT.txt` |
| §6.1 *A second example*: lattice computation of `A_22(69)` (index, conductor, layers, type two ways) | `orden_A22_69.sage` | `orden_A22_69_OUT.txt` |
| §7 *What the figure already certifies*: Corollary 7 on the rows of Figure 8 | `item1.py` | `item1_OUT.txt` |
| §7: lattices `A_19(117)` and `A_22(225)`; §6.2 *A larger delta*: lattice `A_15(155)` | `lat.sage` (loads `capas.sage`) | `lat_A19_117_OUT.txt`, `lat_A22_225_OUT.txt`, `lat_A15_155_OUT.txt` |
| Proposition 19; Proposition 20, the equality `W_P = W_1(A_0)` (118 cases) | `item2.sage` (loads `capas.sage`) | `item2_OUT.txt` |
| §7: the finite logarithm (16 cases, plus controls that must fail) | `item2_log.py` | `item2_log_OUT.txt` |
| Corollary 17, factor at 2, prime `q'` (438 cases) | `item3.sage` | `item3_OUT.txt` |
| near-miss control for Corollary 17: pairs with `p | h^-`, outside the hypothesis | `item3_ctrl.sage` | `item3_ctrl_OUT.txt` |
| Proposition 11: the 63 lattice pairs with `e = 2` | `item4.sage` (loads `lat.sage`; reads `gor_sweep_OUT.txt`) | `item4_OUT.txt` |
| Proposition 11: the 403 pairs with `q <= 250` and certified `e = 2`, 19 with `kappa > 1`, 9 also by a lattice colon; §6.1: `A_40(451)` at 11 | `item4b.sage` (loads `lat.sage`) | `item4b_OUT.txt`, `item4b_hs_OUT.txt` |
| Proposition 20: descent, ring of square zero (28 + 4 cases) | `check_descenso.py` (`python check_descenso.py p,a,n,M ...`) | `out_batch.txt`, `out_hminus.txt` |
| §7, after Proposition 20: the non-centred segment, where `pi^P` does enter (4 cases) | `control.py` (imports `check_descenso.py`; `python control.py 3,1,5,1 3,2,5,1 3,1,7,1 5,1,7,1`) | `control_OUT.txt` |
| Propositions 19 and 20: ladder and descent, 953 cases (51 skipped for cost) | `descent.sage` (loads `inst.sage`; `sage descent.sage SHARD NSHARDS COSTCAP`, 8 shards concatenated) | `descent_OUT_all.txt` |
| Lemma 14(ii): column identities (odd `r <= 45`) | `duplicacion.py` | `duplicacion_OUT.txt` |
| Lemma 14(ii): lattices, elementary divisors, `Delta` (22 values of `r`); Theorem 12(a) for `q' = 2r` end to end (304 cases) | `l2r.sage` (loads `inst.sage`; reads `stickelberger_smith_OUT.txt`) | `l2r_OUT.txt` |
| Lemma 14(i): shift identity, `L' = L/2`, elementary divisors and ranks (34 even `n <= 70`); Theorem 12(a) on the half-period families with `4 | q'` end to end (158 lines, 156 distinct cases: the two with q = 156, p = 3 are printed twice), and 12 controls with `p | B` | `medio.sage` (loads `inst.sage`; reads `stickelberger_smith_OUT.txt`) | `medio_OUT.txt` |
| §6.1, *One identity, two regimes*: `2 hat_s_n = s_n(2c) - s_n(c)` for `n <= 100`; the first moment as count times midpoint on the four stable families (98 cases); the rank `d - delta_l` for general `l` (52 pairs, over `Q` and mod `p <= 13`); the second-order expansion of `eta_j` (130 pairs) | `cruce.sage` (reads `stickelberger_smith_OUT.txt`) | `cruce_OUT.txt` |
| Corollary 18 (the index): `v_p([O:A]) = 2d-1-W_1` against the lattice indices of `datos_indice.json` (104 pairs; 103 in explicit form; 3 pairs fall outside, all with `p | N`) | `indice.sage` (loads `inst.sage`; reads `datos_indice.json`, `stickelberger_smith_OUT.txt`) | `indice_OUT.txt` |
| §7: the second layer as the square of the first (238 pairs with `p` not dividing `N`); and, beyond the descent, `W_{P+1}` against the module generated by `c -> hat(c)^2` and `W_{P+2} = W_P` (38 pairs; it holds in 26 of them, all 8 with `p = 3` included, counted by `capa2_recuento.py`) | `capa2.sage` (loads `inst.sage`), `capa2_recuento.py` | `capa2_OUT.txt`, `capa2_recuento_OUT.txt` |
| §7: the hypothesis of Corollary 18 is not superfluous — 359 stable subspaces generated by random vectors (all those with `dim V > d/2` square to the whole ring, 7 of the 9 with `dim V <= d/2` do not), the constants, and the ladder layer as a stable line | `cuadrado.sage` (loads `inst.sage`) | `cuadrado_OUT.txt` |
| Corollary 17 for odd `q'`, prime and composite (1038 cases, type in 18; with the 10 of dos_tipo4_OUT.txt, two in common, 26 distinct) | `dos.sage` (loads `capas.sage`, `inst.sage`) | `dos_OUT.txt` |
| Corollary 17: near-miss controls with `p | h^-` (7 pairs, 3 differ) | `dos_ctrl.sage` | `dos_ctrl_OUT.txt` |
| Corollary 17: the type for `delta <= 4` (10 cases) | `dos_tipo4.sage` | `dos_tipo4_OUT.txt` |
| §6.2 *Composite q'*: lattices `A_7(105)`, `A_31(315)` (and `A_25(357)` without the colon type) | `latchk.sage` (loads `lat.sage`; `sage latchk.sage m q p`) | `lat_A7_105_OUT.txt`, `lat_A31_315_OUT.txt`, `lat_A25_357_OUT.txt` |
| the layer instrument `inst.sage` against `capas.sage` and against `perfil_W.py` | `xval_inst.sage`, `xval_inst_dimsB.py` | `xval_inst_OUT.txt`, `xval_inst_dimsB_OUT.txt` |
| numbers the other outputs did not print (added in version 2): the 21 `q'` whose square unit matrix is singular over `Q` (caption of `fig_stickel`); the split of the 168 profiles into first-layer bands (`fig_zonas` and its paragraph); the whole of Sinnott's formula on the 66 `n` | `cifras_sin_salida.py` (reads `perfil_W_OUT.txt`, `stickelberger_smith_OUT.txt`) | `cifras_sin_salida_OUT.txt` |
| the attribution to Kostant (version 2): with the form normalised by `(theta,theta) = 1/h^v`, his element `exp(2 pi i 2 rho)` is `g_{1/(2m+2)}` in `Sp(2m)`, `exp(2 pi i rho^v/(2m+2))` is not, and every character takes a value in `{0, +-1}` there (`m <= 5`, 701711 dominant weights) | `kostant_check.py` | `kostant_check_OUT.txt` |
| all figures | `figuras.py` | `fig_*.pdf` |
| every claim made by a figure caption, against the data the figure draws (7 checks) | `auditar_figuras.py` | `auditar_figuras_OUT.txt` |

### Details of the SageMath rows

| file | what it does |
|---|---|
| `orden_generado.sage` | `og_orden(L, gens)`: the order `Z[g_1..g_m]` as a lattice, by closure (HNF of the lattice plus its products with the generators, to a fixed point). Loaded by `factorizar_conductor.sage`, `orden_generado_control.sage`, `orden_A22_69.sage`. |
| `orden_generado_control.sage` | On 22 small cases `og_orden` gives the same lattice (same HNF, same index) as Sage's `L.order`; a decoy with one generator removed must change the ring. |
| `factorizar_conductor.sage` | For each cell `(m,q)` (env `CASILLAS="m:q,..."`): the index `[O+:A]`, the conductor as a lattice, and its factorisation, maximising only at the primes of its norm. One table row per cell. |
| `ley_conductor.sage` | Earlier sweep over `2 <= m <= 5` printing the same table row (index, conductor norm, `[A:c]`, Gorenstein, primes `(p,f,e)`), plus its own controls. Env `GRADO_MAX`, `Q_MAX`, `M_MIN`, `M_MAX`. |
| `h2_residuo.sage` | For each cell and each prime `p` of the conductor: `a = v_p([O+:A])`, `b = v_p([A:c])`, the primes above `p` with `(f,e)`, and the image of `A` in the residue fields. One line per pair: `H2|m|q|p|a|b|dim|forma|...`. |
| `capas.sage` | Library: the layers `W_i`, `kappa` and the type modulo `p` of `A_m(q)`, by linear algebra in `F_p[X]/(Psi_{q'}^K)`, `K <= E`. An implementation independent of `perfil_W.py`. |
| `lat.sage` | Library and script: `A_m(q)` as an integer lattice in `O = Z[alpha]`; index, exponent `e` at `p`, layers by lattice intersections, type by the colon lattice; `kappa` and the type modulo `p` by `capas.sage`. `sage lat.sage m q p colon`. |
| `inst.sage` | Library: the layers of `A_m(q)` in `Z[xi]/(p, pi^K) = F_p[y]/(Phi_{q'}(y)^K)`, `K <= E`, with the `e_k` built by a product tree (FLINT); much faster than `capas.sage` for large `m`, and checked against it and against `perfil_W.py`. Also `moment_rank`, the rank of `[S(u^-1 c)]` mod `p`. |
| `gor.sage` | Independent Gorenstein test on lattices (type by the colon `(A : m)`), sweep `sage gor.sage list 100 20` (`q <= 100`, degree `<= 20`); produced `gor_sweep_OUT.txt`, which `item4.sage` reads. |
| `stickelberger_smith.sage` | Elementary divisors of the moment matrix `M_n` against the exact `h^-(Q(zeta_n))` from the analytic class number formula in cyclotomic arithmetic; ranks mod `p`; control against tabulated values. |
| `h_menos_pari.sage`, `h_menos_pari_grh.sage` | `h(Q(zeta_n)) / h(Q(zeta_n)^+)` from PARI's class group (`proof=True`, resp. `proof=False`), against the analytic `h^-`. |

Runs of the longer SageMath files:

- `item4b_OUT.txt`: `sage item4b.sage 250 40` (lattice colon check for degree `<= 40`; about two minutes).
  It was archived before the final `HS12.2.3 on e=2 cases` tally line was added to the script, so that
  line is missing there; the rest is reproduced exactly.
- `item4b_hs_OUT.txt`: `sage item4b.sage 250 0` (no lattice check), with the tally line; the same 403 cases.
- `lat_A15_155_OUT.txt`: `sage lat.sage 15 155 5 colon` (about 6 minutes); `lat_A19_117_OUT.txt`:
  `sage lat.sage 19 117 3 colon`; `lat_A22_225_OUT.txt`: `sage lat.sage 22 225 5 colon` (about 7 minutes).

### Lattice sweeps: drivers and logs

Drivers (one Sage process per cell, resumable; a log records every session with its script and cell list):

| file | what it does |
|---|---|
| `correr_por_casilla.sh` | First sessions of `celdas_m*_OUT.txt`: `sage < factorizar_conductor.sage`, one cell at a time. |
| `correr_reanudable.sh` + `casillas_hechas.py` | Resumable runner: skips cells already done in the log, appends each cell only when it finishes. |
| `cola_16sep.sh` | Second sessions of `celdas_m6/7/8_OUT.txt` (retries of cells that had run out of memory). |
| `h2_barrido.sh` + `h2_lista.txt` | `h2_OUT.txt`: `h2_residuo.sage` over the cell list in `h2_lista.txt`. |
| `r2_barrido.sh` | `r2_OUT.txt`: `h2_residuo.sage` over 15 further cells. |

| log | produced by |
|---|---|
| `celdas_m2_OUT.txt`, `celdas_m3_OUT.txt` | `correr_por_casilla.sh` + `factorizar_conductor.sage` |
| `celdas_m6_OUT.txt`, `celdas_m7_OUT.txt`, `celdas_m8_OUT.txt` | the same, then `cola_16sep.sh` (`correr_reanudable.sh` + `factorizar_conductor.sage`) |
| `h2_OUT.txt` | `h2_barrido.sh` (`h2_residuo.sage`) |
| `r2_OUT.txt` | `r2_barrido.sh` (`h2_residuo.sage`) |
| `ley_conductor_repro_OUT.txt` | `ley_conductor.sage`, `m <= 5`, `q <= 70`, degree `<= 12` (defaults) |
| `ley_conductor_ancho_OUT.txt` | `ley_conductor.sage`, `m <= 5`, `q <= 120`, degree `<= 20` |

The first sessions of the `celdas_m*` logs were run with a version of `factorizar_conductor.sage` that
built `A` with Sage's `L.order` instead of `og_orden`; the two constructions give the same lattice
(`orden_generado_control_OUT.txt`).

`perfil_W.py` reads exactly the logs named in `reparto_p2.LOGS_RETICULO` (`celdas_m6/7/8`,
`ley_conductor_ancho`, `ley_conductor_repro`) and `LOGS_H2` (`h2_OUT.txt`, `r2_OUT.txt`); removing any
one of them lowers the count. `datos_indice.py` reads `celdas_m*_OUT.txt`, `h2_OUT.txt`, `r2_OUT.txt`.
`capas_locales.py` and `figuras.py` (Figure 4) read `perfil_W_OUT.txt`; `figuras.py` (Figure 7),
`item3.sage` and `item3_ctrl.sage` read `stickelberger_smith_OUT.txt`.

Modules imported by the Python scripts: `reparto.py` (reads table rows from the logs), `reparto_p2.py`
(layer arithmetic, and the explicit list of logs), `perfil_W.py` (`dims_B`, the layers modulo `p`;
imported by `capas_locales.py`, `item1.py`, `figuras.py`), `mc2_combinatorio.py`, `bernoulli_W1.py`
(imported by `figuras.py`).

`figuras.py` writes `fig_estructura` (Figure 1), `fig_historia` (2), `fig_paisaje2d` (3),
`fig_simetria` (4), `fig_sierra` (5), `fig_clase` (6), `fig_stickel` (7), `fig_perfil` (8), each also
with an `_es` variant for the Spanish edition; `fig_paisaje` is a 3D view not used in the note.

## How to run

Python, from this directory:

```
python datos_indice.py
python perfil_W.py
python catalogo_e.py
python capas_locales.py
python catalogo_N.py
python bernoulli_W1.py
python bernoulli_2m1.py
python h_menos_exacto.py
python rango_compuesto.py
python psi_y_orbitas.py
python item1.py
python item2_log.py
python figuras.py [png_preview_dir]
```

Each prints to standard output; compare with the corresponding `*_OUT.txt` (lines with timings differ).
`datos_indice.py` rewrites `datos_indice.json`, and `figuras.py` writes the `fig_*.pdf` files here.
`capas_locales.py` takes about a minute.

SageMath, from this directory (each under two minutes):

```
sage orden_generado_control.sage
sage stickelberger_smith.sage
sage h_menos_pari.sage
sage h_menos_pari_grh.sage
sage orden_A22_69.sage
sage item2.sage
sage item3.sage
sage item3_ctrl.sage
sage item4.sage
sage item4b.sage 250 0
sage lat.sage 19 117 3 colon
```

`item4b.sage 250 40` takes about two minutes, `lat.sage 15 155 5 colon` and `lat.sage 22 225 5 colon`
about seven. The
lattice sweeps take hours and need several GB of memory per cell; the drivers show the exact
invocations, e.g. `CASILLAS="3:21" sage < factorizar_conductor.sage` for a single cell, and
`sage gor.sage list 100 20` for `gor_sweep_OUT.txt`.
