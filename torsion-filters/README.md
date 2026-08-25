# Torsion filters at a root-of-unity orbit — code and data

Companion material for *Schur polynomials twisted by roots of unity and reciprocal pairs: torsion
filters, fusion quotients, and total unimodularity at odd order* (Carles Marín),
[arXiv:2608.18302](https://arxiv.org/abs/2608.18302) ·
[doi:10.48550/arXiv.2608.18302](https://doi.org/10.48550/arXiv.2608.18302). This bundle is the one
that travels with **v2**, announced 24 August 2026; v1 was titled *…and an equal-rank reduction at
odd order*, and the identifier and the DOI are unchanged.

This is the second paper of a series. The first, *Schur polynomials twisted by roots of unity and
reciprocal pairs: exactly three factors, and where they vanish*
([arXiv:2608.09619](https://arxiv.org/abs/2608.09619)), settles the case of one free reciprocal
pair; this one takes up what happens with two or more.

## What is here

| directory | what it holds |
|---|---|
| `gates/` | the computations behind the claims — 266 scripts, 239 of them with the archived run that produced their numbers |
| `figures/` | the 15 scripts that draw the figures of the paper: the thirteen figures themselves, and the two that redraw them with Spanish labels |
| `audits/` | the 25 checks run on the manuscript itself rather than on the mathematics |
| `MANIFEST.md` | every script, one line on what it does, and whether its run is archived |

## How to read a claim

The paper's verification section states each measured claim with its range, its result, and the
**decoy** that was run against it — a deliberately wrong variant which had to fail, with its failure
count printed beside the claim. A statement whose decoy also passes is reported there as untested,
not as confirmed. The scripts carry the same discipline in their headers: each says what it tests,
what would falsify it, and what the decoy is.

To check a number, find the script in `MANIFEST.md`, read its header, and either read the archived
`*_OUT.txt` beside it or run it again.

## Running things

* Scripts ending in `.py` are plain Python 3 and need only the standard library. The exceptions are
  the figure scripts, which need `matplotlib`, and a few that use `pymupdf` to measure the rendered
  page.
* Scripts ending in `.sage` need SageMath. They were run in the `sage-normaliz` image; any recent
  Sage should do.
* Output files are named after their script: `foo.py` writes `foo_OUT.txt` and, where the result is
  machine-readable, `foo_DUMP.json`.

## What is deliberately not here

A handful of scripts stay in the working directory and do not travel: their subject is the *process*
— correspondence, and a third party's PDF-reading tool — rather than a statement of the paper. None
of them is the sole support of a number the paper quotes; that was checked against the source, not
assumed. The paper's own verification section quotes a gate count from the version in which they
still travelled, so it reads three higher than this bundle.

## What this material does not claim

The bundle is complete but it is not a per-row index: the paper's tables are organised by statement
and the scripts by computation, and the two groupings do not coincide. Twenty-four of the gate
scripts are helpers or exploratory probes with no archived run; none of them is the sole support of
a number quoted in the paper.

Most of the claims in the paper are measurements, and the paper says so throughout. Six statements
are conjectures, and the map in its §8.1 says which implies which and where we would attack. Nothing
in this bundle turns a measurement into a theorem.

## Licence and citation

Cite the paper. The scripts may be used freely; they are research code, written to be read and
doubted rather than to be depended on.
