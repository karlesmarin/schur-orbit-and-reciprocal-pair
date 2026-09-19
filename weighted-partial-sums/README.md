# 🔢 Zeros of weighted partial sums of completely multiplicative functions — code and data

Companion material for *Zeros of weighted partial sums of completely multiplicative functions:
Dirichlet characters, the Euler factor at 2 and Conrey's sine series* (Carles Marín), deposited on
Zenodo: [doi:10.5281/zenodo.22847374](https://doi.org/10.5281/zenodo.22847374). That is the
**concept DOI** — it always resolves to the current version — and it carries the note in both
editions, English (19 pp) and Spanish (20 pp).

Both PDFs are here too — `segmento_II.pdf` and `segmento_II_es.pdf`, byte for byte the files
deposited on Zenodo — together with the tool, the scripts behind every number the note calls
*measured*, their archived output, and the Lean file.

It lives in this repository, next to [`segmento/`](../segmento) and [`conductor/`](../conductor),
because the question it answers came out of them: which components of the moment matrix of the ring
of character values vanish on its diagonal.

## 🎯 What the note proves

For a completely multiplicative $\varepsilon:\mathbb N\to\{\pm1\}$ put
$S_\varepsilon(m)=\sum_{n\le m}n\,\varepsilon(n)$. Parity forces $m\equiv0,3\pmod4$ at any zero.

1. **Every admissible length.** For every $m\ge3$ with $m\equiv0$ or $3\pmod4$ there is a
   completely multiplicative $\varepsilon$ with $S_\varepsilon(m)=0$ — by 104 interval certificates
   in integer arithmetic up to $10^8$, and by Sárközy's arithmetic progressions in subset sums
   beyond.
2. **How often one function vanishes.** Every finite prefix of signs extends to a zero at all large
   admissible $m$: the functions with infinitely many zeros form a dense $G_\delta$ of probability
   zero.
3. **Which functions do.** Quadratic characters give explicit points of that set, through a digit
   recursion
   $S_\chi(qk+r)=s\,q\,S_\chi(k)+k\bigl(T_q+q\,C(r)\bigr)+T(r)$ that propagates their zeros in base
   $q$. For odd characters the digits that propagate zeros are exactly the intervals on which
   Conrey's sine series of the Legendre symbol vanishes identically.
4. **Where the question came from.** For every odd Dirichlet character modulo an odd $f$,
   $\sum_{c\le(f-1)/2}c\,\chi(c)=\frac f2B_{1,\chi}\bigl(\chi(2)^{-1}-1\bigr)$: the Euler factor at
   $2$ decides which components vanish on the diagonal.

## 📦 What is here

| | |
|---|---|
| `segmento_II.pdf`, `segmento_II_es.pdf` | the note, English and Spanish, as deposited |
| `segmento_II.py` | the tool: one file, Python standard library only, exact integers |
| `SegmentoII.lean` | the digit recursion, the propagation of its zeros and the pairs lemma, in Lean 4 / Mathlib |
| `II_*.py`, `*.sage`, `*_OUT.txt`, `II_certificados_A_bloques.csv` | the scripts behind every *measured* number, each with its archived run, and the 104 block certificates |
| `figuras_II.py` | draws the five figures, in English or Spanish (`python figuras_II.py es`) |

[`TOOL_README.md`](TOOL_README.md) — the README of the ancillary archive — maps each script to the
statement it checks.

## 🧪 How to doubt a number

```
python segmento_II.py autotest
```

recomputes from nothing the worked examples of the note — the zero cells of $q=103$, the empty
digit set of $q=59$, the loop of $q=131$, the families at $(q^j-1)/2$, a sample of certificates —
exercises every guard, and ends with `SELF-TEST: all correct`. `python segmento_II.py certificado`
re-verifies all 104 blocks in integer arithmetic. Every other script says in its header what it
computes; read the `*_OUT.txt` beside it or run it again. The decoys are meant to fail, and do.

## 🚫 What is deliberately not here

The correspondence around this line of work does not travel with it, held out **by location, not by
a list of names**. A check refuses to release a bundle in which anything imports or reads a file
that is not in it; this bundle passes it.

## 📄 Licence and citation

Cite the note — concept DOI [10.5281/zenodo.22847374](https://doi.org/10.5281/zenodo.22847374). The
note is under CC BY 4.0; the code, the data and the Lean file under CC0 1.0.

Carles Marín, independent researcher — `karlesmarin@gmail.com` —
[ORCID 0009-0007-5637-9688](https://orcid.org/0009-0007-5637-9688). Claude (Anthropic) was used as a
research assistant; every statement was checked by the author, who takes full responsibility.
