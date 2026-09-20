# -*- coding: utf-8 -*-
r"""columnas.py -- the columns of SU(n)_{p-n} by the semigroup of their spectrum (the data of Figure 2).

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

For p = 17 and 19 and every 2 <= n <= p-2, every centred spectrum T (|T| = n) is classified by its
exact semigroup S_T (P3a.s_exacto): N_0 (maximal), <2,3>, other symmetric, not symmetric (not
Gorenstein); each spectrum stands for n columns.  Checks, in exact arithmetic:
  * level-rank: the proportions for n and p-n agree class by class (cross products);
  * the maximal counts against the closed formula N_0(p, n);
  * the non-symmetric fraction at n = 8, p = 17, against 1/17.
  python columnas.py          the table and the checks
  python columnas.py pdf      also draws fig_columnas.pdf and fig_columnas_es.pdf (needs matplotlib)
Public domain (CC0).
"""
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations

from P3a import maximales, s_exacto, semigrupo

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
CLASES = ("N0", "<2,3>", "sim", "nosim")


def clase(T, p):
    S = semigrupo(s_exacto(T, p))
    if not S["huecos"]:
        return "N0"
    if S["huecos"] == [1]:
        return "<2,3>"
    return "sim" if S["simetrico"] else "nosim"


def censo(p):
    d = {}
    for n in range(2, p - 1):
        c = Counter()
        for T in combinations(range(p), n):
            if sum(T) % p == 0:
                c[clase(T, p)] += n
        d[n] = c
    return d


def dibuja(datos):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colores = {"N0": "#1F4E79", "<2,3>": "#4E8C4A", "sim": "#C9B458", "nosim": "#B03030"}
    textos = {"en": {"N0": "maximal ($S_T=\\mathbb{N}_0$)", "<2,3>": "$\\langle2,3\\rangle$",
                     "sim": "other symmetric", "nosim": "not symmetric (not Gorenstein)",
                     "x": "$n$  (level $k=p-n$)", "y": "fraction of columns"},
              "es": {"N0": "maximal ($S_T=\\mathbb{N}_0$)", "<2,3>": "$\\langle2,3\\rangle$",
                     "sim": "otros simétricos", "nosim": "no simétrico (no Gorenstein)",
                     "x": "$n$  (nivel $k=p-n$)", "y": "fracción de columnas"}}
    for idioma, sufijo in (("en", ""), ("es", "_es")):
        tx = textos[idioma]
        fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.8))
        for ax, p in zip(axes, sorted(datos)):
            d = datos[p]
            ns = sorted(d)
            abajo = [0.0] * len(ns)
            for c in CLASES:
                vals = [d[n][c] / sum(d[n].values()) for n in ns]
                ax.bar(ns, vals, bottom=abajo, color=colores[c], width=0.8, label=tx[c],
                       edgecolor="white", linewidth=0.4)
                abajo = [a + v for a, v in zip(abajo, vals)]
            ax.set_xticks(ns)
            ax.set_xlabel(tx["x"])
            ax.set_title("$\\mathrm{SU}(n)_{%d-n}$, $p=%d$" % (p, p), fontsize=10)
            ax.set_ylim(0, 1)
        axes[0].set_ylabel(tx["y"])
        axes[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
        fig.tight_layout()
        fig.savefig("fig_columnas%s.pdf" % sufijo, bbox_inches="tight")
        print("written fig_columnas%s.pdf" % sufijo)


datos = {}
for p in (17, 19):
    d = censo(p)
    datos[p] = d
    sim = all(d[n][c] * sum(d[p - n].values()) == d[p - n][c] * sum(d[n].values())
              for n in d for c in CLASES)
    formula = all(d[n]["N0"] == maximales(p, n) for n in d)
    print("p=%d  level-rank (n <-> p-n) proportions equal: %s ; maximal counts = N_0(p, n): %s"
          % (p, sim, formula))
    for n in sorted(d):
        print("   n=%2d %s" % (n, {c: d[n][c] for c in CLASES if d[n][c]}))
    sys.stdout.flush()
d = datos[17][8]
fr = Fraction(d["nosim"], sum(d.values()))
print("non-symmetric fraction at n=8, p=17: %s = %.4f  (1/17 = %.4f)" % (fr, float(fr), 1 / 17))
if "pdf" in sys.argv[1:]:
    dibuja(datos)
