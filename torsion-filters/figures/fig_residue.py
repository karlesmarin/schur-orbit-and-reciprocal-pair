"""
Figure: where the open problem actually is, and the discriminants that do not find it.

Authors: Carles Marin, Claude (AI assistant).

Three panels, and each answers a different question about the residual vanishing locus.

LEFT -- how much is left. A funnel on a logarithmic scale: of all the forms in range, the classical
core criterion accounts for most of the vanishing, the reduction to the unspecialised object accounts
for a few more, and what survives both is a set small enough to list. The point of the panel is the
last bar: the open problem is not "most of the zeros", it is a dozen forms.

CENTRE -- the strongest signal we have, and its limit. Every t-core carried by the occupied
population, ordered by size, with the vanishing forms stacked on top. Thirty-three of the
thirty-six cores contain no vanishing form at all, which is a real and sharp restriction; but the
three that do are far from pure, which is why the core is a necessary condition and not a criterion.

RIGHT -- the best discriminant, and the stratum where it dies. Inside one core the split is perfect:
the six vanishing forms are exactly the six whose quotient component at the self-paired residue is
non-empty. Inside the next core the same condition holds for all the vanishing forms and for eight
of the twenty others. Both strata are drawn, at the same scale, so the failure is as visible as the
success.

Everything is computed from the populations; the counts printed are the counts measured.

Palette: the two validated categorical hues (blue = vanishing, orange = the failed discriminant),
neutrals for the population. Hatching doubles the colour so the figure survives greyscale printing.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
ZERO = "#2a78d6"      # blue   -- vanishing
FAIL = "#eb6834"      # orange -- the discriminant's false positives
POP = "#d6d5cf"       # the population

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

D = json.load(open("quotient_split_DUMP.json"))          # t=6, r=2, W=13: 491 ocupadas, 12 nulas

fig = plt.figure(figsize=(9.6, 3.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.0], wspace=0.33)

# ------------------------------------------------------------------ IZQUIERDA: el embudo --------
ax = fig.add_subplot(gs[0, 0])
# los TRES paneles sobre la MISMA poblacion, t=6: mezclar t=4 en uno solo y t=6 en los otros dos
# hacia que las cifras no se pudieran comparar de un panel a otro.
ETAPAS = [("all forms", 715, POP),
          ("vanishing", 236, MUTED),
          (r"not seen by the $t$-core", 12, ZERO)]
ys = np.arange(len(ETAPAS))[::-1]
for y, (lab, n, col) in zip(ys, ETAPAS):
    ax.barh(y, n, height=0.62, color=col, edgecolor=MUTED, linewidth=0.5, zorder=3)
    ax.text(n * 1.16, y, "%d" % n, va="center", fontsize=8.5,
            color=INK if col is not POP else SECOND)
ax.set_yticks(ys); ax.set_yticklabels([e[0] for e in ETAPAS], fontsize=8)
ax.set_xscale("log"); ax.set_xlim(7, 2000)
ax.set_xlabel(r"forms ($t=6$, $r=2$, $\beta_i\leq13$)")
ax.set_title(r"how much is left", fontsize=9, color=INK)
ax.grid(True, axis="x", color=GRID, linewidth=0.5, zorder=0)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.6)

# ------------------------------------------------------------------ CENTRO: por t-core ----------
ax = fig.add_subplot(gs[0, 1])
porcore = {}
for x in D:
    c = tuple(x["core"])
    a, b = porcore.get(c, (0, 0))
    porcore[c] = (a + 1, b + int(x["nula"]))
orden = sorted(porcore.items(), key=lambda kv: (-kv[1][1], -kv[1][0]))
X = np.arange(len(orden))
tot = [v[0] for _, v in orden]
nul = [v[1] for _, v in orden]
ax.bar(X, tot, color=POP, edgecolor=MUTED, linewidth=0.4, zorder=3, width=0.82)
ax.bar(X, nul, color=ZERO, edgecolor=ZERO, linewidth=0.4, zorder=4, width=0.82)
con = sum(1 for n in nul if n)
ax.axvline(con - 0.5, color=FAIL, linewidth=1.0, linestyle=(0, (4, 3)), zorder=6)
ax.text(con + 1.2, max(tot) * 1.045, r"%d cores with no vanishing form at all" % (len(orden) - con),
        fontsize=8, color=SECOND, va="bottom")
ax.set_ylim(0, max(tot) * 1.22)
for i in range(con):
    ax.text(i, tot[i] + 0.6, r"$%d/%d$" % (nul[i], tot[i]), ha="center", fontsize=7.5, color=ZERO)
ax.set_xticks([]); ax.set_xlabel(r"the %d $t$-cores of the occupied population" % len(orden))
ax.set_ylabel(r"forms")
ax.set_title(r"the vanishing concentrates on three cores", fontsize=9, color=INK)
ax.grid(True, axis="y", color=GRID, linewidth=0.5, zorder=0)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.6)

# ------------------------------------------------------------------ DERECHA: el separador -------
ax = fig.add_subplot(gs[0, 2])
ESTRATOS = [((1, 1, 1), "core $(1,1,1)$"), ((2, 1), "core $(2,1)$")]
w = 0.34
for k, (C, lab) in enumerate(ESTRATOS):
    est = [x for x in D if tuple(x["core"]) == list(C) or tuple(x["core"]) == C]
    nz = [x for x in est if x["nula"]]
    no = [x for x in est if not x["nula"]]
    a = sum(1 for x in nz if x["quot"][0])          # nulas con la componente q=0 no vacia
    b = sum(1 for x in no if x["quot"][0])          # NO nulas con la misma condicion: los falsos
    ax.bar(k - w / 2, len(nz), width=w, color=POP, edgecolor=MUTED, linewidth=0.5, zorder=3)
    ax.bar(k - w / 2, a, width=w, color=ZERO, edgecolor=ZERO, linewidth=0.5, zorder=4)
    ax.bar(k + w / 2, len(no), width=w, color=POP, edgecolor=MUTED, linewidth=0.5, zorder=3)
    ax.bar(k + w / 2, b, width=w, color=FAIL, edgecolor=FAIL, linewidth=0.5, zorder=4,
           hatch="///")
    ax.text(k - w / 2, len(nz) + 0.35, r"$%d/%d$" % (a, len(nz)), ha="center", fontsize=8, color=ZERO)
    ax.text(k + w / 2, len(no) + 0.35, r"$%d/%d$" % (b, len(no)), ha="center", fontsize=8,
            color=FAIL if b else SECOND)
ax.set_xticks(range(len(ESTRATOS)))
ax.set_xticklabels([e[1] for e in ESTRATOS], fontsize=8.5)
ax.set_ylabel(r"forms")
ax.set_title(r"$\lambda^{(0)}\neq\varnothing$: perfect, then not", fontsize=9, color=INK)
ax.set_ylim(0, 24)
ax.grid(True, axis="y", color=GRID, linewidth=0.5, zorder=0)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.6)

leyenda = [Patch(facecolor=ZERO, edgecolor=ZERO, label=r"vanishing forms"),
           Patch(facecolor=FAIL, edgecolor=FAIL, hatch="///", label=r"false positives of the discriminant"),
           Patch(facecolor=POP, edgecolor=MUTED, label=r"the rest of the population")]
fig.legend(handles=leyenda, frameon=False, fontsize=8, ncol=3,
           loc="lower center", bbox_to_anchor=(0.5, -0.10))

fig.savefig("fig_residue.pdf", bbox_inches="tight", pad_inches=0.02)
print("fig_residue.pdf escrito")
print("  cores: %d, con alguna nula: %d" % (len(orden), con))
for c, (n, z) in orden[:4]:
    print("    %-14s %3d formas, %2d nulas" % (str(c), n, z))
