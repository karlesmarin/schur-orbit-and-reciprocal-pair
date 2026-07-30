"""
Figure: the independence locus of Theorem 5.1, computed rather than asserted.

Authors: Carles Marin, Claude (AI assistant).

For every two-row partition in range we TEST whether Phi_t(lambda;z) = s_lambda(z,1/z) by
evaluating both sides, and then colour the point by which of the two families of the theorem it
belongs to. The odd-t panel is the control: no extra family can exist there, and none appears.

Palette: the two validated categorical hues (blue = t-cores, the family already covered by
[AK25, Thm 5.3]; orange = the extra family created by the reciprocal specialization), plus a
neutral for the partitions where independence fails. Marker shape doubles the colour so the
figure survives greyscale printing.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpmath import mp, mpf, sinh

from theorem_full import closed_form          # noqa: E402
from ak53_consistency import core             # noqa: E402

mp.dps = 30
THETA = mpf("0.41")

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
CORE = "#2a78d6"      # blue   -- t-cores
EXTRA = "#eb6834"     # orange -- the extra family
NONE = "#d6d5cf"      # neutral

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

PANELS = [3, 4, 6]
MAXL1 = 26

fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.5))

for ax, t in zip(axes, PANELS):
    xs_c, ys_c, xs_e, ys_e, xs_n, ys_n = [], [], [], [], [], []
    n_c = n_e = 0
    for l1 in range(0, MAXL1 + 1):
        for l2 in range(0, l1 + 1):
            lam = tuple(v for v in (l1, l2) if v > 0)
            val = closed_form(lam, t, THETA)
            ref = sinh((l1 - l2 + 1) * THETA) / sinh(THETA)
            indep = abs(abs(val) - ref) < mpf("1e-18")
            if not indep:
                xs_n.append(l2); ys_n.append(l1); continue
            if core(lam, t, t + 2) == lam:
                xs_c.append(l2); ys_c.append(l1); n_c += 1
            else:
                xs_e.append(l2); ys_e.append(l1); n_e += 1

    ax.scatter(xs_n, ys_n, s=5, c=NONE, marker=".", linewidths=0, zorder=1)
    ax.scatter(xs_c, ys_c, s=30, c=CORE, marker="o", edgecolor=SURFACE, linewidth=0.5,
               zorder=3, label=rf"$t$-cores  (${n_c}$)")
    ax.scatter(xs_e, ys_e, s=52, c=EXTRA, marker="^", edgecolor=SURFACE, linewidth=0.6,
               zorder=4, label=rf"extra family  (${n_e}$)")

    ax.set_title(rf"$t={t}$" + ("   (odd: no extra family)" if t % 2 else ""),
                 fontsize=10, color=INK, pad=6)
    ax.set_xlabel(r"$\lambda_2$", color=INK)
    if t == PANELS[0]:
        ax.set_ylabel(r"$\lambda_1$", color=INK)
    hi = max(ys_c + ys_e) + 2
    ax.set_xlim(-0.9, max(xs_c + xs_e) + 2.2); ax.set_ylim(-0.9, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(colors=MUTED, labelsize=7.5)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.grid(True, color=GRID, linewidth=0.4, zorder=0)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(loc="lower right", fontsize=7.2, frameon=True, facecolor=SURFACE,
              edgecolor=GRID, labelcolor=SECOND, handletextpad=0.4, borderpad=0.4)

fig.tight_layout(pad=0.7)
out = os.path.join(os.path.dirname(__file__), "fig_locus")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=200)
print("wrote", out + ".pdf/.png")
