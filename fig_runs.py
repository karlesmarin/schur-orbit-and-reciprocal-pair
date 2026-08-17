"""
Figure: the words sigma_m(k) behind Proposition 6.4, and the involution read off them.

Authors: Carles Marin, Claude (AI assistant).

Proposition 6.4 is a statement about one word per k, and the word is invisible in the text. Drawn,
the three parts are one picture:

  (i)   every bar has height exactly +1 or -1 -- the terms are a SET, so the cancellation the
        enumeration needs is a matching and not something carrying multiplicity;
  (ii)  in the right panel every k is odd and no two bars are adjacent along m: the runs have
        length one and there is nothing to cancel;
  (iii) in the left panel every k is even and the colour alternates along each run, which is what
        makes "pair m with m+1" sign-reversing. Where a run ends, a gap ends it.

Palette: the same sequential blue as the other figures for +1, the one warm accent for -1, so that
the alternation is legible without a legend lookup.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from involution_runs import sequences, runs_of

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECOND = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BLUE = "#256abf"
ACCENT = "#eb6834"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 9,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

PANELS = [((10, 8), "every $k$ even: the runs alternate"),
          ((19, 2), "every $k$ odd: no two are adjacent")]

fig = plt.figure(figsize=(9.4, 4.2))

for panel, (lam, sub) in enumerate(PANELS):
    seqs = sequences(lam)
    ks = sorted(seqs)
    ax = fig.add_subplot(1, 2, panel + 1, projection="3d")
    ax.set_facecolor(SURFACE)

    nruns = 0
    for k in ks:
        seq = seqs[k]
        nruns += len(runs_of(seq))
        for m, s in enumerate(seq):
            if not s:
                ax.scatter([m], [k], [0], c=GRID, s=7, depthshade=False, zorder=2)
                continue
            col = BLUE if s > 0 else ACCENT
            ax.plot([m, m], [k, k], [0, s], color=col, lw=1.5, zorder=3)
            ax.scatter([m], [k], [s], c=col, s=26, edgecolor=SURFACE, linewidth=0.6,
                       depthshade=False, zorder=4)

    mmax = max(len(seqs[k]) for k in ks)
    ax.set_xlim(-0.5, mmax - 0.5)
    ax.set_ylim(min(ks) - 1, max(ks) + 1)
    ax.set_zlim(-1.25, 1.25)
    ax.set_xticks(list(range(0, mmax)))
    ax.set_xlabel(r"$m$", labelpad=-2, color=INK, fontsize=10)
    ax.set_ylabel(r"$k$", labelpad=2, color=INK, fontsize=10)
    ax.set_zlabel(r"$\sigma_m(k)$", labelpad=-2, color=INK, fontsize=9, rotation=0)
    ax.set_zticks([-1, 0, 1])
    ax.set_yticks(ks[::2] if len(ks) > 6 else ks)
    ax.tick_params(axis="x", colors=MUTED, labelsize=7, pad=-3)
    ax.tick_params(axis="y", colors=MUTED, labelsize=7, pad=-2)
    ax.tick_params(axis="z", colors=MUTED, labelsize=7, pad=-1)
    ax.view_init(elev=22, azim=-62)
    ax.set_box_aspect((1.25, 1, 0.5), zoom=1.15)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(SURFACE)
        axis.pane.set_edgecolor(GRID)
        axis._axinfo["grid"]["color"] = GRID
    ax.text2D(0.5, 0.97, r"$\lambda=(%s)$" % ",".join(str(p) for p in lam),
              transform=ax.transAxes, ha="center", color=INK, fontsize=11)
    ax.text2D(0.5, 0.90, sub, transform=ax.transAxes, ha="center",
              color=SECOND, fontsize=8.5)
    print("lambda=%-8s rows k: %d, runs: %d, all |sigma| = 1: %s"
          % (str(lam), len(ks), nruns,
             all(abs(s) <= 1 for k in ks for s in seqs[k])))

fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.03, wspace=0.04)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_runs")
fig.savefig(out + ".pdf")
fig.savefig(out + ".png", dpi=200)
print("wrote", out + ".pdf/.png")
