"""
Figure: the image of the compression map  lambda -> (d1,d2,d3).

Authors: Carles Marin, Claude (AI assistant).

Every partition with at most t+2 rows is sent by the main theorem to three integers, so the
whole family collapses onto a sparse set of lattice points. The figure is built to make three
facts legible at once:

  (1) COMPRESSION  -- colour and area are the size of the fibre (how many partitions land on
      that one point), so the picture shows the collapse rather than asserting it;
  (2) d1, d2 are multiples of t -- visible as the coarse grid spacing on the two floor axes,
      and the spacing changes with t between the two panels;
  (3) d3 is never a multiple of t outside the overlapping profile -- the missing layers, and
      the vanishing plane d3 = 0 is drawn as the floor.

Palette: sequential single-hue blue ramp (validated), light -> dark = small -> large fibre;
the vanishing plane is the one warm accent, used nowhere else.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from fig_data import image  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECOND = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
ACCENT = "#eb6834"
ACCENT_DK = "#b2451f"
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seqblue", RAMP)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 9,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

MAXSIZE = 20
PANELS = [(3, MAXSIZE), (5, MAXSIZE)]

fig = plt.figure(figsize=(9.4, 4.3))
scs = []

for k, (t, M) in enumerate(PANELS):
    pts, nlam, nz = image(t, M)
    xs = np.array([p[0] for p in pts], dtype=float)
    ys = np.array([p[1] for p in pts], dtype=float)
    zs = np.array([p[2] for p in pts], dtype=float)
    ns = np.array([v["n"] for v in pts.values()], dtype=float)

    ax = fig.add_subplot(1, 2, k + 1, projection="3d")
    ax.set_facecolor(SURFACE)

    xhi, yhi, zhi = xs.max() + t * 0.6, ys.max() + t * 0.6, zs.max() * 1.06
    floor = [[(0, 0, 0), (xhi, 0, 0), (xhi, yhi, 0), (0, yhi, 0)]]
    ax.add_collection3d(Poly3DCollection(floor, facecolor=ACCENT, alpha=0.13,
                                         edgecolor=ACCENT, linewidths=0.9, zorder=0))

    for x, y, z in zip(xs, ys, zs):                      # stems: depth alone is unreadable
        ax.plot([x, x], [y, y], [0, z], color=GRID, lw=0.55, zorder=1)

    norm = Normalize(vmin=1, vmax=ns.max())
    sizes = 12 + 105 * np.sqrt(ns / ns.max())
    sc = ax.scatter(xs, ys, zs, c=ns, cmap=CMAP, norm=norm, s=sizes,
                    edgecolor=SURFACE, linewidth=0.7, depthshade=False, zorder=3)
    scs.append(sc)

    ax.set_xlim(0, xhi); ax.set_ylim(0, yhi); ax.set_zlim(0, zhi)
    ax.set_xticks([v for v in range(t, int(xhi) + 1, t)])
    ax.set_yticks([v for v in range(t, int(yhi) + 1, t)])
    ax.set_zticks([v for v in range(0, int(zhi) + 1, max(1, int(zhi) // 4))])
    ax.set_xlabel(r"$d_1$", labelpad=-2, color=INK, fontsize=10)
    ax.set_ylabel(r"$d_2$", labelpad=1, color=INK, fontsize=10)
    ax.set_zlabel(r"$d_3$", labelpad=-8, color=INK, fontsize=10, rotation=0)
    ax.tick_params(axis="x", colors=MUTED, labelsize=7, pad=-3)
    ax.tick_params(axis="y", colors=MUTED, labelsize=7, pad=-2)
    ax.tick_params(axis="z", colors=MUTED, labelsize=7, pad=-1)
    ax.view_init(elev=19, azim=-60)
    ax.set_box_aspect((1, 1, 0.62), zoom=1.18)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(SURFACE)
        axis.pane.set_edgecolor(GRID)
        axis._axinfo["grid"]["color"] = GRID
        axis._axinfo["grid"]["linewidth"] = 0.45
        axis.line.set_color(MUTED)
    fig.text(0.265 + 0.485 * k, 0.955,
             rf"$t={t}$:   {nz} partitions $\longrightarrow$ {len(pts)} points",
             ha="center", fontsize=10, color=INK)
    fig.text(0.265 + 0.485 * k, 0.905,
             rf"floor grid spacing $= t = {t}$",
             ha="center", fontsize=8, color=SECOND)

cax = fig.add_axes([0.365, 0.105, 0.27, 0.026])
cb = fig.colorbar(scs[0], cax=cax, orientation="horizontal")
cb.set_label("partitions collapsing onto one point", fontsize=8, color=SECOND, labelpad=4)
cb.ax.tick_params(labelsize=7, colors=MUTED, length=2)
cb.outline.set_edgecolor(GRID)

fig.subplots_adjust(left=0.0, right=1.0, top=1.02, bottom=0.17, wspace=0.0)
out = os.path.join(os.path.dirname(__file__), "fig_image")
fig.savefig(out + ".pdf")
fig.savefig(out + ".png", dpi=200)
print("wrote", out + ".pdf/.png")
