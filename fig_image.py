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
PANELS = [(3, MAXSIZE), (4, MAXSIZE)]
ZEROS = {}          # concentric shapes per panel: they live ON the plane d3 = 0
DRAWN = {}          # the compression the caption quotes: partitions in, points out

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

    # the concentric shapes, which vanish and so are not in the image of the compression map:
    # they are exactly the points of the plane d3 = 0, and by Proposition 3.5 there are none
    # at all when t is odd.  Drawn in the accent, the one warm colour, and counted in ZEROS.
    zpts, _, _ = image(t, M, include_zero=True)
    zk = [p for p in zpts if p[2] == 0]
    ZEROS[t] = (len(zk), sum(zpts[p]["n"] for p in zk))
    DRAWN[t] = (nz, len(pts))
    if zk:
        ax.scatter([p[0] for p in zk], [p[1] for p in zk], [0] * len(zk),
                   c=ACCENT_DK, s=26, marker="D", edgecolor=SURFACE, linewidth=0.6,
                   depthshade=False, zorder=4)

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
for t_, (nz_, npt_) in sorted(DRAWN.items()):
    print("t=%d: the compression drawn: %d partitions land on %d points" % (t_, nz_, npt_))
for t_, (npts, nlam_) in sorted(ZEROS.items()):
    print("t=%d: concentric points on the plane d3=0: %d (from %d partitions)"
          % (t_, npts, nlam_))
# the two counts the caption quotes.  "Without the identification" means the number of ORDERED
# triples (d1,d2,d3) actually taken -- not the symmetrised set closed up under the exchange, which
# would count triples nothing maps to.
#
# WHICH ORDER.  Theorem 3.1 fixes it: "let r_A <= r_B be the residues involved".  theorem_full.setup
# returns the two classes sorted by COLUMN, which is a different rule, and the two rules give
# different counts -- 108/80 by column against 126/92 by residue.  The value does not care, since
# the closed form is symmetric in d1 and d2, but this count does, so it has to be taken in the
# order the theorem states.  Reordering here rather than in setup() keeps the sign code, which is
# order-invariant by construction, untouched.
from theorem_full import setup                                   # noqa: E402
from law_control import partitions as _parts                     # noqa: E402
for t_, M_ in PANELS:
    pts_, _, _ = image(t_, M_)
    ordered, ordered_by_column = set(), set()
    for n_ in range(0, M_ + 1):
        for lam_ in _parts(n_, t_ + 2):
            st_ = setup(lam_, t_)
            if st_ is None:
                continue
            b_, Ac_, Bc_ = st_
            byc_ = (b_[Ac_[0]] - b_[Ac_[1]], b_[Bc_[0]] - b_[Bc_[1]],
                    abs(b_[Ac_[0]] + b_[Ac_[1]] - b_[Bc_[0]] - b_[Bc_[1]]))
            if b_[Ac_[0]] % t_ > b_[Bc_[0]] % t_:                # impose r_A <= r_B
                Ac_, Bc_ = Bc_, Ac_
            a1_, a2_, b1_, b2_ = b_[Ac_[0]], b_[Ac_[1]], b_[Bc_[0]], b_[Bc_[1]]
            d_ = (a1_ - a2_, b1_ - b2_, abs(a1_ + a2_ - b1_ - b2_))
            if d_[2]:
                ordered.add(d_)
                ordered_by_column.add(byc_)
    print("t=%d: points up to the exchange %d, ordered triples actually taken %d"
          % (t_, len(pts_), len(ordered)))
    print("t=%d: CONTROL, the same count with A and B ordered by column instead: %d"
          % (t_, len(ordered_by_column)))
print("wrote", out + ".pdf/.png")
