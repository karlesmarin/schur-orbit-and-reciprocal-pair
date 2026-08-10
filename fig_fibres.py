"""
Figure: the fibres of the evaluation invariant  I_t(lambda) = (d1,d2,d3,eps).

Authors: Carles Marin, Claude (AI assistant).

The companion of fig_image.py, and its complement.  That figure plots what the invariant RECORDS:
the set of triples actually taken.  This one plots what the invariant FORGETS: over each cell of the
(d1,d2) floor stands the column of all partition SIZES that land there, so a column is a set of
partitions of different sizes carrying one and the same value of Phi_t.

Reading the picture:

  (1) a vertical column is one fibre -- every marker in it has the same value, so the height of the
      column is the range of sizes the alphabet cannot distinguish;
  (2) the columns do not thin out as |lambda| grows, which is the compression getting stronger, not
      weaker, with size: the number of values is bounded by the lattice, the number of partitions is
      not;
  (3) the accent column is the fibre of the empty partition.  At t=3 it is I=(3,3,2,+1), and it
      contains lambda = () together with shapes of every size up to the range: Phi_3 cannot tell the
      empty partition from a partition of 20 boxes.

Several d3 can share one (d1,d2) cell, so colour carries d3 and a cell may hold more than one fibre;
the caption says so, and the accent column is a single fibre, sign included.

Palette: the sequential blue ramp of fig_image.py, here indexed by d3, with the same single warm
accent -- reused for the one distinguished fibre rather than for the vanishing plane.
"""
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

from law_control import partitions
from theorem_full import setup, lambda11, phi_bialternant

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
THETA = 0.3714159          # a generic point: no sinh in the formula vanishes there


def fibres(t, maxsize):
    """invariant (d1<=d2, d3, eps) -> sorted list of sizes |lambda| landing on it.

    Also returns, per invariant, the set of values Phi_t actually takes on that fibre, evaluated
    from the bialternant at a generic theta.  The picture asserts that a column is ONE value; that
    assertion is a computation, so it is made here rather than in the caption.
    """
    N = t + 2
    out = defaultdict(list)
    vals = defaultdict(list)
    for n in range(0, maxsize + 1):
        for lam in partitions(n, N):
            st = setup(lam, t)
            if st is None:
                continue
            beta, Ac, Bc = st
            a1, a2 = beta[Ac[0]], beta[Ac[1]]
            b1, b2 = beta[Bc[0]], beta[Bc[1]]
            d3 = abs(a1 + a2 - b1 - b2)
            if d3 == 0:                       # the concentric locus: the value is zero, no fibre
                continue
            l11 = lambda11(beta, Ac, Bc, N)
            eps = (-1) ** (t + (t + 2) * (t + 3) // 2) * l11 * (1 if a1 + a2 - b1 - b2 > 0 else -1)
            d1, d2 = sorted((a1 - a2, b1 - b2))
            out[(d1, d2, d3, eps)].append(n)
            vals[(d1, d2, d3, eps)].append(phi_bialternant(lam, t, THETA))
    return {k: sorted(v) for k, v in out.items()}, dict(vals)


def empty_fibre(t):
    """the invariant of lambda = (), the fibre singled out in the accent colour."""
    N = t + 2
    beta, Ac, Bc = setup((), t)
    a1, a2 = beta[Ac[0]], beta[Ac[1]]
    b1, b2 = beta[Bc[0]], beta[Bc[1]]
    l11 = lambda11(beta, Ac, Bc, N)
    eps = (-1) ** (t + (t + 2) * (t + 3) // 2) * l11 * (1 if a1 + a2 - b1 - b2 > 0 else -1)
    d1, d2 = sorted((a1 - a2, b1 - b2))
    return (d1, d2, abs(a1 + a2 - b1 - b2), eps)


fig = plt.figure(figsize=(11.6, 4.9))
REPORT = []

for k, (t, M) in enumerate(PANELS):
    F, VALS = fibres(t, M)
    star = empty_fibre(t)
    ax = fig.add_subplot(1, 2, k + 1, projection="3d")
    ax.set_facecolor(SURFACE)

    d3max = max(key[2] for key in F)
    norm = Normalize(vmin=1, vmax=d3max)

    # every fibre is a vertical stem over its (d1,d2) cell; the stem carries one marker per size
    for key, sizes in sorted(F.items()):
        d1, d2, d3, eps = key
        is_star = (key == star)
        col = ACCENT if is_star else CMAP(norm(d3))
        z0, z1 = sizes[0], sizes[-1]
        ax.plot([d1, d1], [d2, d2], [z0, z1],
                color=ACCENT_DK if is_star else col,
                lw=1.9 if is_star else 0.55,
                alpha=1.0 if is_star else 0.42, zorder=6 if is_star else 2)
        ax.scatter([d1] * len(sizes), [d2] * len(sizes), sizes,
                   s=26 if is_star else 5.5, c=[col] * len(sizes),
                   edgecolors=ACCENT_DK if is_star else "none",
                   linewidths=0.5 if is_star else 0,
                   depthshade=False, zorder=7 if is_star else 3)

    # d1 and d2 are multiples of t (Proposition on the quotient), so the floor ticks are put on the
    # lattice rather than left to the auto-locator: the spacing is a fact about the data
    x1 = max(key[0] for key in F)
    y1 = max(key[1] for key in F)
    ax.set_xticks(list(range(t, x1 + 1, t)))
    ax.set_yticks(list(range(t, y1 + 1, t)))
    ax.set_xlim(0, x1 + t * 0.6)
    ax.set_ylim(0, y1 + t * 0.6)
    ax.set_zlim(-1, M + 1)
    ax.set_xlabel(r"$d_1$", labelpad=1, color=SECOND)
    ax.set_ylabel(r"$d_2$", labelpad=1, color=SECOND)
    ax.set_zlabel(r"$|\lambda|$", labelpad=1, color=SECOND)
    ax.tick_params(labelsize=7, colors=MUTED, pad=0.6)
    ax.view_init(elev=17, azim=-58)
    ax.set_box_aspect((1, 1, 0.72))
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(SURFACE)
        axis.pane.set_edgecolor(GRID)
        axis._axinfo["grid"]["color"] = GRID
        axis._axinfo["grid"]["linewidth"] = 0.45
        axis.line.set_color(MUTED)

    sizes_star = F[star]
    fig.text(0.265 + 0.485 * k, 0.955,
             rf"$t={t}$:   {len(F)} fibres over $|\lambda|\leq {M}$",
             ha="center", fontsize=10, color=INK)
    fig.text(0.265 + 0.485 * k, 0.905,
             rf"accent: $I_{{{t}}}=({star[0]},{star[1]},{star[2]},{'+' if star[3] > 0 else '-'}1)$,"
             rf" {len(sizes_star)} partitions of sizes {sizes_star[0]}–{sizes_star[-1]}",
             ha="center", fontsize=8, color=ACCENT_DK)

    REPORT.append((t, M, F, star, VALS))

sm = plt.cm.ScalarMappable(cmap=CMAP, norm=Normalize(vmin=1, vmax=max(
    key[2] for _, _, F, _, _ in REPORT for key in F)))
sm.set_array([])
cax = fig.add_axes([0.365, 0.085, 0.27, 0.026])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_label(r"the coupling term $d_3$", fontsize=8, color=SECOND, labelpad=4)
cb.ax.tick_params(labelsize=7, colors=MUTED, length=2)
cb.outline.set_edgecolor(GRID)

fig.subplots_adjust(left=0.0, right=1.0, top=1.02, bottom=0.155, wspace=0.0)
out = os.path.join(os.path.dirname(__file__), "fig_fibres")
fig.savefig(out + ".pdf")
fig.savefig(out + ".png", dpi=200)

# every number the caption quotes, printed so it can be archived and checked
TOL = 1e-7
for t, M, F, star, VALS in REPORT:
    sizes = [len(v) for v in F.values()]
    spans = [v[-1] - v[0] for v in F.values()]
    s = F[star]
    print("t=%d |lam|<=%d: fibres %d ; singletons %d ; largest fibre %d ; mean fibre %.2f"
          % (t, M, len(F), sum(1 for x in sizes if x == 1), max(sizes), sum(sizes) / len(sizes)))
    # `sizes %d..%d` is a RANGE, not a claim that every size in it occurs, and the two are not the
    # same statement: at t=3 the sizes 1 and 4 are absent from this fibre.  Print the gap, so that a
    # sentence about it can be checked against the run rather than read off the endpoints.
    absent = [n for n in range(s[0], s[-1] + 1) if n not in set(s)]
    print("t=%d the empty fibre I=(%d,%d,%d,%+d): %d partitions, sizes %d..%d, "
          "%d distinct sizes, absent: %s"
          % (t, star[0], star[1], star[2], star[3], len(s), s[0], s[-1],
             len(set(s)), absent if absent else "none"))
    print("t=%d largest size-span in one fibre: %d ; mean span %.2f"
          % (t, max(spans), sum(spans) / len(spans)))

    # CONTROL 1.  A column is one value.  Within every fibre the bialternant must be constant.
    worst = max(max(abs(v - vs[0]) for v in vs) for vs in VALS.values())
    nlam = sum(len(v) for v in VALS.values())
    print("t=%d control, the value is constant on each fibre: %d partitions, worst spread %.2e -> %s"
          % (t, nlam, worst, "0 fail" if worst < TOL else "*** FAIL ***"))

    # CONTROL 2.  A control that can fail: are distinct invariants distinct values?  They are not,
    # and the reason is legible -- the closed form is symmetric in d1, d2, d3, so the value sees
    # only the MULTISET and the sign.  Counting the three classes separates the three notions.
    def cl(x):
        return round(x.real, 7), round(x.imag, 7)
    distinct_vals = {cl(vs[0]) for vs in VALS.values()}
    multisets = {(tuple(sorted(k[:3])), k[3]) for k in VALS}
    print("t=%d fibres %d ; distinct multisets+sign %d ; distinct values %d"
          % (t, len(F), len(multisets), len(distinct_vals)))
    print("t=%d the invariant is complete but not minimal: %s"
          % (t, "values = multisets" if len(distinct_vals) == len(multisets)
             else "values %d != multisets %d" % (len(distinct_vals), len(multisets))))
print("wrote", out + ".pdf/.png")
