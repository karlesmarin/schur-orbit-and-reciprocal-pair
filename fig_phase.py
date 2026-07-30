"""
Figure: where one reciprocal pair stops being typical.

Authors: Carles Marin, Claude (AI assistant).

The evaluation at z_i = 1 is an element of order two, and it is the cheapest thing to compute about
the object. At one pair it decides everything: the object vanishes identically exactly when it
vanishes there, because it is +- a genuine character and the value at the identity is a dimension.
At two pairs and beyond the object is properly virtual, and the endpoint stops deciding: there are
shapes whose value at z = 1 is zero and whose character is not.

The figure counts those shapes. Each panel is one r; every partition with at most N parts is a point
at (|lambda|, ell(lambda)); the rare class -- vanishing at the endpoint but not identically -- is
drawn large and named, because it is the whole content. It is empty at r = 1 over 973 shapes and
non-empty from r = 2 on. The surviving implication, exact vanishing implies endpoint vanishing, is
what makes the endpoint a valid sieve and no more.

Data: fig_data_new.json, computed in exact arithmetic by fig_data_new.sage. Nothing is redrawn by
hand.

Palette: blue for the identically vanishing locus, orange for the endpoint-only class, neutral for
the rest. Marker shape doubles the colour for greyscale.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
ZERO = "#2a78d6"     # blue   -- identically zero
ENDP = "#eb6834"     # orange -- zero at the endpoint only
NONE = "#bfbeb6"     # neutral

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

D = json.load(open("fig_data_new.json"))["loci"]

fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.3))
for ax, key in zip(axes, ["r1", "r2", "r3"]):
    blk = D[key]
    N = blk["N"]; r = (N - 2) // 2
    rows = blk["rows"]
    nx, ny = [], []
    zx, zy = [], []
    ex, ey, elab = [], [], []
    for x in rows:
        if x["z_exact"]:
            zx.append(x["size"]); zy.append(x["ell"])
        elif x["z_end"]:
            ex.append(x["size"]); ey.append(x["ell"]); elab.append(x["lam"])
        else:
            nx.append(x["size"]); ny.append(x["ell"])
    ax.scatter(nx, ny, s=7, c=NONE, marker='.', linewidths=0, alpha=0.85, zorder=1)
    ax.scatter(zx, zy, s=20, facecolors='none', edgecolors=ZERO, marker='s',
               linewidths=0.85, zorder=3)
    ax.scatter(ex, ey, s=64, c=ENDP, marker='*', linewidths=0, zorder=5)
    # Autoscale fits the data, but the annotations are offset in POINTS, which autoscale cannot
    # foresee: a label on the rightmost point lands outside the axes. Give x explicit headroom and
    # flip the label inwards for points near the right edge.
    allx = nx + zx + ex
    maxx = max(allx) if allx else 1
    ax.set_xlim(-1, maxx + 1)
    # Several stars can share a row (at r=3 both (2,2,2,1,1) and (5,5,4,1,1) sit at ell=5), so a
    # fixed offset makes the labels collide. Alternate above and below, and anchor inwards near the
    # edges so no label leaves the axes.
    # At r=2 three of the four stars are bunched together, so no automatic offset keeps four labels
    # apart: they overprint each other and their markers. The full list of shapes is given in the
    # caption and in the accompanying remark, so the picture names one representative per panel --
    # the smallest by |lambda| -- and leaves the rest to the text.
    if ex:
        i = min(range(len(ex)), key=lambda j: (ex[j], ey[j]))
        ax.annotate("$(%s)$" % ",".join(str(t) for t in elab[i]), (ex[i], ey[i]),
                    textcoords="offset points", xytext=(0, 10), ha="center", va="bottom",
                    fontsize=7.0, color=ENDP, zorder=6)
    ax.set_title(r"$r=%d$,  $N=%d$:  %s" % (
        r, N, ("no such shape in %d" % len(rows)) if not ex else
        ("%d of %d shapes" % (len(ex), len(rows)))), fontsize=8.8, color=INK, pad=5)
    ax.set_xlabel(r"$|\lambda|$", color=INK)
    if key == "r1":
        ax.set_ylabel(r"$\ell(\lambda)$", color=INK)
    ax.set_ylim(-0.7, N + 1.4)          # headroom for the label placed above a star
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))   # |lambda| is an integer
    ax.grid(True, color=GRID, lw=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_color(MUTED); sp.set_linewidth(0.7)
    ax.tick_params(colors=SECOND, labelsize=7.6)

h = [plt.Line2D([], [], marker='s', markerfacecolor='none', markeredgecolor=ZERO, lw=0,
                markersize=5.0, markeredgewidth=0.85, label="vanishes identically"),
     plt.Line2D([], [], marker='*', color=ENDP, lw=0, markersize=8.5,
                label="vanishes at the endpoint only"),
     plt.Line2D([], [], marker='.', color=NONE, lw=0, markersize=6, label="neither")]
fig.legend(handles=h, loc='lower center', ncol=3, frameon=False, fontsize=7.8,
           bbox_to_anchor=(0.5, -0.035))
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig("fig_phase.pdf"); fig.savefig("fig_phase.png", dpi=170)
print("wrote fig_phase.pdf / .png")
for key in ["r1", "r2", "r3"]:
    rows = D[key]["rows"]
    print("  %s: %d shapes, %d identically zero, %d endpoint-only" % (
        key, len(rows), sum(1 for x in rows if x["z_exact"]),
        sum(1 for x in rows if x["z_end"] and not x["z_exact"])))
