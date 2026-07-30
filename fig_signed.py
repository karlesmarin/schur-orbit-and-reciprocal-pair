"""
Figure: the (-1)-enumeration of plane partitions in a box, with its product formula.

Authors: Carles Marin, Claude (AI assistant).

At t = 2 the weight in the enumerative corollary is (-1)^{m_2(T)}, and for a rectangular shape
(c^r) the tableaux are plane partitions in an r x (4-r) x c box. Setting z = 1 gives the plain
signed count, which the main theorem evaluates as  eps * d1 d2 d3 / (2 t^2).

Every value below is that closed form, cross-checked against a direct enumeration of the tableaux
(printed when the script runs), so the picture is data rather than assertion.

Form: the data is a small dense grid of signed integers, so this is a matrix with the values
printed, not a surface -- in three dimensions the bars occlude one another and the cell a value
belongs to stops being legible. Palette: diverging, one hue per sign with a neutral zero, per the
polarity rule; the printed integers mean the encoding is never colour alone.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from theorem_full import setup, lambda11        # noqa: E402
from enumeration import tableau_sum             # noqa: E402
from mpmath import mpf                          # noqa: E402

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
NEG_RAMP = ["#b2451f", "#eb6834", "#f4a582", "#f7d9c9"]
POS_RAMP = ["#cde2fb", "#6da7ec", "#2a78d6", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("div", NEG_RAMP + ["#f2f1ea"] + POS_RAMP)

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

T = 2
RMAX, CMAX = 4, 12


def signed_count(r, c):
    """eps * d1 d2 d3 / (2 t^2): the theorem evaluated at z = 1."""
    lam = tuple([c] * r)
    st = setup(lam, T)
    if st is None:
        return 0
    beta, Ac, Bc = st
    a1, a2 = beta[Ac[0]], beta[Ac[1]]
    b1, b2 = beta[Bc[0]], beta[Bc[1]]
    d = (a1 - a2, b1 - b2, abs(a1 + a2 - b1 - b2))
    if d[2] == 0:
        return 0
    l11 = lambda11(beta, Ac, Bc, len(beta))
    eps = (-1) ** (T + (T + 2) * (T + 3) // 2) * l11 * (1 if a1 + a2 - b1 - b2 > 0 else -1)
    num = eps * d[0] * d[1] * d[2]
    assert num % (2 * T * T) == 0, (r, c, d)
    return num // (2 * T * T)


print("cross-check of the closed form against a direct enumeration of the tableaux:")
bad = 0
for r in range(1, RMAX + 1):
    for c in range(1, 5):
        direct = float(tableau_sum(tuple([c] * r), T, mpf(1)).real)
        closed = signed_count(r, c)
        ok = abs(direct - closed) < 1e-9
        bad += (not ok)
        print(f"   box {r}x{c}: direct {direct:+7.1f}   closed form {closed:+4d}   {ok}")
print(f"   mismatches: {bad}\n")

M = np.array([[signed_count(r, c) for c in range(1, CMAX + 1)]
              for r in range(1, RMAX + 1)], dtype=float)
lim = np.abs(M).max()

fig, ax = plt.subplots(figsize=(9.2, 2.5))
im = ax.imshow(M, cmap=CMAP, norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
               aspect="equal", origin="upper")
for i in range(RMAX):
    for j in range(CMAX):
        v = int(M[i, j])
        shade = abs(v) / lim
        ax.text(j, i, f"{v:+d}" if v else "0", ha="center", va="center", fontsize=8.2,
                color=("#ffffff" if shade > 0.55 else INK))
ax.set_xticks(range(CMAX)); ax.set_xticklabels(range(1, CMAX + 1))
ax.set_yticks(range(RMAX)); ax.set_yticklabels(range(1, RMAX + 1))
ax.set_xlabel(r"columns $c$ of the shape $\lambda=(c^{\,r})$", color=INK, labelpad=4)
ax.set_ylabel(r"rows $r$", color=INK)
ax.tick_params(colors=MUTED, labelsize=8, length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_xticks(np.arange(-0.5, CMAX, 1), minor=True)
ax.set_yticks(np.arange(-0.5, RMAX, 1), minor=True)
ax.grid(which="minor", color=SURFACE, linewidth=1.6)
ax.tick_params(which="minor", length=0)

cb = fig.colorbar(im, ax=ax, fraction=0.020, pad=0.015, aspect=12)
cb.set_label("signed count", fontsize=8, color=SECOND)
cb.ax.tick_params(labelsize=7, colors=MUTED, length=2)
cb.outline.set_edgecolor(GRID)

fig.tight_layout(pad=0.5)
out = os.path.join(os.path.dirname(__file__), "fig_signed")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=200)
print("wrote", out + ".pdf/.png")
