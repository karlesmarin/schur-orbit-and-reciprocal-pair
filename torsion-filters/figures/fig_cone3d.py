"""
Figure: the interlacing cone, and the checkerboard inside it.

Authors: Carles Marin, Claude (AI assistant).

This is the second figure that earns three dimensions, and it earns them for a plain reason: at rank
two the whole domain of the operator IS three-dimensional. A pair (Lambda, mu) with Lambda of two
parts and mu of one is a point of (Lambda_1, Lambda_2, mu_1), the double interlacing carves a cone
out of that space, and the parity condition selects a sublattice inside the cone. Both statements
that the text proves are therefore visible at once, on the same picture.

  the CONE      Lambda_1 >= mu_1 >= 0 and the interlacing bound: the wedge drawn in outline;
  the CHECKER   E != 0 forces |Lambda| + |mu| even, so the survivors sit on alternate lattice
                planes -- and adding one box moves you off the surviving set, which is the content
                of the checkerboard lemma;
  the SIGN      blue for +1, orange for -1, and the alternation along the cone is the transfer
                structure: each step multiplies by one local factor.

Everything is computed from the operator; nothing is placed by hand. The pale points are the lattice
points of the cone that the parity kills, drawn so that the surviving set can be seen as a sublattice
rather than as a scatter.

Palette: the two validated categorical hues for the two signs, house neutrals for the killed points
and the cone. Marker shape doubles the colour so the figure survives greyscale printing.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"
PLUS = "#2a78d6"; MINUS = "#eb6834"; DEAD = "#d6d5cf"; CONE = "#b9b7ae"

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

INF = 10 ** 9
R = 2          # Lambda con 2 partes, mu con 1: el dominio es 3-dimensional


def entrelaza(mu, Lam):
    Lp = list(Lam) + [0]
    return all(Lp[i] >= mu[i] >= Lp[i + 2] for i in range(len(Lam) - 1))


def rs(mu, Lam):
    L = list(Lam) + [0]
    M = [INF] + list(mu) + [0] * (len(Lam) + 2)
    return [min(L[i], M[i]) - max(L[i + 1], M[i + 1]) for i in range(len(Lam))]


def E(mu, Lam):
    if not entrelaza(mu, Lam):
        return None                      # fuera del cono
    r = rs(mu, Lam)
    if any(x % 2 for x in r):
        return 0                         # dentro del cono, matado por la paridad
    return (-1) ** (sum(r) // 2)


W = 12
P, M, D = [], [], []
for L1 in range(W + 1):
    for L2 in range(L1 + 1):
        for u1 in range(W + 1):
            v = E((u1,), (L1, L2))
            if v is None:
                continue
            (P if v == 1 else M if v == -1 else D).append((L1, L2, u1))

fig = plt.figure(figsize=(9.6, 4.3))

for k, (elev, azim, titulo) in enumerate([(20, -60, r"the cone, and the surviving sublattice"),
                                          (8, -78, r"edge on: the surviving planes alternate")]):
    ax = fig.add_subplot(1, 2, k + 1, projection="3d")
    if D:
        X, Y, Z = zip(*D)
        ax.scatter(X, Y, Z, s=3.5, c=DEAD, marker=".", linewidths=0, depthshade=False, alpha=0.55)
    for L, col, mk, sz in [(P, PLUS, "o", 15), (M, MINUS, "s", 13)]:
        if L:
            X, Y, Z = zip(*L)
            ax.scatter(X, Y, Z, s=sz, c=col, marker=mk, linewidths=0, depthshade=True)
    # las aristas del cono: Lambda_1 = Lambda_2, y mu_1 = Lambda_1
    ax.plot([0, W], [0, W], [0, 0], color=CONE, linewidth=0.9)
    ax.plot([0, W], [0, 0], [0, W], color=CONE, linewidth=0.9)
    ax.plot([0, W], [0, W], [0, W], color=CONE, linewidth=0.9)
    ax.set_xlabel(r"$\Lambda_1$", labelpad=-4)
    ax.set_ylabel(r"$\Lambda_2$", labelpad=-4)
    ax.set_zlabel(r"$\mu_1$", labelpad=-4)
    ax.set_xlim(0, W); ax.set_ylim(0, W); ax.set_zlim(0, W)
    ax.set_xticks([0, 6, 12]); ax.set_yticks([0, 6, 12]); ax.set_zticks([0, 6, 12])
    ax.tick_params(labelsize=7, pad=-2)
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect((1, 1, 1))
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_facecolor(SURFACE); pane.pane.set_edgecolor(MUTED); pane.pane.set_alpha(1.0)
        pane._axinfo["grid"]["color"] = "#eeede8"; pane._axinfo["grid"]["linewidth"] = 0.4
    ax.set_title(titulo, fontsize=9, color=INK, pad=0)

fig.text(0.5, 0.055,
         r"$%d$ lattice points satisfy the interlacing; $%d$ survive the parity "
         r"($%d$ with $+1$, $%d$ with $-1$), and $%d$ are killed"
         % (len(P) + len(M) + len(D), len(P) + len(M), len(P), len(M), len(D)),
         ha="center", fontsize=8.5, color=SECOND)
fig.text(0.5, 0.008,
         r"the survivors lie on alternate planes $|\Lambda|+|\mu|$ even: adding a single box "
         r"always leaves the set",
         ha="center", fontsize=8.5, color=SECOND)

leyenda = [Line2D([], [], color=PLUS, marker="o", linestyle="none", markersize=5, label=r"$E=+1$"),
           Line2D([], [], color=MINUS, marker="s", linestyle="none", markersize=5, label=r"$E=-1$"),
           Line2D([], [], color=DEAD, marker=".", linestyle="none", markersize=6,
                  label=r"in the cone, killed by the parity"),
           Line2D([], [], color=CONE, linewidth=1.0, label=r"edges of the interlacing cone")]
fig.legend(handles=leyenda, frameon=False, fontsize=8, ncol=4,
           loc="lower center", bbox_to_anchor=(0.5, 0.10))
fig.subplots_adjust(left=0.01, right=0.99, top=1.0, bottom=0.20, wspace=0.03)

fig.savefig("fig_cone3d.pdf", bbox_inches="tight", pad_inches=0.02)
print("fig_cone3d.pdf escrito")
print("  en el cono: %d   sobreviven: %d (+1: %d, -1: %d)   matados por paridad: %d"
      % (len(P) + len(M) + len(D), len(P) + len(M), len(P), len(M), len(D)))
par = all((a + b + c) % 2 == 0 for a, b, c in P + M)
print("  control: |Lambda|+|mu| par en TODOS los supervivientes : %s" % ("si" if par else "*** NO ***"))
