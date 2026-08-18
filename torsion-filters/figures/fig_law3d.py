"""
Figure: the law for the highest weight, as three solids and one vector subtraction.

Authors: Carles Marin, Claude (AI assistant).

This is the figure that earns three dimensions, because its objects are solids and not points in a
box. Ostrowski's additivity says the Newton polytope of a product is the Minkowski sum of the
factors, and Phi = N_beta / N_delta, so

    Newt(N_beta)  =  Newt(Phi)  (+)  Newt(N_delta) .

All three are orbit polytopes of the Weyl group of the free symplectic factor, so each is determined
by its dominant vertex, and the identity descends to a subtraction of vertices -- which is the law.

The three panels are the three solids, at r = 3, drawn from the ACTUAL supports: the convex hull of
the exponents that really occur, not the hull the law would predict. Drawing conv(W(C_r) mu_max)
instead would be drawing the statement rather than the data. The dominant vertex of each is marked,
and the arithmetic under the panels is the law: the vertex of the middle solid is the vertex of the
right one minus the vertex of the left one, and that difference is the sum of the positive roots of
the free factor plus a uniform offset t-1.

Palette: house neutrals for the solids, the two validated categorical hues for the vertices that
enter the subtraction. Wireframes rather than opaque faces so that the interior lattice points stay
visible in print.
"""
import json
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"
BLUE = "#2a78d6"; ORANGE = "#eb6834"; FACE = "#e6e5de"; EDGE = "#b9b7ae"

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

D = json.load(open("polytopes_DUMP.json"))
r, t, N = D["r"], D["t"], D["N"]

PANELES = [("Phi", r"$\mathrm{Newt}(\Phi_{t,r})$", D["top_Phi"], BLUE),
           ("Ndelta", r"$\mathrm{Newt}(N_\delta)$", D["top_Ndelta"], ORANGE),
           ("Nbeta", r"$\mathrm{Newt}(N_\beta)$", D["top_Nbeta"], INK)]

fig = plt.figure(figsize=(9.6, 3.9))
LIM = max(max(abs(v) for v in p) for p in D["Nbeta"]) + 1

CONTROLES = []
for k, (clave, titulo, top, colv) in enumerate(PANELES):
    P = np.array(D[clave], dtype=float)
    ax = fig.add_subplot(1, 3, k + 1, projection="3d")
    try:
        h = ConvexHull(P)
        caras = [P[s] for s in h.simplices]
        col = Poly3DCollection(caras, facecolors=FACE, edgecolors=EDGE, linewidths=0.35, alpha=0.42)
        ax.add_collection3d(col)
        nv = len(h.vertices)
        VERTS = {tuple(int(round(x)) for x in P[i]) for i in h.vertices}
    except Exception:
        nv = 0
        VERTS = set()
    ax.scatter(P[:, 0], P[:, 1], P[:, 2], s=1.6, c=MUTED, marker=".", linewidths=0,
               depthshade=False, alpha=0.45)
    # el vertice dominante y toda su orbita bajo los cambios de signo y las permutaciones
    orb = {tuple(s * v for s, v in zip(sg, pm))
           for pm in itertools.permutations(top) for sg in itertools.product([1, -1], repeat=r)}
    # C1  el pie afirma que cada solido es un politopo de orbita.  Tener el CARDINAL correcto no lo
    # demuestra: hay que comprobar que el conjunto de vertices ES la orbita W(C_r) del dominante.
    CONTROLES.append((clave, VERTS == orb, len(VERTS), len(orb)))
    O = np.array(sorted(orb), dtype=float)
    ax.scatter(O[:, 0], O[:, 1], O[:, 2], s=9, c=colv, marker="o", linewidths=0, depthshade=False,
               alpha=0.55)
    ax.scatter([top[0]], [top[1]], [top[2]], s=42, c=colv, marker="o", linewidths=0,
               depthshade=False, zorder=10)
    ax.text2D(0.5, 0.995, r"dominant vertex $%s$" % str(tuple(top)).replace(" ", ""),
              transform=ax.transAxes, fontsize=8.5, color=colv, ha="center", va="top")

    ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM); ax.set_zlim(-LIM, LIM)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=19, azim=-58)
    ax.set_xticks([-10, 0, 10]); ax.set_yticks([-10, 0, 10]); ax.set_zticks([-10, 0, 10])
    ax.tick_params(labelsize=7, pad=-2)
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_facecolor(SURFACE); pane.pane.set_edgecolor(MUTED); pane.pane.set_alpha(1.0)
        pane._axinfo["grid"]["color"] = "#eeede8"; pane._axinfo["grid"]["linewidth"] = 0.4
    ax.set_title(titulo, fontsize=9.5, color=INK, pad=8)
    ax.text2D(0.5, -0.02, r"%d exponents, %d vertices" % (len(P), nv), transform=ax.transAxes,
              fontsize=7.5, color=SECOND, ha="center")

sig = D["sigma"]
top_phi, top_del, top_bet = D["top_Phi"], D["top_Ndelta"], D["top_Nbeta"]
fig.text(0.5, 0.055,
         r"$\mathrm{Newt}(N_\beta)=\mathrm{Newt}(\Phi_{t,r})\oplus\mathrm{Newt}(N_\delta)$"
         r"$\qquad\Longrightarrow\qquad$"
         r"$%s=%s-%s$" % (str(tuple(top_phi)).replace(" ", ""),
                          str(tuple(top_bet)).replace(" ", ""),
                          str(tuple(top_del)).replace(" ", "")),
         ha="center", fontsize=10, color=INK)
fig.text(0.5, 0.005,
         r"and the subtracted vertex is $\sigma_r=%s=2\rho_{C_r}+(t-1)$, "
         r"the dominant vertex of the denominator, read in exponent coordinates"
         % str(tuple(sig)).replace(" ", ""),
         ha="center", fontsize=8.5, color=SECOND)
fig.subplots_adjust(left=0.005, right=0.995, top=1.0, bottom=0.155, wspace=0.0)

fig.savefig("fig_law3d.pdf", bbox_inches="tight", pad_inches=0.02)
print("fig_law3d.pdf escrito")
print("  t=%d r=%d beta=%s" % (t, r, D["beta"]))
print("  C1  Vert(P) == W(C_r)-orbita del dominante:")
for clave, ok, nv_, no_ in CONTROLES:
    print("      %-8s %s   (%d vertices, orbita %d)" % (clave, "SI" if ok else "NO", nv_, no_))
print("  C1 global: %d de %d" % (sum(1 for _, ok, _, _ in CONTROLES if ok), len(CONTROLES)))
for clave, _, top, _ in PANELES:
    print("  %-8s %5d exponentes, top %s" % (clave, len(D[clave]), top))
print("  ley: %s = %s - %s" % (top_phi, top_bet, top_del))
