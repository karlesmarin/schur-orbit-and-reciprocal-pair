"""
Figure: the unit conjecture, drawn as what it looks like from inside -- a wandering sum that lands.

Authors: Carles Marin, Claude (AI assistant).

The conjecture says the coefficient at the highest surviving weight is always a unit. Stated as a
number that is a table entry it is easy to underrate, because the reader has no way to see that it
is surprising. What makes it surprising is the size of the terms being added: the branching
multiplicities at that weight are not small, and the filter only supplies signs.

So the figure adds them one at a time. Each path is one form: the running partial sum of
B_{eta,mu} * tau_t(eta) over the surviving eta, taken in decreasing order of eta. The paths wander
-- tens away from zero -- and every one of them ends on +1 or -1. Nothing in the two factors
separately forces that; the branching does not know about the root of unity, and the filter does not
know about the branching.

LEFT: the paths, with the band |A| = 1 marked, and the endpoints on it.
RIGHT: the same fact as two distributions. How far the partial sums travel, against where they end.
The left distribution is spread over an order of magnitude; the right one is a single value.

Everything is computed; the paths are the actual summands in the actual order.

Palette: the two validated categorical hues by the sign of the landing (blue = +1, orange = -1),
house neutrals for the band and the axes. Endpoint marker shape doubles the colour so the figure
survives greyscale printing, as in the companion paper.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
PLUS = "#2a78d6"; MINUS = "#eb6834"; BAND = "#e6e5de"

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

D = json.load(open("wall_table_DUMP.json"))
FILAS = D["filas"]

CAMINOS = []
for f in FILAS:
    viv = [x for x in f["bloque_mu_max"] if x["tau"]]          # ya vienen en orden de eta decreciente
    if not viv:
        continue
    term = [x["B"] * x["tau"] for x in viv]
    parc = np.cumsum(term)
    CAMINOS.append({"beta": f["beta"], "term": term, "parc": parc,
                    "A": int(parc[-1]), "maxB": max(abs(t) for t in term),
                    "maxparc": int(max(abs(v) for v in parc))})

assert CAMINOS, "el volcado no trae bloques de mu_max"
for c in CAMINOS:
    assert abs(c["A"]) == 1, "una forma no acaba en +-1: %s -> %d" % (c["beta"], c["A"])

fig = plt.figure(figsize=(9.6, 3.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1.0], wspace=0.26)

# ------------------------------------------------------------------ IZQUIERDA: las caminatas ----
ax = fig.add_subplot(gs[0, 0])
LMAX = max(len(c["parc"]) for c in CAMINOS)
ax.axhspan(-1, 1, color=BAND, zorder=1)
ax.axhline(0, color=MUTED, linewidth=0.6, zorder=2)
for c in CAMINOS:
    col = PLUS if c["A"] == 1 else MINUS
    y = np.concatenate([[0], c["parc"]])
    x = np.arange(len(y))
    ax.plot(x, y, color=col, linewidth=1.0, alpha=0.72, zorder=3, solid_capstyle="round")
    ax.plot([x[-1]], [y[-1]], marker="o" if c["A"] == 1 else "s", color=col, markersize=4.6,
            linewidth=0, zorder=5)
ax.set_xlabel(r"terms added, in decreasing order of $\eta$")
ax.set_ylabel(r"partial sum of $B_{\eta,\mu_{\max}}\,\tau_t(\eta)$")
ax.set_title(r"every path wanders, every path lands on $\pm1$", fontsize=9, color=INK)
ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_xlim(0, LMAX)
lim = max(c["maxparc"] for c in CAMINOS) * 1.12
ax.set_ylim(-lim, lim)
ax.text(0.985, 0.965,
        "the band $|A_{\\mu_{\\max}}|\\leq1$ is the thin strip at zero;\n"
        "the inset is the same strip, magnified", transform=ax.transAxes,
        fontsize=8, color=SECOND, ha="right", va="top", linespacing=1.4)

# el recuadro: la banda a escala propia, que es donde esta el enunciado
axi = ax.inset_axes([0.065, 0.07, 0.27, 0.30])
axi.axhspan(-1, 1, color=BAND, zorder=1)
axi.axhline(0, color=MUTED, linewidth=0.5, zorder=2)
# alineados por PASOS DESDE EL FINAL, no por el indice absoluto: los caminos miden 6 y 9, y con el
# indice absoluto los finales caian en abscisas distintas y el recuadro salia un enrejado.
for c in CAMINOS:
    col = PLUS if c["A"] == 1 else MINUS
    y = np.concatenate([[0], c["parc"]])
    axi.plot([-2, -1, 0], y[-3:], color=col, linewidth=0.9, alpha=0.75, zorder=3)
    axi.plot([0], [y[-1]], marker="o" if c["A"] == 1 else "s", color=col,
             markersize=4.6, linewidth=0, zorder=5)
axi.set_ylim(-3.2, 3.2); axi.set_yticks([-1, 0, 1])
axi.set_xlim(-2.1, 0.45)
axi.set_xticks([-2, -1, 0]); axi.set_xticklabels([r"$-2$", r"$-1$", "end"], fontsize=6.5)
axi.tick_params(labelsize=7)
axi.set_facecolor(SURFACE)
for sp in axi.spines.values():
    sp.set_color(MUTED); sp.set_linewidth(0.6)
axi.set_title(r"the landing", fontsize=7.5, color=SECOND, pad=2)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.6)

# ------------------------------------------------------------------ DERECHA: las dos escalas ----
ax = fig.add_subplot(gs[0, 1])
viaje = [c["maxparc"] for c in CAMINOS]
final = [abs(c["A"]) for c in CAMINOS]
rng = np.random.default_rng(20260815)
jx = rng.uniform(-0.13, 0.13, len(CAMINOS))
ax.scatter(0 + jx, viaje, s=26, c=[PLUS if c["A"] == 1 else MINUS for c in CAMINOS],
           marker="o", linewidths=0, alpha=0.85, zorder=3)
ax.scatter(1 + jx, final, s=26, c=[PLUS if c["A"] == 1 else MINUS for c in CAMINOS],
           marker="D", linewidths=0, alpha=0.85, zorder=3)
for j, c in zip(jx, CAMINOS):
    ax.plot([0 + j, 1 + j], [c["maxparc"], abs(c["A"])], color=MUTED, linewidth=0.55,
            alpha=0.5, zorder=2)
ax.set_xticks([0, 1])
ax.set_xticklabels([r"$\max$ partial sum", r"$|A_{\mu_{\max}}|$"], fontsize=8.5)
ax.set_xlim(-0.45, 1.45)
ax.set_yscale("log"); ax.set_ylim(0.55, 1500)
ax.set_ylabel(r"absolute value (log)")
ax.set_title(r"an order of magnitude, then one value", fontsize=9, color=INK)
ax.grid(True, axis="y", color=GRID, linewidth=0.5, zorder=0)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.6)
ax.text(0.5, 0.03, r"$n=%d$ forms;  largest single term $|B\tau|=%d$"
        % (len(CAMINOS), max(c["maxB"] for c in CAMINOS)),
        transform=ax.transAxes, fontsize=8, color=SECOND, ha="center")

leyenda = [Line2D([], [], color=PLUS, marker="o", linewidth=1.2, markersize=4.6,
                  label=r"lands on $A_{\mu_{\max}}=+1$"),
           Line2D([], [], color=MINUS, marker="s", linewidth=1.2, markersize=4.6,
                  label=r"lands on $A_{\mu_{\max}}=-1$")]
fig.legend(handles=leyenda, frameon=False, fontsize=8, ncol=2,
           loc="lower center", bbox_to_anchor=(0.5, -0.075))

fig.savefig("fig_collapse.pdf", bbox_inches="tight", pad_inches=0.02)
print("fig_collapse.pdf escrito")
print("  %d caminos, longitudes %s" % (len(CAMINOS), sorted(len(c['parc']) for c in CAMINOS)))
print("  termino individual mas grande |B tau| = %d" % max(c["maxB"] for c in CAMINOS))
print("  suma parcial mas lejana        = %d" % max(c["maxparc"] for c in CAMINOS))
print("  finales: %s" % sorted(c["A"] for c in CAMINOS))
