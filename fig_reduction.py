"""
Figure: the type-D reduction on this alphabet, and why the criterion is a finite check.

Authors: Carles Marin, Claude (AI assistant).

The expansion s_lambda = sum c_nu o_nu is an identity of symmetric functions and holds in every
range; specializing it at the alphabet moves the whole range-dependence into the values o_nu(A).
Those values are computed here, label by label, and they fall into four classes:

  standard with nu'_1 < N/2   -- the surviving BASIS, and it is linearly independent
  standard with nu'_1 = N/2   -- self-associate: the value is zero
  standard with nu'_1 > N/2   -- folds onto the basis as +- o_{nu*}
  non-standard                -- either zero, or +-1 times a single basis element

No general modification rule is used: King's rules are the classical source and Fauser, Jarvis, King
and Wybourne record that no general formalism exists, so the reduction is obtained rather than
recalled. What the figure shows is that everything lands on a small independent basis, which is why
the vanishing becomes the finite condition C_mu = 0 for all mu.

Data: fig_data_new.json, computed in exact arithmetic. Bar heights are counts of labels, and the
split between vanishing and folding inside each class is measured, not assumed.

Palette: blue for what survives to the basis, orange for what folds onto it, neutral for what dies.
Hatching doubles the colour for greyscale.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
KEEP = "#2a78d6"     # blue   -- the basis
FOLD = "#eb6834"     # orange -- folds onto the basis
DIE = "#cfcec7"      # neutral -- vanishes

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

D = json.load(open("fig_data_new.json"))["reduction"]

ORDER = ["basis", "selfassoc", "assoc", "nonstd"]
NAMES = {"basis": "standard,\n$\\nu'_1<N/2$", "selfassoc": "standard,\n$\\nu'_1=N/2$",
         "assoc": "standard,\n$\\nu'_1>N/2$", "nonstd": "non-standard,\n$\\nu'_1+\\nu'_2>N$"}

fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.5))
for ax, key in zip(axes, ["r1", "r2"]):
    blk = D[key]
    N = blk["N"]; half = blk["half"]; tab = blk["table"]
    zero, live = [], []
    for cls in ORDER:
        rows = [t for t in tab if t["cls"] == cls]
        zero.append(sum(1 for t in rows if t["zero"]))
        live.append(sum(1 for t in rows if not t["zero"]))
    xs = range(len(ORDER))
    ax.bar(xs, zero, width=0.62, color=DIE, edgecolor=MUTED, linewidth=0.6, zorder=3,
           label="value is zero")
    cols = [KEEP, KEEP, FOLD, FOLD]
    ax.bar(xs, live, width=0.62, bottom=zero, color='none', edgecolor=None, zorder=3)
    for i, (c, v, z) in enumerate(zip(cols, live, zero)):
        if v:
            ax.bar([i], [v], width=0.62, bottom=[z], color=c, zorder=4)
    top = max(z + v for z, v in zip(zero, live))
    for i, (z, v) in enumerate(zip(zero, live)):
        if z:
            ax.text(i, z / 2, str(z), ha='center', va='center', fontsize=7.6, color=SECOND)
        if v:
            ax.text(i, z + v + top * 0.022, str(v), ha='center', va='bottom', fontsize=7.6,
                    color=(KEEP if i < 2 else FOLD))
    ax.set_ylim(0, top * 1.13)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([NAMES[c] for c in ORDER], fontsize=7.3, color=INK)
    ax.set_title(r"$r=%d$,  $N=%d$:  %d labels, basis of %d"
                 % ((N - 2) // 2, N, len(tab),
                    sum(1 for t in tab if t["cls"] == "basis" and not t["zero"])),
                 fontsize=8.8, color=INK, pad=6)
    if key == "r1":
        ax.set_ylabel("number of labels $\\nu$", color=INK)
    ax.grid(True, axis='y', color=GRID, lw=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_color(MUTED); sp.set_linewidth(0.7)
    ax.tick_params(colors=SECOND, labelsize=7.6)

h = [plt.Rectangle((0, 0), 1, 1, facecolor=KEEP, edgecolor='none',
                   label="survives: an element of the basis"),
     plt.Rectangle((0, 0), 1, 1, facecolor=FOLD, edgecolor='none',
                   label=r"folds onto the basis, as $\pm\,o_\mu(A)$"),
     plt.Rectangle((0, 0), 1, 1, facecolor=DIE, edgecolor=MUTED, linewidth=0.6,
                   label=r"$o_\nu(A)=0$")]
fig.legend(handles=h, loc='lower center', ncol=3, frameon=False, fontsize=7.8,
           bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig("fig_reduction.pdf"); fig.savefig("fig_reduction.png", dpi=170)
print("wrote fig_reduction.pdf / .png")
for key in ["r1", "r2"]:
    tab = D[key]["table"]
    for cls in ORDER:
        rows = [t for t in tab if t["cls"] == cls]
        print("  %s %-10s total %3d  zero %3d  live %3d"
              % (key, cls, len(rows), sum(1 for t in rows if t["zero"]),
                 sum(1 for t in rows if not t["zero"])))
