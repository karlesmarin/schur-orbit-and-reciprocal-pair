"""
Figure: the zero locus of Theorem 8.1, counted, for one, two and three reciprocal pairs.

Authors: Carles Marin, Claude (AI assistant).

A first attempt plotted each vanishing partition at (|lambda|, ell(lambda)). That projection is not
injective -- hundreds of shapes collapse onto a handful of marks -- so it showed presence and hid
multiplicity. Here the vertical axis is the honest one: for each size |lambda| we count how many
partitions with at most N parts are killed, split by the branch of the theorem that kills them.
Branch (a) is the constant-parity beta set, branch (b) is self-complementarity of odd width. The
faint step line is the total number of shapes of that size, on the same axis, so that the rarity of
the locus is visible rather than asserted.

What the figure is for: the zeros are not spread out. They sit at sizes in arithmetic progression,
the two branches occupy different sizes, and the whole locus is a thin set inside a fast-growing
background. Everything is recomputed here from the definitions.

Palette: the two validated categorical hues (blue = branch (a), orange = branch (b)) plus the
neutral for the background count. Hatching doubles the colour so the figure survives greyscale.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
BRA = "#2a78d6"      # blue   -- branch (a)
BRB = "#eb6834"      # orange -- branch (b)
NONE = "#c9c8c1"     # neutral -- the background count

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})


def partitions(n, maxlen, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, maxlen - 1, p):
            yield (p,) + rest


def classify(lam, N):
    L = list(lam) + [0] * (N - len(lam))
    beta = [L[i] + (N - 1 - i) for i in range(N)]
    if len({b % 2 for b in beta}) == 1:
        return 'a'
    w = L[0] + L[N - 1]
    if w % 2 == 1 and all(L[i] + L[N - 1 - i] == w for i in range(N)):
        return 'b'
    return None


PANELS = [1, 2, 3]
SMAX = {1: 30, 2: 30, 3: 30}

fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.25), sharex=True)
for ax, r in zip(axes, PANELS):
    N = 2 * r + 2
    xs = list(range(0, SMAX[r] + 1))
    ca, cb, tot = [], [], []
    for S in xs:
        a = b = t = 0
        for lam in partitions(S, N):
            t += 1
            k = classify(lam, N)
            if k == 'a':
                a += 1
            elif k == 'b':
                b += 1
        ca.append(a); cb.append(b); tot.append(t)
    ax.step(xs, tot, where='mid', color=NONE, lw=1.1, zorder=1)
    ax.fill_between(xs, tot, step='mid', color=NONE, alpha=0.5, lw=0, zorder=0)
    ax.bar(xs, cb, width=0.72, color=BRB, zorder=3, label="branch (b)")
    ax.bar(xs, ca, width=0.72, bottom=cb, color='none', edgecolor=BRA, hatch='///',
           linewidth=0.8, zorder=4, label="branch (a)")
    ax.set_yscale('symlog', linthresh=1)
    ax.set_title(r"$r=%d$,  $N=%d$" % (r, N), fontsize=9, color=INK, pad=5)
    ax.set_xlabel(r"$|\lambda|$", color=INK)
    if r == 1:
        ax.set_ylabel("number of partitions", color=INK)
    ax.grid(True, axis='y', color=GRID, lw=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_color(MUTED); sp.set_linewidth(0.7)
    ax.tick_params(colors=SECOND, labelsize=7.6)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=7))
    nz = sum(ca) + sum(cb)
    ax.text(0.03, 0.955, r"$%d$ zeros of $%d$ shapes" % (nz, sum(tot)), transform=ax.transAxes,
            ha='left', va='top', fontsize=7.4, color=SECOND)
    ax.text(0.03, 0.875, r"branch (b) only at $|\lambda|=w(r{+}1)$, $w$ odd", transform=ax.transAxes,
            ha='left', va='top', fontsize=7.0, color=BRB)

h = [plt.Rectangle((0, 0), 1, 1, facecolor=BRB, edgecolor='none',
                   label=r"branch (b): self-complementary, odd width"),
     plt.Rectangle((0, 0), 1, 1, facecolor='none', edgecolor=BRA, hatch='///', linewidth=0.8,
                   label=r"branch (a): $\beta$ of constant parity"),
     plt.Rectangle((0, 0), 1, 1, facecolor=NONE, alpha=0.6, edgecolor=NONE,
                   label=r"all shapes with at most $N$ parts")]
fig.legend(handles=h, loc='lower center', ncol=3, frameon=False, fontsize=7.8,
           bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.055, 1, 1))
fig.savefig("fig_zeros.pdf"); fig.savefig("fig_zeros.png", dpi=170)
print("wrote fig_zeros.pdf / .png")
