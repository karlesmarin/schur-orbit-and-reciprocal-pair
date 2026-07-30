"""
Figure: why self-complementarity kills the character, in one picture.

Authors: Carles Marin, Claude (AI assistant).

The bialternant is expanded along the two frozen rows coming from the letters 1 and -1. Their 2x2
minor is zero unless the two chosen columns carry beta's of OPPOSITE parity, so the surviving terms
are indexed by pairs (e,o) with e in E and o in O. The reflection x -> c - x, with c = w + N - 1,
acts on those terms; when beta = c - beta it maps the term set to itself, and because c is even it
preserves the parity classes, which is exactly what makes the parity sign flip while the alternant
does not. The involution is free -- a fixed point would need two beta's of equal parity in an
opposite-parity pair -- so the terms cancel in pairs and the sum is zero.

The left panel is a self-complementary shape of odd width: every arc lands on another element of
beta, and the arcs join same-parity entries. The right panel is a shape one box away from it: the
arcs miss beta, there is no pairing, and nothing cancels. Both beta sets and both arc systems are
computed from the definitions.

Palette: blue for even entries of beta, orange for odd, neutral for the integers not in beta.
Filled versus open markers double the colour for greyscale.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
EVEN = "#2a78d6"     # blue   -- even beta
ODD = "#eb6834"      # orange -- odd beta
NONE = "#cfcec7"     # neutral

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

N = 6                     # r = 2
RHO = list(range(N - 1, -1, -1))


def beta_of(lam):
    L = list(lam) + [0] * (N - len(lam))
    return [L[i] + RHO[i] for i in range(N)], L


def selfcomp(L):
    w = L[0] + L[N - 1]
    return (w if all(L[i] + L[N - 1 - i] == w for i in range(N)) else None)


CASES = [((3, 3, 2, 1, 0, 0), "self-complementary, $w=3$ odd: the character vanishes"),
         ((3, 3, 2, 2, 0, 0), "one box away: no pairing, and it does not")]

fig, axes = plt.subplots(2, 1, figsize=(9.6, 4.5))
for ax, (lam, title) in zip(axes, CASES):
    beta, L = beta_of(lam)
    w = selfcomp(L)
    c = (beta[0] + beta[-1])
    hi = max(beta) + 1
    ax.plot([-0.6, hi + 0.6], [0, 0], color=MUTED, lw=0.7, zorder=1)
    for x in range(0, hi + 1):
        if x not in beta:
            ax.plot([x], [0], marker='.', color=NONE, ms=4.5, zorder=2)
    for b in beta:
        col = EVEN if b % 2 == 0 else ODD
        ax.plot([b], [0], marker='o', color=col, ms=7.5, zorder=4)
        ax.annotate(str(b), (b, 0), textcoords="offset points", xytext=(0, -14),
                    ha='center', fontsize=7.4, color=col)
    # reflection arcs b -> c-b.  Collect the pairs first and de-duplicate, so that an arc whose
    # image falls OUTSIDE beta is still drawn -- drawing only when c-b > b silently hid exactly the
    # arcs that fail, which is the case the figure exists to show.
    ok = all((c - b) in beta for b in beta)
    pairs = sorted({(min(b, c - b), max(b, c - b)) for b in beta if c - b != b})
    fixed = [b for b in beta if c - b == b]
    for lo, hi_ in pairs:
        lands = (lo in beta) and (hi_ in beta)
        mid = (lo + hi_) / 2.0
        span = hi_ - lo
        ax.add_patch(Arc((mid, 0), span, span * 0.85, theta1=0, theta2=180,
                         color=(SECOND if lands else MUTED), lw=1.0,
                         linestyle='-' if lands else (0, (2, 2)), zorder=3))
        for end in (lo, hi_):
            if end not in beta:
                ax.plot([end], [0], marker='x', color=MUTED, ms=6.5, zorder=5)
                ax.annotate(r"$\notin\beta$", (end, 0), textcoords="offset points",
                            xytext=(0, -14), ha='center', fontsize=7.0, color=MUTED)
    for b in fixed:
        ax.plot([b], [0], marker='D', markerfacecolor='none', markeredgecolor=SECOND,
                ms=9, mew=0.9, zorder=5)
    ax.set_title(r"$\lambda=(%s)$,  $\beta=(%s)$,  $c=%d$ --- %s"
                 % (",".join(str(t) for t in lam if t > 0),
                    ",".join(str(b) for b in beta), c, title),
                 fontsize=8.5, color=INK, pad=7)
    ax.set_xlim(-1.0, hi + 1.0); ax.set_ylim(-0.9, (max(beta) - min(beta)) * 0.48 + 0.5)
    ax.axis('off')
    ax.text(hi + 0.4, -0.55, r"$x\mapsto c-x$" + ("  closes on $\\beta$" if ok else "  misses $\\beta$"),
            ha='right', va='center', fontsize=7.6, color=SECOND if ok else MUTED)

h = [plt.Line2D([], [], marker='o', color=EVEN, lw=0, ms=6, label=r"$\beta_j$ even"),
     plt.Line2D([], [], marker='o', color=ODD, lw=0, ms=6, label=r"$\beta_j$ odd"),
     plt.Line2D([], [], color=SECOND, lw=1.0, label=r"$x\mapsto c-x$, landing in $\beta$"),
     plt.Line2D([], [], color=MUTED, lw=1.0, ls=(0, (2, 2)), label=r"landing outside $\beta$")]
fig.legend(handles=h, loc='lower center', ncol=4, frameon=False, fontsize=7.8,
           bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.055, 1, 1))
fig.savefig("fig_involution.pdf"); fig.savefig("fig_involution.png", dpi=170)
print("wrote fig_involution.pdf / .png")
for lam, _ in CASES:
    beta, L = beta_of(lam)
    print("  lam=%s beta=%s selfcomp width=%s" % (lam, beta, selfcomp(L)))
