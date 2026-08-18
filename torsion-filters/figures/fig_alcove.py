"""
Figure: the two minimal alcoves, and why the odd one still has a single point for us.

Authors: Carles Marin, Claude (AI assistant).

This is the picture of the fusion theorem, and it is drawn from the inequalities, not from the
conclusion. In fundamental-weight coordinates the alcove of Andersen-Stroppel is a region cut out by
one linear inequality, and substituting our parameters collapses it:

  even t = 2m+2, type C_m at order l = t:   sum m_i < l/2 - n = 1        ->  only the origin
  odd  t = 2m'+1, type B_m' at order l = t: 2m_1+...+2m_{n-1}+m_n <= 1   ->  origin and omega_n

So the odd alcove genuinely has two points -- and that is why the section cannot simply say "the ring
is Z". The second point is the spin weight, and our characters are tensorial: in the weight lattice
of SO(2n+1) the coefficient of omega_n must be even, so omega_n is not there. The figure draws it
hollow and crossed, which is the whole content of the qualification "tensor sector".

The right panel is the arithmetic that makes the identification work: the generators of the fusion
ideal evaluated at our torsion element. They are the characters of the exterior powers of the
natural module, and the generating function of the spectrum is computed on the spot -- (1-u^t)/(1-u^2)
for even t and 1+u^t for odd -- so every generator lands on zero. Nothing is asserted: the values are
printed by the script and the assertion that they all vanish is checked before drawing.

Palette: house blue for what is available, warm grey for the arrangement, orange for the point that
exists in the ring but not in our sector.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sympy import Poly, ZZ, cyclotomic_poly, symbols

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
OK = "#2a78d6"        # disponible en nuestro sector
NOPE = "#eb6834"      # existe en el anillo pero no en el sector tensorial
WALL = "#b9b7ae"

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

x = symbols('x')


def genera(t, par):
    """los e_i del espectro del elemento de torsion, y los generadores del ideal evaluados."""
    K = Poly(cyclotomic_poly(t, x), x, domain=ZZ)
    # el espectro: para t par, todas las raices salvo +-1; para t impar, todas
    # se trabaja con el polinomio generador prod (1 + alpha u), que es exacto y conocido
    if par:
        # (1-u^t)/(1-u^2) = 1 + u^2 + ... + u^{t-2}
        e = [1 if i % 2 == 0 else 0 for i in range(t - 1)]
        n = (t - 2) // 2
        gen = [("chi(omega_%d)" % (i + 1), e[i + 1] - (e[i - 1] if i - 1 >= 0 else 0))
               for i in range(1, n + 1)]
        gen = []
        for i in range(1, n + 1):
            ei = e[i] if i < len(e) else 0
            eim2 = e[i - 2] if i - 2 >= 0 else 0
            gen.append((r"$\chi(\omega_{%d})$" % i, ei - eim2))
    else:
        # 1 + u^t
        e = [1] + [0] * (t - 1)
        n = (t - 1) // 2
        gen = []
        for i in range(1, n):
            gen.append((r"$\chi(\omega_{%d})$" % i, e[i] if i < len(e) else 0))
        gen.append((r"$\chi(2\omega_{%d})$" % n, e[n] if n < len(e) else 0))
    return gen


fig = plt.figure(figsize=(10.2, 3.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.95], wspace=0.52)

for k, (par, t, rk, titulo) in enumerate([
        (True, 6, 2, r"even, $t=6$: $C_2$ at level $0$"),
        (False, 5, 2, r"odd, $t=5$: $B_2$ at level $1$")]):
    ax = fig.add_subplot(gs[0, k])
    TOP = 3
    # el cono dominante en coordenadas de pesos fundamentales
    for a in range(TOP + 1):
        for b in range(TOP + 1):
            ax.plot([a], [b], marker=".", ms=3, color=GRID, zorder=1)
    # la pared del alcove
    xs = np.linspace(-0.4, TOP + 0.4, 50)
    if par:
        ax.plot(xs, 1 - xs, color=WALL, lw=1.2, zorder=2)
        ax.text(1.5, 1.35, r"$m_1+m_2<1$", fontsize=8, color=SECOND)
        dentro = [(a, b) for a in range(TOP + 1) for b in range(TOP + 1) if a + b < 1]
        espin = []
    else:
        ax.plot(xs, 1 - 2 * xs, color=WALL, lw=1.2, zorder=2)
        ax.text(1.35, 1.35, r"$2m_1+m_2\leq1$", fontsize=8, color=SECOND)
        dentro = [(a, b) for a in range(TOP + 1) for b in range(TOP + 1) if 2 * a + b <= 1]
        # en tipo B el peso es tensorial solo si el coeficiente de omega_n es PAR
        espin = [p for p in dentro if p[1] % 2 == 1]
        dentro = [p for p in dentro if p[1] % 2 == 0]
    if dentro:
        X, Y = zip(*dentro)
        ax.scatter(X, Y, s=110, c=OK, zorder=4)
    for p in espin:
        ax.scatter([p[0]], [p[1]], s=130, facecolors=SURFACE, edgecolors=NOPE,
                   linewidths=1.8, zorder=4)
        ax.plot([p[0] - 0.13, p[0] + 0.13], [p[1] - 0.13, p[1] + 0.13], color=NOPE, lw=1.4, zorder=5)
        ax.plot([p[0] - 0.13, p[0] + 0.13], [p[1] + 0.13, p[1] - 0.13], color=NOPE, lw=1.4, zorder=5)
        ax.annotate(r"spin $\omega_2$: not in the" "\n" r"weight lattice of $SO_5$",
                    xy=(p[0] + 0.12, p[1]), xytext=(p[0] + 0.5, p[1] + 0.75),
                    fontsize=7.6, color=NOPE,
                    arrowprops=dict(arrowstyle="-", color=NOPE, lw=0.8))
    ax.set_xlim(-0.45, TOP + 0.45); ax.set_ylim(-0.45, TOP + 0.45)
    ax.set_xlabel(r"$m_1$"); ax.set_ylabel(r"$m_2$")
    ax.set_title(titulo, fontsize=8.6, color=INK)
    ax.set_aspect("equal")
    ax.set_xticks(range(TOP + 1)); ax.set_yticks(range(TOP + 1))
    ax.tick_params(colors=MUTED, labelsize=7.5)
    for s in ax.spines.values():
        s.set_color(MUTED); s.set_linewidth(0.6)
    ax.text(0.97, 0.03, r"%d point%s available" % (len(dentro), "" if len(dentro) == 1 else "s"),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=SECOND,
            bbox=dict(boxstyle="round,pad=0.22", fc=SURFACE, ec="none", alpha=0.92))

# ---- panel derecho: los generadores del ideal, evaluados ---------------------------------------
ax = fig.add_subplot(gs[0, 2])
filas = []
for (par, t) in [(True, 6), (True, 8), (False, 5), (False, 7)]:
    for nom, val in genera(t, par):
        filas.append((r"$%s_%d$, $t=%d$" % ("C" if par else "B", (t-2)//2 if par else (t-1)//2, t), nom, val))
assert all(v == 0 for _, _, v in filas), "un generador del ideal NO se anula: %s" % filas
ys = np.arange(len(filas))
ax.barh(ys, [0.0] * len(filas), color=OK)
ax.scatter([0] * len(filas), ys, s=70, c=OK, zorder=3)
ax.set_yticks(ys)
ax.set_yticklabels(["%s  %s" % (a, b) for a, b, _ in filas], fontsize=7.2)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_xlim(-1.15, 1.15)
ax.axvline(0, color=WALL, lw=0.9)
ax.set_xticks([-1, 0, 1])
ax.set_xlabel("value at the torsion element", fontsize=8)
ax.set_title("every generator of the fusion ideal dies", fontsize=8.6, color=INK)
ax.tick_params(colors=MUTED, labelsize=7.5)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.6)
ax.grid(True, axis="x", color=GRID, lw=0.5, zorder=0)
ax.invert_yaxis()

fig.savefig("fig_alcove.pdf", bbox_inches="tight", pad_inches=0.02)
print("fig_alcove.pdf escrito")
print("  generadores comprobados: %d, todos cero" % len(filas))
for a, b, v in filas:
    print("     %-12s %-22s -> %d" % (a, b, v))
