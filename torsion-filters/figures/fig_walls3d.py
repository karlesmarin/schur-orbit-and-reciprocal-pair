"""
Figure: the walls of the torsion filter in type C_3 at t=8, and the ones the dual group cannot see.

Authors: Carles Marin, Claude (AI assistant).

Nothing here is drawn from the statement. Every point is a weight eta of C_3, placed at its shifted
coordinate a = eta + rho inside the cube of residues modulo t, and its colour is decided by two
computations run on the spot:

  * whether the filter survives, from the Weyl numerator det(xi^{i a_j} - xi^{-i a_j}) evaluated
    exactly in the cyclotomic field -- no character library, no lookup;
  * whether the point is regular in the group itself (roots e_i - e_j, e_i + e_j and the LONG root
    2 e_i) and whether it is regular in the dual (where the long root has become the short e_i).

The picture is the paper's second parity switch made visible. The grey points die on a wall both
readings see. The orange points die on the wall 2a = 0 that only the original group sees, and they
lie -- with no help from us -- on the three planes a_i = t/2: those are the weights at which the
centraliser criterion available in the literature gives the wrong answer.

A note on which population this is, because it is not the one quoted in the verification table. The
table sweeps partitions eta up to size 3t and counts 102 discrepancies at t = 8; this figure sweeps
the shifted coordinate a directly over strictly decreasing triples with max a <= 2t, which is 560
points and 96 discrepancies. The two are different windows on the same set, and the figure prints
its own counts when it runs, so the caption can quote them and not the table's.

The three planes are drawn faintly AFTER the points, so that if a single orange point fell off one
of them it would be visible. None does.

Palette: the house blue for what survives, the house orange for the discrepancy, warm grey for the
rest; the same three used throughout both papers.
"""
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from sympy import Poly, ZZ, cyclotomic_poly, symbols

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; MUTED = "#898781"
SURV = "#2a78d6"; DISC = "#eb6834"; DEAD = "#c9c6bf"

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

t, m = 8, 3
x = symbols('x')
phi = Poly(cyclotomic_poly(t, x), x, domain=ZZ)
_p = {}


def xp(k):
    k %= t
    if k not in _p:
        _p[k] = Poly(x ** k, x, domain=ZZ).rem(phi)
    return _p[k]


def det_poly(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    tot = Poly(0, x, domain=ZZ)
    for j in range(n):
        if M[0][j].is_zero:
            continue
        sub = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        term = (M[0][j] * det_poly(sub)).rem(phi)
        tot = tot + term if j % 2 == 0 else tot - term
    return tot.rem(phi)


def vive(a):
    M = [[(xp((i + 1) * a[j]) - xp(-(i + 1) * a[j])).rem(phi) for j in range(m)]
         for i in range(m)]
    return not det_poly(M).is_zero


# los pesos: a = eta + rho con eta particion de a lo sumo 3 partes, rho = (3,2,1),
# recorridos hasta que las clases modulo t se agotan
PTS = []
for a1 in range(1, 3 * t):
    for a2 in range(1, a1):
        for a3 in range(1, a2):
            a = (a1, a2, a3)
            if max(a) > 2 * t:
                continue
            c = tuple(v % t for v in a)
            v = vive(a)
            pares = all((a[i] - a[j]) % t != 0 and (a[i] + a[j]) % t != 0
                        for i, j in itertools.combinations(range(m), 2))
            largas = all((2 * a[i]) % t != 0 for i in range(m))
            cortas = all(a[i] % t != 0 for i in range(m))
            reg_G = pares and largas
            reg_D = pares and cortas
            PTS.append((c, v, reg_G, reg_D))

surv = [p[0] for p in PTS if p[1]]
disc = [p[0] for p in PTS if (not p[1]) and p[3]]      # muerto, pero el dual lo da por regular
dead = [p[0] for p in PTS if (not p[1]) and not p[3]]

# control impreso: los discrepantes tienen que estar TODOS en algun plano a_i = t/2
fuera = [c for c in disc if all(v != t // 2 for v in c)]
coinc = sum(1 for p in PTS if p[1] == p[2])
print("pesos dibujados            : %d" % len(PTS))
print("sobreviven                 : %d" % len(surv))
print("mueren y el dual no lo ve  : %d   (los del texto)" % len(disc))
print("mueren y el dual si lo ve  : %d" % len(dead))
print("(R) en G acierta           : %d de %d" % (coinc, len(PTS)))
print("discrepantes FUERA de a_i=t/2 : %d   (tiene que ser 0)" % len(fuera))
assert not fuera, "un punto naranja no esta en el plano a_i = t/2"
assert coinc == len(PTS), "(R) falla en algun peso dibujado"

fig = plt.figure(figsize=(7.4, 5.4))
ax = fig.add_subplot(111, projection="3d")


def sc(P, col, s, mark, lab, alpha=1.0, ec=None):
    if not P:
        return
    A = np.array(P, dtype=float)
    ax.scatter(A[:, 0], A[:, 1], A[:, 2], s=s, c=col, marker=mark, label=lab,
               depthshade=False, alpha=alpha, edgecolors=ec if ec else "none", linewidths=0.6)


sc(dead, DEAD, 16, "o", "on a wall both readings see", 0.55)
sc(surv, SURV, 30, "o", r"$\tau_t\neq0$", 1.0)
sc(disc, DISC, 74, "D", r"dies on $2a_i\equiv0$; the dual reads it as regular", 1.0, INK)

# los tres planos a_i = t/2, DESPUES de los puntos
g = np.linspace(0, t - 1, 2)
G1, G2 = np.meshgrid(g, g)
half = t / 2.0
for k in range(3):
    if k == 0:
        ax.plot_surface(np.full_like(G1, half), G1, G2, color=DISC, alpha=0.07, shade=False)
    elif k == 1:
        ax.plot_surface(G1, np.full_like(G1, half), G2, color=DISC, alpha=0.07, shade=False)
    else:
        ax.plot_surface(G1, G2, np.full_like(G1, half), color=DISC, alpha=0.07, shade=False)

ax.set_xlabel(r"$a_1 \ \mathrm{mod}\ t$", labelpad=1)
ax.set_ylabel(r"$a_2 \ \mathrm{mod}\ t$", labelpad=1)
ax.set_zlabel(r"$a_3 \ \mathrm{mod}\ t$", labelpad=1)
for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    axis.set_ticks(range(0, t, 2))
    axis.label.set_color(SECOND if (SECOND := "#52514e") else INK)
ax.tick_params(colors=MUTED, labelsize=7.5)
ax.set_xlim(-0.4, t - 0.6); ax.set_ylim(-0.4, t - 0.6); ax.set_zlim(-0.4, t - 0.6)
ax.view_init(elev=20, azim=-58)
ax.set_box_aspect((1, 1, 1))
for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
    pane.set_facecolor(SURFACE); pane.set_edgecolor(MUTED); pane.set_alpha(0.25)
ax.grid(True, color=MUTED, alpha=0.18, linewidth=0.5)

leg = ax.legend(loc="upper left", bbox_to_anchor=(-0.06, 0.99), frameon=False,
                fontsize=8.2, handletextpad=0.5, borderaxespad=0.0, labelspacing=0.55)
for txt in leg.get_texts():
    txt.set_color(INK)

ax.set_title(r"$C_3$ at $t=8$: the wall $2a_i\equiv0$ is invisible in the dual",
             color=INK, fontsize=10.5, pad=2)
fig.tight_layout()
fig.savefig("fig_walls3d.pdf", bbox_inches="tight")
print("escrito fig_walls3d.pdf")
