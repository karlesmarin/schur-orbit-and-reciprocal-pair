"""
Figure: the filter's locus is a union of Galois orbits, and the two types differ by one of them.

Authors: Carles Marin, Claude (AI assistant).

Three things that are one thing, drawn so the eye can check them.

LEFT.  The regular locus in exponent coordinates, for t = 7 at rank two, coloured by the orbit of
(Z/t)^x acting by x -> kx. It is a union of orbits, never a partial one. That is the visible form of
a fact about values: the filter takes values in {0,+-1} c Q, so it is fixed by
Gal(Q(zeta_t)/Q) = (Z/t)^x, and its vanishing locus has to be stable. The practical consequence is
the label: the filter only has to be evaluated ONCE PER COLOUR, a saving of phi(t) as soon as the
action is free.

MIDDLE.  The same picture for the two rules that the paper actually uses, at t = 5 and rank two.
Filled squares are the weights the type-B rule keeps, open circles the ones the type-C rule keeps.
They are different sets -- that is the disagreement measured in Figure 2. The arrows are the
translation by (t-1)/2 in every coordinate, and every arrow lands on a circle: the two sets are one
set, moved. The translation is multiplication by 2^{-1} in exponent coordinates, so it exists exactly
when 2 is a unit, that is exactly when t is odd.

RIGHT.  Why that is a statement about parity and not about t = 5. The sizes of the two loci at rank
two, for 3 <= t <= 14. They agree for every odd t and differ for every even t (t = 4 excepted, where
both are empty and nothing is being compared). No translation can fix a difference of cardinality,
so the even case fails for a reason no cleverness removes.

Everything is computed from the definitions below; nothing is placed by hand. Two controls print on
every run: that no non-unit preserves the locus, and that no translation other than (t-1)/2 works.
"""
import itertools
from math import gcd

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
PLUS = "#2a78d6"      # blue
MINUS = "#eb6834"     # orange
WALL = "#b9b7ae"
GHOST = "#cfccc2"
# hues for the Galois orbits: the two validated categorical ones first, then neutrals
ORBIT = [PLUS, MINUS, "#4a9a6a", "#8a6bb0", "#b8952f", "#5f5b54"]

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})


def regular(v, t):
    """no coordinate zero, and the 2n classes +-v_i pairwise distinct."""
    cl = []
    for x in v:
        cl.append(x % t)
        cl.append((-x) % t)
    return (0 not in cl) and len(set(cl)) == len(cl)


def locus_exponentes(t, n):
    """el lugar regular en coordenadas de EXPONENTE: ahi es donde actuan las unidades."""
    return set(v for v in itertools.product(range(t), repeat=n) if regular(v, t))


# Los dos lugares siguientes estan en coordenadas de PESO, que es donde vive la traslacion y donde
# se dibuja la Figura 2.  Mezclarlas con las de exponente fue el primer error de este guion, y el
# control lo cazo: la traslacion buena parecia fallar en los cuatro t impares a la vez.
def locus_C(t, n):
    """type C, weight coordinates: a_i = eta_i + rho_{C,i}  con rho_{C,i} = n - i + 1."""
    return set(e for e in itertools.product(range(t), repeat=n)
               if regular([e[i] + (n - i) for i in range(n)], t))


def locus_B(t, n):
    """type B, weight coordinates: A_i = 2 eta_i + 2(n - i) + 1."""
    return set(e for e in itertools.product(range(t), repeat=n)
               if regular([2 * e[i] + 2 * (n - i - 1) + 1 for i in range(n)], t))


def orbitas(S, t):
    """the orbits of (Z/t)^x acting coordinatewise by multiplication."""
    U = [k for k in range(1, t) if gcd(k, t) == 1]
    vistos, out = set(), []
    for v in sorted(S):
        if v in vistos:
            continue
        orb = set(tuple((k * x) % t for x in v) for k in U)
        out.append(sorted(orb))
        vistos |= orb
    return out


# ----------------------------------------------------------------- controles
def controles():
    fallos = []
    # C1: ninguna NO-unidad preserva el lugar
    malos = 0
    for t in range(3, 15):
        for n in (1, 2):
            S = locus_exponentes(t, n)
            if not S:
                continue
            for k in range(1, t):
                if gcd(k, t) == 1:
                    continue
                if set(tuple((k * x) % t for x in v) for v in S) == S:
                    malos += 1
    if malos:
        fallos.append("una NO-unidad preserva el lugar (%d casos)" % malos)
    # C2: para t impar, SOLO la traslacion (t-1)/2 lleva el lugar B al lugar C
    otras = 0
    for t in (5, 7, 9, 11):
        for n in (2,):
            B, C = locus_B(t, n), locus_C(t, n)
            if not B:
                continue
            for k in range(t):
                ok = set(tuple((x + k) % t for x in u) for u in B) == C
                if ok and k != (t - 1) // 2:
                    otras += 1
                if not ok and k == (t - 1) // 2:
                    fallos.append("la traslacion buena FALLA en t=%d n=%d" % (t, n))
    if otras:
        fallos.append("otra traslacion tambien funciona (%d casos)" % otras)
    print("   control  ninguna no-unidad preserva el lugar         : %s" % ("ok" if not malos else "FALLA"))
    print("   control  ninguna otra traslacion funciona            : %s" % ("ok" if not otras else "FALLA"))
    if fallos:
        raise SystemExit("CONTROL FALLIDO: " + "; ".join(fallos))


controles()

fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.75))

# ============================================================ LEFT: Galois orbits
tL, nL = 7, 2
ax = axes[0]
S = locus_exponentes(tL, nL)
orbs = orbitas(S, tL)
for i, orb in enumerate(orbs):
    xs = [p[0] for p in orb]
    ys = [p[1] for p in orb]
    ax.scatter(xs, ys, s=62, color=ORBIT[i % len(ORBIT)], zorder=3,
               edgecolors=SURFACE, linewidths=0.8)
    # el representante de la orbita, marcado: es el unico punto que hay que calcular
    ax.scatter([orb[0][0]], [orb[0][1]], s=170, facecolors="none", zorder=4,
               edgecolors=ORBIT[i % len(ORBIT)], linewidths=1.4)
todos = set(itertools.product(range(tL), repeat=2))
fuera = todos - S
ax.scatter([p[0] for p in fuera], [p[1] for p in fuera], s=16, color=GHOST, zorder=2)
ax.set_xticks(range(tL)); ax.set_yticks(range(tL))
ax.set_xlim(-0.7, tL - 0.3); ax.set_ylim(-0.7, tL - 0.3)
ax.set_xlabel(r"$a_1$ mod $t$"); ax.set_ylabel(r"$a_2$ mod $t$")
ax.set_title(r"$t=7$: the locus is $%d$ Galois orbits" % len(orbs), fontsize=9.5, color=INK)
ax.text(0.5, -0.30, r"$%d$ points, $%d$ orbits: one evaluation per colour, a saving of $\varphi(t)=%d$"
        % (len(S), len(orbs), len([k for k in range(1, tL) if gcd(k, tL) == 1])),
        transform=ax.transAxes, ha="center", fontsize=8, color=SECOND)
ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

# ============================================================ MIDDLE: the translation
tM, nM = 5, 2
ax = axes[1]
B, C = locus_B(tM, nM), locus_C(tM, nM)
k = (tM - 1) // 2
ax.scatter([p[0] for p in C], [p[1] for p in C], s=210, facecolors="none",
           edgecolors=MINUS, linewidths=1.6, zorder=3)
ax.scatter([p[0] for p in B], [p[1] for p in B], s=58, color=PLUS, marker="s", zorder=4)
for (a, b) in B:
    a2, b2 = (a + k) % tM, (b + k) % tM
    # se dibuja la flecha solo si no cruza el borde, para no ensuciar con saltos del toro
    if a + k < tM and b + k < tM:
        ax.annotate("", xy=(a2, b2), xytext=(a, b), zorder=5,
                    arrowprops=dict(arrowstyle="->", color=SECOND, linewidth=0.9,
                                    shrinkA=5, shrinkB=9))
    else:
        ax.scatter([a2], [b2], s=26, color=SECOND, marker="x", zorder=5, linewidths=1.1)
todos = set(itertools.product(range(tM), repeat=2))
fuera = todos - B - C
ax.scatter([p[0] for p in fuera], [p[1] for p in fuera], s=16, color=GHOST, zorder=2)
ax.set_xticks(range(tM)); ax.set_yticks(range(tM))
ax.set_xlim(-0.7, tM - 0.3); ax.set_ylim(-0.7, tM - 0.3)
ax.set_xlabel(r"$\eta_1$ mod $t$"); ax.set_ylabel(r"$\eta_2$ mod $t$")
ax.set_title(r"$t=5$: type $B$ moved by $\frac{t-1}{2}$ is type $C$", fontsize=9.5, color=INK)
ax.text(0.5, -0.30, r"every arrow lands on a circle; $\times$ marks a jump across the torus",
        transform=ax.transAxes, ha="center", fontsize=8, color=SECOND)
ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_axisbelow(True)
ax.legend(handles=[Line2D([], [], marker="s", color=PLUS, linestyle="", markersize=6,
                          label=r"type $B$ keeps"),
                   Line2D([], [], marker="o", markerfacecolor="none", color=MINUS,
                          linestyle="", markersize=9, label=r"type $C$ keeps")],
          loc="upper left", frameon=False, fontsize=7.5, handletextpad=0.4,
          bbox_to_anchor=(-0.02, 1.02))

# ============================================================ RIGHT: the parity in the counts
ax = axes[2]
TS = list(range(3, 15))
nb = [len(locus_B(t, 2)) for t in TS]
nc = [len(locus_C(t, 2)) for t in TS]
ax.plot(TS, nb, marker="s", color=PLUS, linewidth=1.3, markersize=5,
        label=r"$\#\{\eta: A_B(\eta)\ \mathrm{regular}\}$")
ax.plot(TS, nc, marker="o", color=MINUS, linewidth=1.3, markersize=5,
        markerfacecolor="none", label=r"$\#\{\eta: a_C(\eta)\ \mathrm{regular}\}$")
for t, a, b in zip(TS, nb, nc):
    if t % 2 == 0 and (a or b):
        ax.annotate("", xy=(t, a), xytext=(t, b), zorder=2,
                    arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=3.0, alpha=0.35))
ax.set_xticks(TS)
ax.set_xlabel(r"$t$"); ax.set_ylabel("size of the locus, rank $2$")
ax.set_title("equal for odd $t$; unequal in all tested even cases", fontsize=9.5, color=INK)
ax.text(0.5, -0.30,
        r"weight coordinates: $t=4$ is the one even case with nothing to compare, both sets empty",
        transform=ax.transAxes, ha="center", fontsize=8, color=SECOND)
ax.legend(loc="upper left", frameon=False, fontsize=7.5)
ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_axisbelow(True)
for a in axes:
    for s in a.spines.values():
        s.set_color(MUTED)
        s.set_linewidth(0.6)

fig.tight_layout()
fig.subplots_adjust(bottom=0.22)
fig.savefig("fig_galois.pdf", bbox_inches="tight")
print("   fig_galois.pdf escrito")
print("   izquierda t=%d: %d puntos en %d orbitas" % (tL, len(S), len(orbs)))
print("   centro    t=%d: |B| = %d, |C| = %d, traslacion por %d" % (tM, len(B), len(C), k))
print("   derecha   pares con cardinales distintos: %s"
      % [t for t, a, b in zip(TS, nb, nc) if t % 2 == 0 and a != b])
