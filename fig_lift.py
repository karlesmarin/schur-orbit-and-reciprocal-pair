"""
Figure: what the numerator lift does to the count of fixed points.

Authors: Carles Marin, Claude (AI assistant).

Section 6.1 turns on one contrast, and until now it lived only in prose. By Lemma 6.6, any
weight-preserving sign-reversing involution on J_lambda leaves at least ||Phi_t(lambda;z)||_1 fixed
points, and on the lift J^_lambda at least ||D_t Phi_t(lambda;z)||_1. The first is unbounded --- it
reaches 30 at t = 2 --- while Laplace offers at most eight monomials; the second is exactly
2|U(lambda)|, so it takes only the values 0, 6 and 8. The lift does not move the requirement: it
collapses it onto the Laplace count.

That is what this figure draws, as two distributions over the same population of shapes. What to
look at is not the shape of either histogram but the horizontal extent of the two: the blue one runs
off to the right of the line at 8, the orange one cannot.

Everything is recomputed here from the definitions, and by a route that shares no code with the
gates: J_lambda is built as in equation (20) --- mu with empty t-core and every t-quotient component
of at most one row, then nu with nu/mu and lambda/nu horizontal strips --- and Phi is the signed sum
over it, so no root-of-unity arithmetic enters at all. The script prints its own counts, which are
the ones the verification table quotes, so that the figure can be checked against
`pI_numerator_lift.py` rather than trusted.

Palette: the two validated categorical hues, blue for the unlifted set and orange for the lift,
with hatching on the blue so the figure survives greyscale.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
UNLIFT = "#2a78d6"     # blue   -- the requirement on J_lambda
LIFT = "#eb6834"       # orange -- the requirement on the lift

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

LMAX = 7


def particiones(maxpart, maxlen):
    out = [()]

    def rec(pref, resto, tope):
        for p in range(min(tope, maxpart), 0, -1):
            nu = pref + [p]
            out.append(tuple(nu))
            if resto > 1:
                rec(nu, resto - 1, p)

    rec([], maxlen, maxpart)
    return sorted(set(out))


def beta_set(lam, n):
    l = list(lam) + [0] * (n - len(lam))
    return [l[i] + n - 1 - i for i in range(n)]


def part_from_beta(B):
    S = sorted(B, reverse=True)
    n = len(S)
    return tuple(x for x in [S[i] - (n - 1 - i) for i in range(n)] if x > 0)


def core_quot_sgn(lam, n, t):
    """t-core, t-quotient and the sign of the sorting permutation, by removing t-hooks."""
    B = set(beta_set(lam, n))
    sg = 1
    seguir = True
    while seguir:
        seguir = False
        for x in sorted(B, reverse=True):
            if x - t >= 0 and (x - t) not in B:
                sg *= (-1) ** len([y for y in B if x - t < y < x])
                B.discard(x)
                B.add(x - t)
                seguir = True
                break
    core = part_from_beta(B)
    orig = sorted(beta_set(lam, n), reverse=True)
    quot = []
    for r in range(t):
        xs = sorted([x for x in orig if x % t == r], reverse=True)
        m = len(xs)
        quot.append(tuple(x for x in [(xs[k] - r) // t - (m - 1 - k) for k in range(m)] if x > 0))
    return core, quot, sg


def sub(mu, lam):
    m = list(mu) + [0] * 10
    l = list(lam) + [0] * 10
    return all(m[i] <= l[i] for i in range(10))


def hstrip(sup, inf):
    a = list(sup) + [0] * 10
    b = list(inf) + [0] * 10
    if any(b[i] > a[i] for i in range(10)):
        return False
    return all(a[i + 1] <= b[i] for i in range(9))


def phi_de(lam, t, n, todas):
    """Phi_t(lambda; z) as an integer Laurent polynomial, by the sum over J_lambda of (20)."""
    acc = {}
    for mu in todas:
        if not sub(mu, lam):
            continue
        core, quot, sg = core_quot_sgn(mu, n, t)
        if len(core) != 0 or any(len(q) > 1 for q in quot):
            continue
        for nu in todas:
            if not (sub(mu, nu) and sub(nu, lam)):
                continue
            if not (hstrip(nu, mu) and hstrip(lam, nu)):
                continue
            e = 2 * sum(nu) - sum(mu) - sum(lam)
            acc[e] = acc.get(e, 0) + sg
    return {k: v for k, v in acc.items() if v}


def por_Dt(p, t):
    """multiply by D_t = (z^t - 1)(z^-t - 1)(z - 1/z)."""
    for factor in ({t: 1, 0: -1}, {-t: 1, 0: -1}, {1: 1, -1: -1}):
        o = {}
        for a in p:
            for b in factor:
                o[a + b] = o.get(a + b, 0) + p[a] * factor[b]
        p = {k: v for k, v in o.items() if v}
    return p


def norma1(p):
    return sum(abs(v) for v in p.values())


# ---------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.15))
resumen = []
for ax, t in zip(axes, (2, 3)):
    n = t + 2
    todas = particiones(LMAX, n)
    # SIN la particion vacia, para que la poblacion sea la misma que la de la corrida archivada
    # de `pI_numerator_lift.py`.  Incluyendola serian 330 y 792, y Phi(vacia) = 1 caeria en el
    # cubo 1 del histograma azul y en el 8 del naranja.  Que la figura y la tabla cuenten
    # poblaciones distintas es la clase de desfase que nadie ve hasta que alguien suma.
    lams = [l for l in todas if l]
    FJ, FK = Counter(), Counter()
    for lam in lams:
        phi = phi_de(lam, t, n, todas)
        FJ[norma1(phi)] += 1
        FK[norma1(por_Dt(phi, t))] += 1
    resumen.append((t, len(lams), dict(sorted(FJ.items())), dict(sorted(FK.items()))))

    xs = list(range(0, max(FJ) + 1))
    ax.bar(xs, [FJ.get(x, 0) for x in xs], width=0.78, color="none", edgecolor=UNLIFT,
           hatch="///", linewidth=0.8, zorder=5)
    ax.bar(xs, [FK.get(x, 0) for x in xs], width=0.78, color=LIFT, zorder=4)
    ax.axvline(8.5, color=MUTED, lw=0.9, ls=(0, (4, 3)), zorder=2)
    # SIN anotacion junto a la linea: arriba chocaba con los dos recuentos de la esquina y abajo se
    # montaba sobre las barras del panel de t=3.  Lo que la linea significa lo dice el pie, que es
    # donde las demas figuras de este articulo ponen su explicacion.  Un rotulo sobre tinta es
    # exactamente lo que `check_text_over_art` esta para cazar.
    fuera = sum(v for k, v in FJ.items() if k > 8)
    ax.set_title(r"$t=%d$,  %d shapes with $\lambda_1\leq %d$" % (t, len(lams), LMAX),
                 fontsize=9, color=INK, pad=5)
    ax.set_xlabel("minimum number of fixed points", color=INK)
    if t == 2:
        ax.set_ylabel("number of shapes", color=INK)
    ax.set_yscale("symlog", linthresh=1)
    ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_color(MUTED); sp.set_linewidth(0.7)
    ax.tick_params(colors=SECOND, labelsize=7.6)
    ax.set_xticks(list(range(0, max(FJ) + 1, 2 if t == 2 else 2)))
    ax.text(0.985, 0.955, r"$%d$ of $%d$ lie beyond $8$" % (fuera, len(lams)),
            transform=ax.transAxes, ha="right", va="top", fontsize=7.4, color=UNLIFT)
    ax.text(0.985, 0.865, r"the lift takes only $%s$"
            % ",\\,".join(str(k) for k in sorted(FK)), transform=ax.transAxes,
            ha="right", va="top", fontsize=7.4, color=LIFT)

h = [plt.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor=UNLIFT, hatch="///", linewidth=0.8,
                   label=r"on $J_\lambda$: $\|\Phi_t(\lambda;z)\|_1$"),
     plt.Rectangle((0, 0), 1, 1, facecolor=LIFT, edgecolor="none",
     # SIN el «= 2|U|» en la leyenda.  Lo tenia, y era FALSO en las veinte concentricas: alli
     # |U| = 4, luego 2|U| = 8, mientras la norma es 0.  La Prop. 6.10 pide Phi != 0 y la leyenda
     # no lo decia.  El pie lo dice ahora, y la leyenda se queda con lo que se dibuja.
                   label=r"on the lift $\widehat J_\lambda$: $\|D_t\Phi_t(\lambda;z)\|_1$")]
fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, fontsize=7.8,
           bbox_to_anchor=(0.5, -0.03))
fig.tight_layout(rect=(0, 0.075, 1, 1))

out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fig_lift")
fig.savefig(out + ".pdf")
fig.savefig(out + ".png", dpi=170)
print("wrote %s.pdf / .png" % out)
print()
print("los numeros que el articulo cita, recomputados por otro camino:")
for t, nl, fj, fk in resumen:
    print("  t=%d  %d formas" % (t, nl))
    print("     F_J = ||Phi||_1        : %s   maximo %d" % (fj, max(fj)))
    print("     F_K = ||D_t.Phi||_1    : %s   maximo %d" % (fk, max(fk)))
    print("     formas con F_J > 8     : %d" % sum(v for k, v in fj.items() if k > 8))
