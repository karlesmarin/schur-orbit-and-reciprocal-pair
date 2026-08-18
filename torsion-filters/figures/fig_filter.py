"""
Figure: the torsion filter drawn as what it is -- an arrangement of walls on the weight lattice.

Authors: Carles Marin, Claude (AI assistant).

The filter says a weight dies when its shifted vector meets a wall, and otherwise survives with a
sign. In rank two the walls are LINES in the eta-plane, and once they are drawn the statement needs
no prose: the survivors are the lattice points that miss every line.

The two panels are NOT the same rule, and an earlier version of this figure drew them as if they
were. For even t the shift is a = eta + rho_C, integral, and the walls are

  a_1 = 0, a_2 = 0, a_1 = t/2, a_2 = t/2, a_1 = +-a_2   (mod t)

For odd t the group is B, rho is half-integral, and the right variable is A = 2(eta + rho_B),
that is A_1 = 2 eta_1 + 3 and A_2 = 2 eta_2 + 1 at rank two; the walls are A_j = 0 and A_1 = +-A_2,
and the two middle families do not exist at all.  Drawing the odd panel with the even rule gives, by
a coincidence that hides the error, the SAME number of survivors -- 49 of 136 -- while differing on
66 of the 136 points, 33 in each direction.

The odd sign is not drawn from a formula either: the natural closed form fails (it is right at
m'=1 and m'=3 and wrong on 13 of 36 weights at m'=2), so the odd panel evaluates the type-B
bialternant exactly in Z[x]/(Phi_t(x)).  Two controls print on every run: the exact value never
leaves {0,+-1}, and its support agrees with the regularity criterion computed independently.

Right panel: the survival rate against the rank, for both parities, in a fixed box.

Everything is computed from the rule; nothing is placed by hand.

Palette: the two validated categorical hues (blue = the fold gives +1, orange = -1), neutrals for
the arrangement. Marker shape doubles the colour so the figure survives greyscale printing, as in the
companion paper.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"; GRID = "#e1e0d9"
PLUS = "#2a78d6"      # blue   -- the fold gives +1
MINUS = "#eb6834"     # orange -- the fold gives -1
WALL = "#b9b7ae"      # the arrangement
GHOST = "#cfccc2"     # the families that do not exist for odd t

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})


def filtro_C(eta, t, m):
    """el filtro PAR, tipo C: la regla probada del Lema 3.1, con su signo cerrado."""
    a = [eta[j] + (m - (j + 1) + 1) for j in range(m)]
    for v in a:
        c = v % t
        if c == 0 or 2 * c == t:
            return 0
    cl, sg = [], 1
    for v in a:
        c = v % t
        if c <= m:
            cl.append(c)
        else:
            cl.append(t - c); sg *= -1
    if len(set(cl)) != m:
        return 0
    perm = [m - cl[j] for j in range(m)]
    inv = sum(1 for i in range(m) for j in range(i + 1, m) if perm[i] > perm[j])
    return sg * (-1) ** inv


# ---------------------------------------------------------------------------------------------
# EL FILTRO IMPAR SE CALCULA, NO SE ADIVINA.  Una version anterior de esta figura dibujaba el panel
# t=5 con la regla SIMPLECTICA de arriba; da los mismos 49 supervivientes por casualidad y difiere
# en 66 de los 136 puntos.  Y la formula cerrada que probamos para el signo impar acierta en
# m'=1 y m'=3 pero solo 23 de 36 en m'=2, asi que tampoco se puede usar.  Se evalua el bialternante
# de tipo B exactamente, en Z[x]/(Phi_t(x)), con  A_j = 2 eta_j + 2(m'-j) + 1.
from sympy import Poly, ZZ as _ZZ, cyclotomic_poly, symbols as _sym

_x = _sym('x')
_anillo = {}


def _ring(t):
    if t not in _anillo:
        phi = Poly(cyclotomic_poly(t, _x), _x, domain=_ZZ)
        _anillo[t] = (phi, {})
    return _anillo[t]


def _xp(t, k):
    phi, cache = _ring(t)
    k %= t
    if k not in cache:
        cache[k] = Poly(_x ** k, _x, domain=_ZZ).rem(phi)
    return cache[k]


def _det2(M, phi):
    return (M[0][0] * M[1][1] - M[0][1] * M[1][0]).rem(phi)


def filtro_B(eta, t, m):
    """el filtro IMPAR, tipo B: cociente de bialternantes, exacto.  Solo rango 2, que es el que
       la figura dibuja."""
    assert m == 2, "el panel impar de esta figura es de rango 2"
    phi, _ = _ring(t)
    A = [2 * eta[j] + 2 * (m - j - 1) + 1 for j in range(m)]
    R = [2 * 0 + 2 * (m - j - 1) + 1 for j in range(m)]        # eta = 0

    def num(b):
        M = [[(_xp(t, (i + 1) * b[j]) - _xp(t, -(i + 1) * b[j])).rem(phi) for j in range(m)]
             for i in range(m)]
        return _det2(M, phi)

    den = num(R)
    if den.is_zero:
        return 0
    q = num(A)
    if q.is_zero:
        return 0
    # el cociente es +-1; se decide comparando q con +-den
    if (q - den).rem(phi).is_zero:
        return 1
    if (q + den).rem(phi).is_zero:
        return -1
    return 2      # no deberia ocurrir; se marca para que el control lo cace


def sobrevive(eta, t, m, es_par):
    """la CONDICION DE REGULARIDAD del Lema 3.4, valida en las dos paridades y en todo rango.
       No da el signo -- y por eso se usa solo donde hace falta la tasa de supervivencia."""
    if es_par:
        b = [eta[j] + (m - j) for j in range(m)]
        largo = [(2 * v) % t != 0 for v in b]
    else:
        b = [2 * eta[j] + 2 * (m - j - 1) + 1 for j in range(m)]
        largo = [v % t != 0 for v in b]
    if not all(largo):
        return False
    if es_par and any(v % t == 0 for v in b):
        return False
    for i in range(m):
        for j in range(i + 1, m):
            if (b[i] - b[j]) % t == 0 or (b[i] + b[j]) % t == 0:
                return False
    return True


def paredes_B(ax, t, TOPE):
    """el arreglo IMPAR, en las coordenadas A_j = 2 eta_j + 2(m'-j) + 1 con m' = 2:
         A_1 = 2 eta_1 + 3,  A_2 = 2 eta_2 + 1.
       A_1 = 0 (mod t)  ->  eta_1 = (t-3)/2 * inverso de 2 ... se resuelve barriendo, que es exacto
       y no se equivoca de inverso."""
    lo, hi = -1.0, TOPE + 1.0
    col, lw, ls, z = WALL, 0.9, "-", 1
    for e1 in range(-1, TOPE + 2):
        if (2 * e1 + 3) % t == 0:
            ax.axvline(e1, color=col, linewidth=lw, linestyle=ls, zorder=z)
    for e2 in range(-1, TOPE + 2):
        if (2 * e2 + 1) % t == 0:
            ax.axhline(e2, color=col, linewidth=lw, linestyle=ls, zorder=z)
    for e1 in range(-1, TOPE + 2):
        for e2 in range(-1, TOPE + 2):
            pass
    # A_1 = +- A_2  ->  2(eta_1 - eta_2) + 2 = 0  y  2(eta_1 + eta_2) + 4 = 0  (mod t)
    for d in range(-TOPE - t, TOPE + t + 1):
        if (2 * d + 2) % t == 0:
            ax.plot([lo, hi], [lo - d, hi - d], color=col, linewidth=lw, linestyle=ls, zorder=z)
    for sm in range(-1, 2 * TOPE + t + 2):
        if (2 * sm + 4) % t == 0:
            ax.plot([lo, hi], [sm - lo, sm - hi], color=col, linewidth=lw, linestyle=ls, zorder=z)


def paredes(ax, t, TOPE, ghost=False):
    """dibuja el arreglo PAR.  ghost=True: las dos familias de t/2, que en t impar NO existen."""
    lo, hi = -1.0, TOPE + 1.0
    col = GHOST if ghost else WALL
    lw = 0.85 if ghost else 0.9
    ls = (0, (3, 3)) if ghost else "-"
    z = 1
    if not ghost:
        for c in range(-2, TOPE + t, t):            # a_1 = 0
            if lo <= c <= hi:
                ax.axvline(c, color=col, linewidth=lw, linestyle=ls, zorder=z)
        for c in range(-1, TOPE + t, t):            # a_2 = 0
            if lo <= c <= hi:
                ax.axhline(c, color=col, linewidth=lw, linestyle=ls, zorder=z)
        for d in range(-1 - 3 * t, TOPE + 3 * t, t):     # a_1 = +a_2
            ax.plot([lo, hi], [lo - d, hi - d], color=col, linewidth=lw, linestyle=ls, zorder=z)
        for sm in range(-3, TOPE * 2 + 3 * t, t):        # a_1 = -a_2
            ax.plot([lo, hi], [sm - lo, sm - hi], color=col, linewidth=lw, linestyle=ls, zorder=z)
    else:
        h = t // 2
        for c in range(h - 2, TOPE + t, t):
            if lo <= c <= hi:
                ax.axvline(c, color=col, linewidth=lw, linestyle=ls, zorder=z)
        for c in range(h - 1, TOPE + t, t):
            if lo <= c <= hi:
                ax.axhline(c, color=col, linewidth=lw, linestyle=ls, zorder=z)


fig = plt.figure(figsize=(9.6, 3.45))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.02], wspace=0.28)

TOPE = 15
for k, (t, titulo, es_par) in enumerate([(6, r"even, $t=6$   (three wall families)", True),
                                         (5, r"odd, $t=5$   (two: the middle pair is absent)", False)]):
    m = (t - 2) // 2 if es_par else (t - 1) // 2
    ax = fig.add_subplot(gs[0, k])
    if not es_par:
        paredes(ax, t, TOPE, ghost=True)      # primero el fantasma, debajo
        paredes_B(ax, t, TOPE)
    else:
        paredes(ax, t, TOPE, ghost=False)
    P, M, D = [], [], []
    for e1 in range(TOPE + 1):
        for e2 in range(e1 + 1):
            v = filtro_C((e1, e2), t, m) if es_par else filtro_B((e1, e2), t, m)
            (P if v == 1 else M if v == -1 else D).append((e1, e2))
    if D:
        X, Y = zip(*D)
        ax.scatter(X, Y, s=3.5, c=MUTED, marker=".", linewidths=0, zorder=2, alpha=0.55)
    for L, col, mk, sz in [(P, PLUS, "o", 26), (M, MINUS, "s", 23)]:
        if L:
            X, Y = zip(*L)
            ax.scatter(X, Y, s=sz, c=col, marker=mk, linewidths=0, zorder=4)
    ax.fill([-1, TOPE + 1, -1], [-1, TOPE + 1, TOPE + 1], color=SURFACE, zorder=5, alpha=0.93)
    ax.plot([-1, TOPE + 1], [-1, TOPE + 1], color=MUTED, linewidth=0.6, zorder=6)
    ax.text(0.245, 0.775, r"not dominant", transform=ax.transAxes, fontsize=7.5,
            color=MUTED, zorder=7, rotation=45, ha="center")
    ax.set_xlabel(r"$\eta_1$"); ax.set_ylabel(r"$\eta_2$")
    ax.set_title(titulo, fontsize=8.5, color=INK)
    ax.set_xlim(-1, TOPE + 1); ax.set_ylim(-1, TOPE + 1)
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_color(MUTED); s.set_linewidth(0.6)
    n = len(P) + len(M) + len(D)
    ax.text(0.965, 0.055, r"%d of %d survive" % (len(P) + len(M), n), transform=ax.transAxes,
            fontsize=8, color=SECOND, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.22", fc=SURFACE, ec="none", alpha=0.92))

# ---- panel derecho ----------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 2])
COTA = 9
for parid, col, mk, lab in [("par", PLUS, "o", r"even $t=2m+2$"), ("impar", MINUS, "s", r"odd $t=2m+1$")]:
    xs, ys = [], []
    for m in range(1, 6):
        t = 2 * m + 2 if parid == "par" else 2 * m + 1
        tot = [0, 0]
        def rec(pref, k, cap):
            if k == m:
                tot[0] += 1
                if sobrevive(tuple(pref), t, m, parid == "par"):
                    tot[1] += 1
                return
            for v in range(cap, -1, -1):
                rec(pref + [v], k + 1, v)
        rec([], 0, COTA)
        xs.append(m); ys.append(100.0 * tot[1] / tot[0])
    ax.plot(xs, ys, marker=mk, color=col, linewidth=1.4, markersize=5, label=lab, zorder=3)
    # Las dos series con la etiqueta ARRIBA se pisan: la curva par va siempre por debajo de la
    # impar, y con el mismo desplazamiento su numero cae justo sobre el marcador de la otra --- el
    # 8 de m=4 quedaba con 0 % de tinta visible, medido por check_text_over_art.py.  Cada serie se
    # rotula hacia su lado libre, y el texto va por encima de lo dibujado.
    # Y el ultimo punto de la par cae tan cerca de cero que el rotulo de abajo se sentaba SOBRE el
    # eje y pegado al "5" de la escala: ahi la unica direccion libre es arriba a la izquierda,
    # porque la curva llega desde arriba y la impar va por encima.
    tope = max(ys)
    for x, y in zip(xs, ys):
        if parid == "par":
            xy, ha = ((-12, 2), "right") if y < 0.10 * tope else ((0, -11), "center")
        else:
            xy, ha = (0, 7), "center"
        ax.annotate("%.0f" % y, (x, y), textcoords="offset points", xytext=xy,
                    ha=ha, fontsize=7, color=col, zorder=6)
ax.set_xlabel(r"rank $m$ of the frozen factor")
# El rotulo NO pasa por LaTeX --- la figura usa mathtext, no usetex ---, asi que un `\%` se
# imprime con la barra.  Estuvo asi en el PDF hasta el 18 de agosto de 2026.
ax.set_ylabel(r"surviving weights (%)")
ax.set_title(r"the collapse with the rank", fontsize=8.5, color=INK)
ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_xticks([1, 2, 3, 4, 5]); ax.set_ylim(bottom=0)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.6)
ax.legend(frameon=False, fontsize=8, loc="upper right")
ax.text(0.035, 0.06, r"box $\eta_1\leq %d$" % COTA, transform=ax.transAxes, fontsize=8, color=MUTED)

leyenda = [Line2D([], [], color=PLUS, marker="o", linestyle="none", markersize=5,
                  label=r"survives, $\tau_t=+1$"),
           Line2D([], [], color=MINUS, marker="s", linestyle="none", markersize=5,
                  label=r"survives, $\tau_t=-1$"),
           Line2D([], [], color=MUTED, marker=".", linestyle="none", markersize=6,
                  label=r"on a wall, $\tau_t=0$"),
           Line2D([], [], color=WALL, linewidth=1.0, label=r"wall"),
           Line2D([], [], color=GHOST, linewidth=1.0, linestyle=(0, (3, 3)),
                  label=r"wall absent for odd $t$")]
fig.legend(handles=leyenda, frameon=False, fontsize=8, ncol=5,
           loc="lower center", bbox_to_anchor=(0.5, -0.155))

fig.savefig("fig_filter.pdf", bbox_inches="tight", pad_inches=0.02)
print("fig_filter.pdf escrito")
# CONTROLES DE LA FIGURA, impresos en cada corrida:
#  - el filtro impar exacto no puede devolver nunca el centinela 2 (|tau| != 1 seria un fallo);
#  - y su soporte tiene que coincidir con la condicion de regularidad, que es otra ruta.
for t in (5, 6):
    es_par = (t % 2 == 0)
    m = (t - 2) // 2 if es_par else (t - 1) // 2
    n = v = mal = disc = 0
    for e1 in range(TOPE + 1):
        for e2 in range(e1 + 1):
            n += 1
            f = filtro_C((e1, e2), t, m) if es_par else filtro_B((e1, e2), t, m)
            if abs(f) > 1:
                mal += 1
            if f:
                v += 1
            if bool(f) != sobrevive((e1, e2), t, m, es_par):
                disc += 1
    print("  t=%d (%s): %d de %d sobreviven (%.1f %%) | |tau|>1: %d | discrepa con regularidad: %d"
          % (t, "par" if es_par else "impar", v, n, 100.0 * v / n, mal, disc))
    assert mal == 0 and disc == 0, "control de la figura roto en t=%d" % t
