# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Figuras de la nota conductor_rho_line.  Todas salen de datos calculados:
#   fig_historia  -- tres corrientes y dónde se cruzan (fechas de publicación verificadas)
#   fig_paisaje   -- 3D: log10 del índice [O+:A_m(q)] por primo, coloreado por la familia h, h+1, h+2
#   fig_perfil    -- 3D: perfil W_i/total en casillas con p | N (escalón de von Staudt en p-1)
#   fig_clase     -- déficit de W_1 frente a (q' primo, p), con h^-(Q(zeta_q')) anotado
#   fig_sierra    -- S(c) de los dos tipos: diente de sierra y diente de sierra desplazado medio periodo
# Paleta: referencia dataviz validada (azul/naranja/aqua, light).  Tinta y rejilla recesivas.
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import Patch
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "gates", "analisis"))

SUP = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASE = "#c3c2b7"
C_H, C_H1, C_H2 = "#2a78d6", "#eb6834", "#1baf7a"
FAM = {"h": C_H, "h+1": C_H1, "h+2": C_H2}
SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.edgecolor": BASE, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.facecolor": SUP, "figure.facecolor": SUP,
    "savefig.facecolor": SUP, "axes.grid": False, "legend.frameon": False,
})


IDIOMA = "en"  # "en" -> fig_x.pdf ;  "es" -> fig_x_es.pdf (edicion en castellano)


BANDA = "#eef2f6"


def T(en, es):
    return es if IDIOMA == "es" else en


def guardar(fig, nombre, preview_dir=None):
    if IDIOMA != "en":
        nombre = nombre + "_" + IDIOMA
    fig.savefig(os.path.join(AQUI, nombre + ".pdf"), bbox_inches="tight")
    if preview_dir:
        fig.savefig(os.path.join(preview_dir, nombre + ".png"), dpi=110, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- historia
def fig_historia(prev):
    carriles = [
        (T("Bernoulli numbers\nand class numbers", "números de Bernoulli\ny números de clase"), [
            (1840, "von Staudt–Clausen"), (1850, "Kummer"), (1890, "Stickelberger"),
            (1955, "Carlitz–Olson\n(Maillet)"), (1962, "Iwasawa"), (1978, "Sinnott")]),
        (T("conductors of\none-dimensional rings", "conductores de anillos\nunidimensionales"), [
            (1877, "Dedekind\n(Führer)"), (1963, "Bass"), (1971, "Herzog–Kunz"), (2006, "Huneke–Swanson")]),
        (T("character values\nat torsion elements", "valores de caracteres\nen elementos de torsión"), [
            (1940, "Littlewood"), (1976, "Kostant"), (2009, T("Amiot\n(musical scales)", "Amiot\n(escalas musicales)")),
            (2019, T("Gannon–Schopieray\n(fields)", "Gannon–Schopieray\n(cuerpos)")), (2020, T("Bächle–Sambale\n(orders, finite G)", "Bächle–Sambale\n(órdenes, G finito)")),
            (2026, T("integer points\non the ρ-line", "puntos enteros\nen la recta ρ"))]),
    ]
    # posicion explicita (x de la etiqueta, desplazamiento vertical) por evento y carril: ningun texto
    # pisa a otro al tamano de pagina.  Las etiquetas desplazadas en x llevan una guia fina.
    pos = {
        (0, 1840): (1840, 0.28), (0, 1850): (1856, -0.28), (0, 1890): (1890, 0.28),
        (0, 1955): (1950, -0.28), (0, 1962): (1962, 0.28), (0, 1978): (1983, -0.28),
        (1, 1877): (1877, 0.28), (1, 1963): (1960, 0.28), (1, 1971): (1974, -0.28), (1, 2006): (1995, 0.28),
        (2, 1940): (1940, 0.28), (2, 1976): (1976, 0.28), (2, 2009): (1990, -0.50),
        (2, 2019): (2000, -1.10), (2, 2020): (2010, 0.62), (2, 2026): (2021, -0.50),
    }
    if IDIOMA == "es":  # las etiquetas en castellano son mas largas
        pos[(2, 2009)] = (1983, -0.50)
        pos[(2, 2026)] = (2025, -0.50)
    colores = [C_H1, C_H2, C_H]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ys = [2.7, 1.35, 0]
    for k, ((nombre, ev), y, col) in enumerate(zip(carriles, ys, colores)):
        ax.plot([1830, 2032], [y, y], color=GRID, lw=1, zorder=0)
        ax.text(1826, y, nombre, ha="right", va="center", color=INK2, fontsize=8.5)
        for anio, etq in ev:
            ax.plot(anio, y, "o", ms=6, color=col, mec=SUP, mew=1.5, zorder=3)
            xl, dy = pos[(k, anio)]
            if abs(xl - anio) > 3:
                ax.plot([anio, xl], [y, y + dy * 0.75], color=BASE, lw=0.6, zorder=1)
            ax.text(xl, y + dy, etq, ha="center", va="bottom" if dy > 0 else "top", fontsize=7.5, color=INK)
    # la convergencia
    ys_star = ys[1]
    ax.plot(2031, ys_star, marker="*", ms=16, color=INK, zorder=4)
    ax.text(2031, ys_star + 0.30, T("this note", "esta nota"), ha="center", va="bottom", fontsize=8.5, color=INK, weight="bold")
    for y in ys:
        if y == ys_star:
            continue
        ax.annotate("", xy=(2030.5, ys_star + (0.12 if y > ys_star else -0.12)), xytext=(2024, y),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1))
    ax.set_xlim(1830, 2036)
    ax.set_ylim(-1.55, 3.35)
    ax.set_yticks([])
    ax.set_xticks([1850, 1900, 1950, 2000])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    guardar(fig, "fig_historia", prev)


# ---------------------------------------------------------------- paisaje 3D
def fig_paisaje(prev):
    datos = json.load(open(os.path.join(AQUI, "datos_indice.json")))
    fig = plt.figure(figsize=(7.2, 4.8))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor(SUP)
    mmax, qmax = 10, 132
    triv_q, triv_m = [], []
    for c in datos["celdas"]:
        m, q = c["m"], c["q"]
        if m > mmax or q > qmax:
            continue
        sop = [d for d in c["primos"] if d["a"] > 0]
        if not sop:
            triv_q.append(q)
            triv_m.append(m)
            continue
        w = 0.9 / len(sop)
        for k, d in enumerate(sop):
            alto = d["a"] * np.log10(d["p"])
            col = FAM.get(d["familia"], MUTED)
            # shade=False: el sombreado de matplotlib oscurece el color y deja de coincidir con la leyenda
            ax.bar3d(q - 0.45 + k * w, m - 0.25, 0, w, 0.5, alto, color=col, edgecolor=INK2, linewidth=0.25,
                     shade=False)
    ax.scatter(triv_q, triv_m, np.zeros(len(triv_q)), s=3, color=BASE, depthshade=False)
    az = [(r["q"], r["m"]) for r in datos["A_igual_Z"] if r["m"] <= mmax and r["q"] <= qmax]
    ax.scatter([a[0] for a in az], [a[1] for a in az], [0] * len(az), marker="x", s=14, color=INK, depthshade=False)
    ax.set_xlabel("q", labelpad=2)
    ax.set_ylabel("m", labelpad=2)
    ax.set_zlabel(r"$\log_{10}[\mathcal{O}:A_m(q)]$ by prime", labelpad=4)
    ax.set_xlim(0, qmax)
    ax.set_ylim(1.5, mmax + 0.5)
    ax.set_yticks(range(2, mmax + 1))
    ax.view_init(elev=24, azim=-62)
    for eje in (ax.xaxis, ax.yaxis, ax.zaxis):
        eje.pane.set_facecolor(SUP)
        eje.pane.set_edgecolor(GRID)
        eje._axinfo["grid"]["color"] = GRID
    leyenda = [Patch(color=C_H, label=r"$q'\mid h$"), Patch(color=C_H1, label=r"$q'\mid h+1$"),
               Patch(color=C_H2, label=r"$q'\mid h+2$"),
               plt.Line2D([], [], marker="o", ls="", color=BASE, ms=3, label=r"index $1$ ($A=\mathcal{O}$)"),
               plt.Line2D([], [], marker="x", ls="", color=INK, ms=5, label=r"$A=\mathbb{Z}$")]
    ax.legend(handles=leyenda, loc="upper left", fontsize=7.5, bbox_to_anchor=(0.0, 0.98))
    guardar(fig, "fig_paisaje", prev)


# ---------------------------------------------------------------- paisaje 2D (paneles)
def fig_paisaje2d(prev):
    """Un panel por m; eje x = q; barra apilada por primo, altura v_p log10 p, color = familia."""
    datos = json.load(open(os.path.join(AQUI, "datos_indice.json")))
    ms = list(range(2, 11))
    qmax = 132
    fig, axs = plt.subplots(len(ms), 1, figsize=(7.2, 8.6), sharex=True)
    ymax = 0
    por_m = {m: [] for m in ms}
    for c in datos["celdas"]:
        if c["m"] in por_m and c["q"] <= qmax:
            por_m[c["m"]].append(c)
            ymax = max(ymax, sum(d["a"] * np.log10(d["p"]) for d in c["primos"]))
    for ax, m in zip(axs, ms):
        calculadas = sorted(c["q"] for c in por_m[m])
        ax.scatter(calculadas, [-0.12] * len(calculadas), s=2, color=BASE, zorder=1)
        for c in por_m[m]:
            base = 0.0
            for d in sorted((d for d in c["primos"] if d["a"] > 0), key=lambda d: d["p"]):
                alto = d["a"] * np.log10(d["p"])
                ax.bar(c["q"], alto, bottom=base, width=0.9, color=FAM.get(d["familia"], MUTED),
                       edgecolor=SUP, linewidth=0.6, zorder=2)
                base += alto
        for q in (2 * m + 1, 2 * m + 2):
            if q <= qmax:
                ax.plot(q, -0.12, marker="x", ms=4, color=INK, zorder=3)
        ax.set_ylim(-0.3, ymax * 1.05)
        ax.set_yticks([0, 2, 4])
        ax.tick_params(axis="y", labelsize=6.5)
        ax.text(1.005, 0.5, "$m=%d$" % m, transform=ax.transAxes, fontsize=8, color=INK2, va="center")
        for s in ("right", "top"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
    axs[-1].set_xlabel("$q$")
    axs[-1].set_xlim(2, qmax + 1)
    fig.text(0.02, 0.5, T(r"$\log_{10}[\mathcal{O}:A_m(q)]$, split by prime", r"$\log_{10}[\mathcal{O}:A_m(q)]$, por primos"), rotation=90, va="center",
             fontsize=9, color=INK2)
    leyenda = [Patch(color=C_H, label=r"$q'\mid h$"), Patch(color=C_H1, label=r"$q'\mid h+1$"),
               Patch(color=C_H2, label=r"$q'\mid h+2$"), Patch(color=MUTED, label=r"$q'\leq 2$"),
               plt.Line2D([], [], marker="o", ls="", color=BASE, ms=3, label=T("computed", "calculado")),
               plt.Line2D([], [], marker="x", ls="", color=INK, ms=5, label=r"$A=\mathbb{Z}$")]
    axs[0].legend(handles=leyenda, ncol=6, loc="lower center", bbox_to_anchor=(0.5, 1.05), fontsize=7)
    fig.subplots_adjust(hspace=0.25, left=0.08, right=0.94)
    guardar(fig, "fig_paisaje2d", prev)


# ---------------------------------------------------------------- perfil 3D
def fig_perfil(prev):
    from perfil_W import dims_B
    from reparto_p2 import phi
    casos = [(3, 5, 2, 7), (3, 7, 2, 10), (3, 11, 2, 16), (3, 13, 2, 19),
             (5, 7, 2, 17), (5, 9, 2, 22), (5, 11, 2, 27), (5, 13, 2, 32),
             (7, 3, 2, 10), (7, 5, 2, 17)]
    K = 9
    # 2D y no 3D: en 3D las etiquetas se amontonaban y el escalon p-1 quedaba tapado.  Celda = W_i en
    # su valor absoluto; color = fraccion de la capa llena; borde naranja = nivel p-1.
    fig, ax = plt.subplots(figsize=(5.6, 3.5))
    etiquetas = []
    M = np.full((len(casos), K), np.nan)
    Wabs = [[None] * K for _ in casos]
    for fila, (p, n, v, m) in enumerate(casos):
        q = p ** v * n
        Kc = min(K, phi(p ** v))
        dB = dims_B(q, p, m, n, Kc)
        W = [dB[0]] + [dB[i] - dB[i - 1] for i in range(1, Kc)]
        full = phi(n) // 2
        etiquetas.append(r"$p=%d$,  $q'=%d$   ($m=%d$, $q=%d$)" % (p, n, m, q))
        for i in range(1, Kc):
            M[fila, i] = W[i] / full
            Wabs[fila][i] = "%d/%d" % (W[i], full)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", [SUP, SEQ[1], SEQ[3], SEQ[4]])
    ax.imshow(np.where(np.isnan(M), 0, M), cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for fila, (p, n, v, m) in enumerate(casos):
        for i in range(1, K):
            if Wabs[fila][i] is None:
                continue
            txt = Wabs[fila][i]
            ax.text(i, fila, txt, ha="center", va="center", fontsize=7,
                    color=SUP if M[fila, i] >= 0.6 else INK2)
        ax.add_patch(matplotlib.patches.Rectangle((p - 1 - 0.5, fila - 0.5), 1, 1, fill=False,
                                                  edgecolor=C_H1, lw=2))
    ax.set_xticks(range(1, K))
    ax.set_xlim(0.5, K - 0.5)
    ax.set_xlabel(T("layer $i$", "capa $i$"))
    ax.set_yticks(range(len(casos)))
    ax.set_yticklabels(etiquetas, fontsize=8)
    ax.set_xticks(np.arange(0.5, K, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(casos), 1), minor=True)
    ax.grid(which="minor", color=GRID, lw=0.6)
    ax.tick_params(which="minor", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.legend(handles=[Patch(facecolor="none", edgecolor=C_H1, lw=2, label=T(r"layer $i=p-1$", r"capa $i=p-1$"))],
              loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=7.5)
    guardar(fig, "fig_perfil", prev)


# ---------------------------------------------------------------- clase
def fig_clase(prev):
    from bernoulli_W1 import W1, h_menos_numerico, es_primo, factores
    qs = [x for x in range(5, 90) if es_primo(x)]
    ps = [3, 5, 7, 11, 13]
    M = np.zeros((len(qs), len(ps)))
    etiquetas = []
    for i, n in enumerate(qs):
        h = h_menos_numerico(n)
        fac = factores(h)
        txt = "1" if h == 1 else r"\cdot ".join(
            ("%d^{%d}" % (f, fac.count(f)) if fac.count(f) > 1 else "%d" % f) for f in sorted(set(fac)))
        etiquetas.append((n, h, txt))
        for j, p in enumerate(ps):
            if p == n:
                M[i, j] = np.nan
                continue
            M[i, j] = (n - 1) // 2 - W1(n, p)
    fig, ax = plt.subplots(figsize=(3.9, 5.6))
    cmap = matplotlib.colors.ListedColormap([SUP, SEQ[2], SEQ[4]])
    ax.imshow(np.nan_to_num(M, nan=0), cmap=cmap, vmin=-0.5, vmax=2.5, aspect="auto")
    for i in range(len(qs)):
        for j in range(len(ps)):
            if np.isnan(M[i, j]):
                ax.text(j, i, "—", ha="center", va="center", color=MUTED, fontsize=7)
            elif M[i, j] > 0:
                ax.text(j, i, "%d" % M[i, j], ha="center", va="center", color=SUP, fontsize=8, weight="bold")
    ax.set_xticks(range(len(ps)))
    ax.set_xticklabels(["p=%d" % p for p in ps])
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(qs)))
    ax.set_yticklabels(["$q'=%d$    $h^-=%s$" % (n, t) for (n, h, t) in etiquetas], fontsize=8)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(ps), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(qs), 1), minor=True)
    ax.grid(which="minor", color=GRID, lw=0.6)
    ax.tick_params(which="minor", length=0)
    guardar(fig, "fig_clase", prev)
    return M, qs, ps


# ---------------------------------------------------------------- el hilo de la nota
def fig_estructura(prev):
    """Diagrama del recorrido: cada seccion responde a la pregunta que deja la anterior; el color marca la linea
    historica de la seccion 2 en la que se apoya; la linea discontinua es el ejemplo guia A_3(35)."""
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    COL_CAR, COL_CON, COL_BER = C_H, C_H2, C_H1  # valores de caracteres, conductores, Bernoulli

    def caja(x, y, w, h, titulo, texto, col):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    facecolor=SUP, edgecolor=BASE, linewidth=0.8, zorder=2))
        ax.add_patch(FancyBboxPatch((x, y), 1.2, h, boxstyle="square,pad=0", facecolor=col,
                                    edgecolor="none", zorder=3))
        ax.text(x + 3, y + h - 2.2, titulo, ha="left", va="top", fontsize=8, color=INK, weight="bold", zorder=4)
        ax.text(x + 3, y + 2.0, texto, ha="left", va="bottom", fontsize=7, color=INK2, zorder=4)

    def flecha(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9, shrinkA=0, shrinkB=0), zorder=1)

    W, H, X = 58, 11.5, 21
    filas = [
        (86, T(r"§3  When is the ring just $\mathbb{Z}$?", r"§3  ¿Cuándo es el anillo $\mathbb{Z}$?"),
         T(r"only at $q=2m+1,\,2m+2$ (and $(6,1)$)", r"solo en $q=2m+1,\,2m+2$ (y $(6,1)$)"), COL_CAR),
        (71, T("§4  Which primes divide the index?", "§4  ¿Qué primos dividen al índice?"),
         T(r"$p\mid q$ with $q'\mid h,\,h+1,\,h+2$", r"$p\mid q$ con $q'\mid h,\,h+1,\,h+2$"), COL_CAR),
        (56, T("§5  What does the defect look like at p?", "§5  ¿Cómo es el defecto en p?"),
         T(r"layers $W_i$: lengths, exponent, Gorenstein, type", r"capas $W_i$: longitudes, exponente, Gorenstein, tipo"), COL_CON),
    ]
    for y, tit, txt, col in filas:
        caja(X, y, W, H, tit, txt, col)
    for (y0, *_), (y1, *_) in zip(filas, filas[1:]):
        flecha(X + W / 2, y0, X + W / 2, y1 + H + 0.8)
    # seccion 6 con dos ramas
    caja(X, 41, W, H, T("§6  Why do class numbers appear?", "§6  ¿Por qué aparecen números de clase?"),
         T("the moment matrix of the first layer", "la matriz de momentos de la primera capa"), COL_BER)
    flecha(X + W / 2, 56, X + W / 2, 41 + H + 0.8)
    wr = 38
    caja(2, 23.5, wr, H, T("6.1  sawtooth", "6.1  diente de sierra"),
         T(r"Sinnott's index, doubling, $h^-$", r"índice de Sinnott, duplicación, $h^-$"), COL_BER)
    caja(60, 23.5, wr, H, T("6.2  shifted sawtooth", "6.2  diente desplazado"),
         T(r"factor at 2: is $-1\in\langle 2\rangle$?", r"factor en 2: ¿$-1\in\langle 2\rangle$?"), COL_BER)
    flecha(X + 14, 41, 2 + wr / 2, 23.5 + H + 0.8)
    flecha(X + W - 14, 41, 60 + wr / 2, 23.5 + H + 0.8)
    caja(X, 6, W, H, T("§7  What lies beyond the first layer?", "§7  ¿Qué hay más allá de la primera capa?"),
         T(r"universal ladder below $p^a$, descent at $p^a$; certified exponents", r"escalera bajo $p^a$, descenso en $p^a$; exponentes certificados"), COL_CON)
    flecha(2 + wr / 2, 23.5, X + 14, 6 + H + 0.8)
    flecha(60 + wr / 2, 23.5, X + W - 14, 6 + H + 0.8)
    ax.text(X + W / 2, 1.0, T("§8 conclusions and limits  ·  §9 how the numbers were computed",
                              "§8 conclusiones y límites  ·  §9 cómo se calcularon los números"),
            ha="center", va="bottom", fontsize=7, color=MUTED)
    # ejemplo guia A_3(35): de la seccion 4 a la 6.2
    xs = [X + W + 3.5, X + W + 3.5, X + W + 3.5, 60 + wr - 3]
    ys = [71 + H / 2, 56 + H / 2, 41 + H / 2, 23.5 + H]
    ax.plot(xs, ys, ls=(0, (3, 2)), color=INK2, lw=0.9, zorder=1)
    for x_, y_ in zip(xs[:3], ys[:3]):
        ax.plot(x_, y_, "o", ms=3.5, color=INK2, zorder=2)
    ax.text(X + W + 5, 71 + H / 2 + 3.5, T("running example\n$A_3(35)$", "ejemplo guía\n$A_3(35)$"), fontsize=7,
            color=INK2, ha="left", va="bottom")
    leyenda = [Patch(color=COL_CAR, label=T("character values", "valores de caracteres")),
               Patch(color=COL_CON, label=T("conductors of one-dimensional rings", "conductores de anillos unidimensionales")),
               Patch(color=COL_BER, label=T("Bernoulli and class numbers", "Bernoulli y números de clase"))]
    ax.legend(handles=leyenda, loc="upper center", bbox_to_anchor=(0.5, 1.06), ncol=3, fontsize=7,
              handlelength=1.2, columnspacing=1.2)
    guardar(fig, "fig_estructura", prev)


# ---------------------------------------------------------------- simetria de Gorenstein
def fig_simetria(prev):
    """Proposiciones 7 y 9: columna i = W_i (abajo) + W_{e-1-i} (arriba); Gorenstein <=> todas llegan a d.
    Datos: filas con reticulo de perfil_W_OUT.txt (W, a, b, forma (f, e))."""
    import re
    fila = re.compile(r"^\s+ok\s+m=\s*(\d+) q=\s*(\d+) p=(\d+) v=\d+ q'=\s*(\d+) E=\s*(\d+) forma=(\[.*?\]) "
                      r"a=(\d+) b=(\d+)\s+W=(\[.*?\])")
    ruta = os.path.join(AQUI, "perfil_W_OUT.txt")
    if not os.path.exists(ruta):
        ruta = os.path.join(AQUI, "..", "gates", "perfil_W_OUT.txt")
    casos = [(3, 35, 5), (7, 45, 3), (10, 63, 3), (6, 56, 2), (9, 40, 2), (16, 99, 3)]
    datos = {}
    for ln in open(ruta, encoding="utf-8"):
        mt = fila.match(ln)
        if mt and (int(mt.group(1)), int(mt.group(2)), int(mt.group(3))) in casos:
            from reparto_p2 import phi
            n = int(mt.group(4))
            datos[(int(mt.group(1)), int(mt.group(2)), int(mt.group(3)))] = (
                eval(mt.group(9)), eval(mt.group(6))[0][1], phi(n) // 2, int(mt.group(7)), int(mt.group(8)))
    assert len(datos) == len(casos), sorted(datos)
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.3))
    for ax, clave in zip(axs.flat, casos):
        W, e, d, a, b = datos[clave]
        simetrico = all(W[i] + W[e - 1 - i] == d for i in range(e))
        assert simetrico == (a == b), (clave, W, e, d, a, b)  # Proposicion 9 contra las longitudes del reticulo
        xs = range(e)
        abajo = [W[i] for i in xs]
        arriba = [W[e - 1 - i] for i in xs]
        ax.bar(xs, abajo, width=0.62, color=C_H, edgecolor=SUP, linewidth=1.2, zorder=2)
        ax.bar(xs, arriba, bottom=abajo, width=0.62, color=C_H1, edgecolor=SUP, linewidth=1.2, zorder=2)
        ax.axhline(d, color=INK, lw=1, ls=(0, (4, 3)), zorder=3)
        ax.set_ylim(0, max(d, max(x + y for x, y in zip(abajo, arriba))) + 0.9)
        ax.set_xticks(list(xs))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.tick_params(labelsize=7)
        m, q, p = clave
        estado = T("Gorenstein", "Gorenstein") if simetrico else T("not Gorenstein", "no Gorenstein")
        ax.set_title(r"$A_{%d}(%d)$, $p=%d$:  $d=%d$, $e=%d$" % (m, q, p, d, e) + "\n" + estado,
                     fontsize=8, color=INK if simetrico else INK2)
        for s in ("right", "top"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
    for ax in axs[1]:
        ax.set_xlabel(T("layer $i<e$", "capa $i<e$"), fontsize=8)
    leyenda = [Patch(color=C_H, label=r"$W_i$"), Patch(color=C_H1, label=r"$W_{e-1-i}$"),
               plt.Line2D([], [], color=INK, lw=1, ls=(0, (4, 3)), label=r"$d$")]
    fig.legend(handles=leyenda, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.03), fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    guardar(fig, "fig_simetria", prev)


# ---------------------------------------------------------------- Stickelberger, q' compuesto
def fig_stickel(prev):
    """Proposicion de Sinnott en la matriz de momentos: defecto d - rango_Fp M_n para n compuesto (o potencia de primo),
    con h^- exacto (stickelberger_smith_OUT.txt).  Filas con daga: la matriz de solo unidades pierde rango sobre Q."""
    from math import gcd
    from fractions import Fraction

    def rango_Q(filas):
        A = [[Fraction(x) for x in f] for f in filas]
        r = 0
        for c in range(len(A[0]) if A else 0):
            piv = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
            if piv is None:
                continue
            A[r], A[piv] = A[piv], A[r]
            for i in range(len(A)):
                if i != r and A[i][c] != 0:
                    f = A[i][c] / A[r][c]
                    A[i] = [x - f * y for x, y in zip(A[i], A[r])]
            r += 1
        return r
    from bernoulli_W1 import rango_mod_p, factores
    ruta = os.path.join(AQUI, "stickelberger_smith_OUT.txt")
    if not os.path.exists(ruta):
        ruta = os.path.join(AQUI, "..", "gates", "stickelberger_smith_OUT.txt")
    hm = {}
    for ln in open(ruta, encoding="utf-8"):
        cs = ln.split()
        if len(cs) > 5 and cs[0].isdigit() and "Delta/h^-" in ln:
            hm[int(cs[0])] = int(cs[4])

    def primo(x):
        return x > 1 and all(x % k for k in range(2, int(x ** 0.5) + 1))

    ns = [n for n in range(8, 90) if n % 4 != 2 and not primo(n) and n in hm]
    ps = [3, 5, 7, 11, 13, 17]
    M = np.zeros((len(ns), len(ps)))
    daga = []
    for i, n in enumerate(ns):
        reps = [u for u in range(1, n) if gcd(u, n) == 1 and u < n - u]
        d = len(reps)
        s = lambda x: 0 if x % n == 0 else 2 * (x % n) - n
        filas = [[s(pow(u, -1, n) * c) for c in range(1, n)] for u in reps]
        unidades = [[s(pow(u, -1, n) * c) for c in range(1, n) if gcd(c, n) == 1] for u in reps]
        daga.append(rango_Q(unidades) < d)
        for j, p in enumerate(ps):
            M[i, j] = np.nan if n % p == 0 else d - rango_mod_p(filas, p)
            if not np.isnan(M[i, j]):
                assert (M[i, j] > 0) == (hm[n] % p == 0), (n, p)  # la proposicion, en la figura
    fig, ax = plt.subplots(figsize=(4.4, 7.4))
    cmap = matplotlib.colors.ListedColormap([SUP, SEQ[2], SEQ[4]])
    ax.imshow(np.nan_to_num(M, nan=0), cmap=cmap, vmin=-0.5, vmax=2.5, aspect="auto")
    for i in range(len(ns)):
        for j in range(len(ps)):
            if np.isnan(M[i, j]):
                ax.text(j, i, "—", ha="center", va="center", color=MUTED, fontsize=7)
            elif M[i, j] > 0:
                ax.text(j, i, "%d" % M[i, j], ha="center", va="center", color=SUP, fontsize=8, weight="bold")
    ax.set_xticks(range(len(ps)))
    ax.set_xticklabels(["p=%d" % p for p in ps])
    ax.xaxis.tick_top()
    etq = []
    for n, dg in zip(ns, daga):
        fac = factores(hm[n])
        txt = "1" if hm[n] == 1 else r"\cdot ".join(
            ("%d^{%d}" % (f, fac.count(f)) if fac.count(f) > 1 else "%d" % f) for f in sorted(set(fac)))
        etq.append(("$\\dagger$ " if dg else "") + "$q'=%d$    $h^-=%s$" % (n, txt))
    ax.set_yticks(range(len(ns)))
    ax.set_yticklabels(etq, fontsize=7.5)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(ps), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ns), 1), minor=True)
    ax.grid(which="minor", color=GRID, lw=0.6)
    ax.tick_params(which="minor", length=0)
    guardar(fig, "fig_stickel", prev)
    return sum(daga), len(ns)


# ---------------------------------------------------------------- sierra
# ---------------------------------------------------------------- zonas de la primera capa
def fig_zonas(prev):
    """Donde cae W_1/d en los perfiles medidos y que enunciado decide el exponente en cada banda."""
    import re
    from reparto_p2 import phi
    ruta = os.path.join(AQUI, "perfil_W_OUT.txt")
    for alt in ("anc", os.path.join("..", "gates")):
        if not os.path.exists(ruta):
            ruta = os.path.join(AQUI, alt, "perfil_W_OUT.txt")
    pat = re.compile(r"m=\s*(\d+) q=\s*(\d+) p=(\d+) v=(\d+) q'=\s*(\d+).*?a=(\d+) b=(\d+)\s+W=\[([0-9, ]*)\]")
    pts = []
    for ln in open(ruta, encoding="utf-8", errors="replace"):
        mm = pat.search(ln)
        if not mm:
            continue
        m, q, p, v, n, a, b = [int(x) for x in mm.groups()[:7]]
        W = [int(x) for x in mm.group(8).split(",") if x.strip()]
        if n in (1, 2, 3, 4, 6) or len(W) < 2:
            continue
        d = phi(n) // 2
        pts.append((W[1] / d, d, p, a, b))
    bandas = [(-0.035, 0.02, "(a)"), (0.02, 0.5, "(b)"), (0.5, 0.985, "(c)"), (0.985, 1.035, "(d)")]
    cuenta = [0] * len(bandas)
    for (w, d, p, a, b) in pts:
        k = 0 if w == 0 else (1 if 2 * w <= 1 else (2 if w < 1 else 3))
        cuenta[k] += 1

    fig, ax = plt.subplots(figsize=(6.0, 2.7))
    dmax = max(d for (w, d, p, a, b) in pts)
    ax.set_ylim(0.5, dmax + 3.0)
    for k, (x0, x1, lab) in enumerate(bandas):
        ax.axvspan(x0, x1, ymax=(dmax + 1.6 - 0.5) / (dmax + 2.5), color=(SUP if k % 2 else BANDA),
                   lw=0, zorder=0)
        ax.plot([x0, x0], [0.5, dmax + 2.6], color=GRID, lw=0.7, zorder=1)
        ax.text((x0 + x1) / 2, dmax + 2.1, "%s  %d" % (lab, cuenta[k]), ha="center", va="center",
                fontsize=8, color=INK)
    ax.plot([-0.035, 1.035], [dmax + 1.6, dmax + 1.6], color=GRID, lw=0.7)
    rng = np.random.default_rng(7)
    for gor, col, lab in ((True, C_H, T("Gorenstein", "de Gorenstein")),
                          (False, C_H1, T("not Gorenstein", "no de Gorenstein"))):
        sel = [(w, d) for (w, d, p, a, b) in pts if (a == b) == gor]
        ax.scatter([w + rng.uniform(-0.006, 0.006) for (w, d) in sel],
                   [d + rng.uniform(-0.22, 0.22) for (w, d) in sel],
                   s=15, color=col, edgecolor=SUP, linewidth=0.4, zorder=3, label=lab)
    ax.set_xlabel(T("first layer $W_1/d$", "primera capa $W_1/d$"))
    ax.set_ylabel(r"$d=\varphi(q')/2$")
    ax.set_xlim(-0.035, 1.035)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([2, 6, 10, 14])
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    ax.legend(loc="lower left", fontsize=7.5, bbox_to_anchor=(0.03, 0.02))
    guardar(fig, "fig_zonas", prev)
    return cuenta


def fig_sierra(prev):
    n = 13
    c = np.arange(1, n)
    s_h = 2 * c - n
    s_h1 = np.where(c <= n // 2, c, c - n)
    fig, ax = plt.subplots(figsize=(5.4, 2.5))
    ax.axhline(0, color=BASE, lw=1)
    ax.plot(c, s_h, "-o", color=C_H, lw=2, ms=5, mec=SUP, mew=1.5)
    ax.plot(c, s_h1, "-o", color=C_H1, lw=2, ms=5, mec=SUP, mew=1.5)
    ax.text(12.3, s_h[-1], T(r"$q'\mid h$ or $h+2$:  $S=N(2c-q')$", r"$q'\mid h$ o $h+2$:  $S=N(2c-q')$"), color=INK, fontsize=7.5, va="center")
    ax.text(12.3, s_h1[-1] - 1.5, T(r"$q'\mid h+1$:  $S=N\hat c$ (shifted by $q'/2$)", r"$q'\mid h+1$:  $S=N\hat c$ (desplazada $q'/2$)"), color=INK, fontsize=7.5, va="center")
    ax.set_xlabel(T("residue $c$ modulo $q'=13$", "residuo $c$ módulo $q'=13$"))
    ax.set_ylabel("$S(c)/N$")
    ax.set_xticks(c)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.6)
    guardar(fig, "fig_sierra", prev)


if __name__ == "__main__":
    prev = sys.argv[1] if len(sys.argv) > 1 else None
    fig_paisaje(prev)  # 3D, ya no la usa la nota; se conserva en ingles
    for IDIOMA in ("en", "es"):
        fig_historia(prev)
        fig_sierra(prev)
        fig_paisaje2d(prev)
        M, qs, ps = fig_clase(prev)
        fig_perfil(prev)
        fig_simetria(prev)
        fig_stickel(prev)
        cuenta = fig_zonas(prev)
        fig_estructura(prev)
    print("fig_clase: celdas con deficit > 0:", [(qs[i], ps[j], int(M[i, j])) for i in range(len(qs))
                                                 for j in range(len(ps)) if not np.isnan(M[i, j]) and M[i, j] > 0])
    print("fig_zonas: reparto por bandas [W_1=0, 0<W_1<=d/2, d/2<W_1<d, W_1=d]:", cuenta)
    print("ok")
