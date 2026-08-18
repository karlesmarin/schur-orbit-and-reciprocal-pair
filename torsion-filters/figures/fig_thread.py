"""
Figure: the paper as one thread that forks, and the pearls it is strung from -- borrowed ones included.

Authors: Carles Marin, Claude (AI assistant).

The companion paper opens with a thread: each question answered by the pearl below it, and each
question picking up a word from the previous answer. This is the same device, with the one change the
subject forces on it -- the thread FORKS. Halfway along, the question "what does the specialisation
do?" is answered by a regularity condition, and from there the paper has to be walked twice, once for
each parity of t, because the two walks are through different groups. They rejoin at the fusion
quotient -- which is where the two filters turn out to be the same kind of object -- and end on the
same open pearl.

The fork is not decoration. It is the paper's thesis in one line: the two strands differ because the
evaluation point sits in different components of the orthogonal group, and because 2 is invertible
modulo an odd t and not modulo an even one.

Some pearls are not ours, and the figure says so. The mechanism by which a weight on a wall dies and
a regular weight folds with a sign is standard representation theory at roots of unity; the odd
filter is a corollary of the principal-element theorem of Nadimpalli-Pattanayak-Prasad; the branching
to a maximal-rank subgroup is Koike-Terada. So the pearls carry the colour of their status, and a
reader can see at a glance which beads were strung here and which were handed to us.

Both strands end on the same open pearl. That is not a rhetorical flourish: the factorisation is
complete on both sides and the last coefficient is the part none of it explains.

But the two parities do not stop at the same distance from it, and since the GKRS reading the gap is
large enough to draw. So the thread grows a short tail on the odd side: there the numerator of the
extremal coefficient is a signed transversal count and is proved, and the only thing left is one
division. The general pearl stays open; the odd tail says how much of it is not.

Palette: the three house status colours, unchanged from the companion (blue = proved, orange =
verified by computation, warm grey = external), plus an unfilled outline for what is open. Shape
doubles the colour so the figure survives greyscale printing.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"
PROVED = "#2a78d6"; VERIF = "#eb6834"; EXT = "#6b6560"; THREAD = "#b9b7ae"

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

# (x, y, pregunta, perla, estado, seccion)
TRONCO_IZQ = [
    (-0.60, 0.0, "when does\nit vanish?",
     "reduce to $t\\in\\{1,2\\}$\nplus a specialisation", "proved", "§2"),
    (1.45, 0.0, "what does the\nspecialisation do?",
     "it is a regularity\ncondition", "proved", "§3"),
]
RAMA_IMPAR = [
    (2.95, 1.35, "regular in which\ngroup?  ($t$ odd)",
     "$B_{R'}\\!\\downarrow B_{m'}\\!\\times\\! D_r$:\nordinary, genuine", "verif", "§5"),
    (5.20, 1.35, "and its filter?", "principal element:\ncited", "ext", "§5"),
]
RAMA_PAR = [
    (3.10, -1.35, "regular in which\ngroup?  ($t$ even)",
     "twining: the input\ngoes virtual", "proved", "§4"),
    (4.85, -1.35, "and its filter?", "$\\tau^C_t$ with its sign,\nproved here", "proved", "§3"),
]
TRONCO_DER = [
    (6.60, 0.0, "where do the two\nfilters live?",
     "minimal fusion: one\ntensor-sector point", "ext", "§6"),
    (8.55, 0.0, "one point: which\nweight is on top?",
     "$\\mathrm{top}\\,\\mathrm{Newt}(N_\\beta)$\n$-\\,\\mathrm{top}\\,\\mathrm{Newt}(N_\\delta)$", "verif", "§7"),
    (10.50, 0.0, "how big is that\ncoefficient?", "primitive:\n$q_t[M]=\\pm[0]$", "open", "§8"),
]
# la cola impar: la pregunta general sigue abierta, pero de este lado ya no lo esta entera.
COLA_IMPAR = [
    (10.10, -3.35, "and for $t$ odd?",
     "$\\nu$: a signed\ntransversal count", "proved", "§5"),
    (7.10, -3.35, "and the division\nby $\\Delta_t$?",
     "a sum along a\nprogression of step $2t$", "proved", "§5"),
    (4.10, -3.35, "and that sum?",
     # El +- NO es adorno: lo PROBADO es c = kappa_{t,r} eps_t det M con kappa en {+-1}; que
     # kappa valga +1 es una MEDICION en doce configuraciones.  Sin el signo, el nodo azul
     # --- "proved" --- afirmaba como teorema lo que es una medida.  La Figura 13 (fig_problems)
     # ya lo escribia bien.
     "$c=\\pm\\epsilon_t\\det M$,\nan explicit $0/{\\pm}1$ matrix", "proved", "§5"),
    (1.10, -3.35, "so what is left?",
     "$M$ unimodular;\n(L2),(L3) still open", "open", "§5"),
]
COL = {"proved": PROVED, "verif": VERIF, "ext": EXT, "open": INK}

# ---------------------------------------------------------------------------------------------
# Los numeros de seccion de las perlas se COMPRUEBAN contra el .aux, no se escriben de memoria.
# Esta figura ya se quedo una vez con "§6 top weight" despues de insertar una seccion nueva, y el
# pie de una figura es justo donde nadie mira.  Si el .aux no esta, se avisa y no se falla.
ETIQUETA = {"§2": "sec:red", "§3": "sec:filter", "§4": "sec:comp", "§5": "sec:odd",
            "§6": "sec:fusion", "§7": "sec:top", "§8": "sec:unit"}


def comprueba_secciones(perlas):
    import os
    import re
    if not os.path.exists("orbit_pair_ii.aux"):
        print("   aviso: no hay .aux, no se comprueban los numeros de seccion")
        return
    aux = open("orbit_pair_ii.aux", encoding="utf-8", errors="replace").read()
    real = {}
    for lab, num in re.findall(r"\\newlabel\{(sec:[a-z0-9]+)\}\{\{([0-9.]+)\}", aux):
        real[lab] = num
    malas = []
    for (_, _, _, _, _, sec) in perlas:
        lab = ETIQUETA.get(sec)
        if lab and lab in real and real[lab] != sec.lstrip("§"):
            malas.append((sec, lab, real[lab]))
    if malas:
        raise AssertionError("la figura cita secciones que el .aux desmiente: %s" % malas)
    print("   numeros de seccion comprobados contra el .aux: %d perlas, 0 discrepancias" % len(perlas))


fig, ax = plt.subplots(figsize=(14.6, 7.6))


def hilo(p, q, curva=0.0):
    """une dos perlas con una curva suave."""
    xs = np.linspace(p[0], q[0], 120)
    s = (xs - p[0]) / max(q[0] - p[0], 1e-9)
    ys = p[1] + (q[1] - p[1]) * (3 * s ** 2 - 2 * s ** 3) + curva * np.sin(np.pi * s)
    ax.plot(xs, ys, color=THREAD, linewidth=1.5, zorder=1, solid_capstyle="round")


hilo(TRONCO_IZQ[0], TRONCO_IZQ[1], 0.05)
hilo(TRONCO_IZQ[1], RAMA_IMPAR[0])
hilo(TRONCO_IZQ[1], RAMA_PAR[0])
hilo(RAMA_IMPAR[0], RAMA_IMPAR[1], 0.04)
hilo(RAMA_PAR[0], RAMA_PAR[1], -0.04)
hilo(RAMA_IMPAR[1], TRONCO_DER[0])
hilo(RAMA_PAR[1], TRONCO_DER[0])
hilo(TRONCO_DER[0], TRONCO_DER[1], 0.05)
hilo(TRONCO_DER[1], TRONCO_DER[2], 0.05)
# el retorno de renglon: la cola no cabe a lo ancho sin encoger la figura entera, asi que baja por
# la derecha y vuelve a la izquierda, como un renglon que continua debajo.
_rx = np.array([10.50, 11.10, 11.35, 11.20, 10.60, 10.10])
_ry = np.array([-0.16, -0.60, -1.40, -2.40, -3.10, -3.35])
_t = np.linspace(0, 1, len(_rx))
_tt = np.linspace(0, 1, 500)
ax.plot(np.interp(_tt, _t, _rx), np.interp(_tt, _t, _ry), color=THREAD, linewidth=1.5,
        zorder=1, solid_capstyle="round")
for _i in range(len(COLA_IMPAR) - 1):
    hilo(COLA_IMPAR[_i], COLA_IMPAR[_i + 1], 0.05)

# Las dos etiquetas de rama NO se colocan a ojo.  Lo intente cuatro veces y el control las rechazo
# cuatro veces: se prueban candidatos a lo largo de la bifurcacion y se queda el primero libre.
RAMAS = []
_CAND_ODD = [(2.42, 0.60), (2.05, 0.45), (1.72, 0.30), (2.60, 0.40), (1.50, 0.55), (2.20, 0.28)]
_CAND_EVEN = [(2.42, -0.62), (2.05, -0.47), (1.72, -0.32), (2.60, -0.42), (1.50, -0.57),
              (2.20, -0.30)]


def coloca_rama(texto, candidatos, giro):
    """Devuelve el Text ya puesto en el primer sitio que no pisa nada de TEXTOS."""
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    puestas = [t_.get_window_extent(renderer=ren) for t_ in TEXTOS]
    for (cx, cy) in candidatos:
        t_ = ax.text(cx, cy, texto, ha="center", va="center", fontsize=11.0, color=SECOND,
                     style="italic", rotation=giro)
        caja = t_.get_window_extent(renderer=ren)
        if not any(caja.x0 < b_.x1 and b_.x0 < caja.x1 and caja.y0 < b_.y1 and b_.y0 < caja.y1
                   for b_ in puestas):
            return t_
        t_.remove()
    raise AssertionError("no hay sitio para la etiqueta de rama %r" % texto)

TODAS = TRONCO_IZQ + RAMA_IMPAR + RAMA_PAR + TRONCO_DER + COLA_IMPAR
TEXTOS = []
for (x, y, preg, perla, estado, sec) in TODAS:
    col = COL[estado]
    relleno = SURFACE if estado == "open" else col
    ax.scatter([x], [y], s=300, facecolors=relleno, edgecolors=col,
               linewidths=1.9 if estado == "open" else 0.0, zorder=4,
               marker="o" if estado != "ext" else "D")
    TEXTOS.append(ax.text(x, y + 0.24, preg, ha="center", va="bottom", fontsize=10.4,
                          color=SECOND, linespacing=1.3))
    TEXTOS.append(ax.text(x, y - 0.26, perla, ha="center", va="top", fontsize=10.8,
                          color=INK if estado != "ext" else EXT, linespacing=1.35))
    # El control media solo preguntas y perlas.  Las etiquetas de rama y los numeros de seccion
    # quedaban FUERA, y son texto igual: un punto ciego del mismo tipo que ya rompio otra figura.
    TEXTOS.append(ax.text(x + 0.215, y + 0.02, sec, ha="left", va="center", fontsize=9.2,
                          color=MUTED))

RAMAS.append(coloca_rama("$t$ odd", _CAND_ODD, 32))
RAMAS.append(coloca_rama("$t$ even", _CAND_EVEN, -32))
TEXTOS.extend(RAMAS)
ax.set_xlim(-1.85, 12.75)
ax.set_ylim(-4.45, 2.45)
ax.axis("off")

leyenda = [Line2D([], [], marker="o", color=PROVED, linestyle="none", markersize=8,
                  label="proved here"),
           Line2D([], [], marker="o", color=VERIF, linestyle="none", markersize=8,
                  label="verified by computation"),
           Line2D([], [], marker="D", color=EXT, linestyle="none", markersize=7,
                  label="not ours: cited"),
           Line2D([], [], marker="o", markerfacecolor=SURFACE, markeredgecolor=INK,
                  color=SURFACE, linestyle="none", markersize=8, markeredgewidth=1.6,
                  label="open")]
ax.legend(handles=leyenda, frameon=False, fontsize=8, ncol=4, loc="lower center",
          bbox_to_anchor=(0.5, -0.03))

def comprueba_solapes(textos):
    """CONTROL de layout, MEDIDO EN PIXELES y no estimado por longitud de cadena.
    Dos etiquetas que se pisan producen una linea ilegible del tipo
    "reduce to t in {1,2}it is a regularity" -- que es exactamente lo que esta figura tuvo durante
    varias versiones, porque un pie de figura es donde nadie mira.  Estimar el ancho contando
    caracteres no vale: "$t\\in\\{1,2\\}$" son 13 caracteres de fuente y cinco de tinta.  Asi que se
    dibuja y se leen las cajas reales."""
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    cajas = [(t.get_text().replace("\n", " "), t.get_window_extent(renderer=ren)) for t in textos]
    malas = []
    for i in range(len(cajas)):
        for j in range(i + 1, len(cajas)):
            (ta, ba), (tb, bb) = cajas[i], cajas[j]
            solape = ba.intersection(ba, bb)
            if solape is not None and solape.width > 1.0 and solape.height > 1.0:
                malas.append((ta[:34], tb[:34], round(solape.width, 1), round(solape.height, 1)))
    if malas:
        raise AssertionError("etiquetas que se pisan, en pixeles: %s" % malas)
    print("   control de solapes: %d etiquetas medidas en pixeles, 0 colisiones" % len(cajas))


comprueba_secciones(TODAS)
comprueba_solapes(TEXTOS)
fig.savefig("fig_thread.pdf", bbox_inches="tight", pad_inches=0.03)
print("fig_thread.pdf escrito, %d perlas, con bifurcacion por paridad" % len(TODAS))
for p in TODAS:
    print("   %-24s -> %-36s [%s]" % (p[2].replace("\n", " "), p[3].replace("\n", " "), p[4]))
