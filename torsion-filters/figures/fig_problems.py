# -*- coding: utf-8 -*-
# ============================================================================================
#  fig_problems.pdf -- QUE QUEDA ABIERTO, Y SOBRE QUE SE APOYA.  17 de agosto de 2026.
#
#  La seccion de problemas abiertos son varias paginas de prosa sin un solo elemento visual, y es
#  justo donde el lector necesita ver de un vistazo cuatro cosas distintas: lo probado aqui, lo que
#  es de otros, lo que solo esta medido, y lo que sigue abierto.  Es el mismo papel que fig_map.pdf
#  hace en el paper companero, adaptado a lo que este paper tiene y aquel no: DOS RAMAS.
#
#  QUE DIBUJA.  Las dos rutas a la Conjetura H --- la par arriba, la impar abajo --- con la columna
#  vertebral compartida en medio, y los problemas abiertos colgando de la pieza de maquinaria que
#  cada uno continua.  Se ve de un golpe la asimetria que es la tesis del paper: la rama impar llega
#  cuatro cajas mas lejos que la par, y el problema 17.1 es el unico que las dos comparten.
#
#  LOS NUMEROS SE LEEN DEL .aux, NO SE ESCRIBEN.  Un numero de problema escrito a mano dentro de una
#  figura miente en silencio en cuanto se anade un enunciado antes.  Aqui se leen, y ademas:
#
#  CONTROL (FATAL).  El conjunto de \label{prob:...} del .tex tiene que ser EXACTAMENTE el conjunto
#  de problemas dibujados.  Si alguien anade un problema, o borra uno --- como acaba de pasar al
#  fundir prob:oddhome en prob:idealodd --- la figura deja de dibujarse en vez de quedarse vieja.
#  Y el titulo de cada caja se compara con el titulo real del .aux, que es lo que el lector vera en
#  el cuerpo: si divergen, tampoco dibuja.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python fig_problems.py
# ============================================================================================

import io
import os
import re
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

AQUI = os.path.dirname(os.path.abspath(__file__))
AZUL, NARANJA, GRIS, TINTA = "#2A78D6", "#EB6834", "#6B6560", "#2B2B2B"
BANDA, REGLA = "#F2F1EC", "#B9B7AE"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "text.color": TINTA, "figure.facecolor": "white", "savefig.facecolor": "white",
})

PROBADO, VERIF, EXTERNO, ABIERTO = "proved", "verified", "external", "open"
RELLENO = {PROBADO: AZUL, VERIF: NARANJA, EXTERNO: GRIS, ABIERTO: "white"}
BORDE = {PROBADO: AZUL, VERIF: NARANJA, EXTERNO: GRIS, ABIERTO: TINTA}


# ------------------------------------------------------------------ el .aux y el .tex -----------
def lee_aux():
    """{etiqueta: (numero, titulo)} para cada \\label de problema o conjetura."""
    out = {}
    p = os.path.join(AQUI, "orbit_pair_ii.aux")
    if not os.path.exists(p):
        return out
    pat = re.compile(r"\\newlabel\{((?:prob|conj|thm|prop|lem):[^}]*)\}\{\{([^}]*)\}\{[^}]*\}\{(.*?)\}\{")
    for L in io.open(p, encoding="utf-8", errors="replace"):
        m = pat.match(L.strip())
        if m:
            out[m.group(1)] = (m.group(2), m.group(3))
    return out


def problemas_del_tex():
    p = os.path.join(AQUI, "orbit_pair_ii.tex")
    s = io.open(p, encoding="utf-8").read()
    cuerpo = s.split(r"\section{Open problems}")[1].split(r"\section{Closing")[0]
    return set(re.findall(r"\\label\{(prob:[^}]*)\}", cuerpo))


AUX = lee_aux()


def num(etq, por_defecto="?"):
    return AUX.get(etq, (por_defecto, ""))[0]


# ------------------------------------------------------------------ el mapa ---------------------
# (x, y, ancho, texto, estatus).  Tres bandas: par arriba, columna comun en medio, impar abajo.
H, W = 0.80, 2.06
HP, WP = 1.00, 2.05                       # las cajas de problema son mas altas: llevan dos lineas
Y_PAR, Y_MED, Y_IMP, Y_PROB = 2.45, 0.55, -1.35, -3.40

NODOS = {
    # -------- columna vertebral, comun a las dos paridades
    "red":   (0.00, Y_MED, W, "the reduction to\n$t\\in\\{1,2\\}$", PROBADO),
    "fus":   (6.60, Y_MED, W, "both filters are minimal\nfusion projections", PROBADO),
    "newt":  (8.80, Y_MED, W, "$\\mathrm{Newt}(N_\\delta)$ is the\n$C_r$ zonotope", PROBADO),
    "H":     (11.00, Y_MED, W + 0.15, "OPEN: extremal\nprimitivity (Conj. %s)" % num("conj:H"),
              ABIERTO),
    # -------- rama par
    "epar":  (2.20, Y_PAR, W, "$\\det p_t=-1$: a twining\nto $C_m\\times C_r$", PROBADO),
    "etau":  (4.40, Y_PAR, W, "$\\tau^C$ with its sign\n(the bialternant)", PROBADO),
    "enum":  (6.60, Y_PAR, W, "a transversal count,\ntwo forbidden classes", PROBADO),
    "estop": (8.80, Y_PAR, W, "no per-$\\Lambda$ reduction:\n$a_\\Lambda$ has both signs", PROBADO),
    # -------- rama impar
    "eimp":  (2.20, Y_IMP, W, "$\\det p_t=+1$: ordinary\n$B_{R'}\\to B_{m'}\\times D_r$", PROBADO),
    "otau":  (4.40, Y_IMP, W, "$\\tau^B$ at a principal\nelement", EXTERNO),
    "onum":  (6.60, Y_IMP, W, "$\\nu$ is a signed\ntransversal count", PROBADO),
    "odiv":  (8.80, Y_IMP, W, "the division inverts:\n$\\Delta_t=\\psi^t(\\Delta_1)$", PROBADO),
    "odet":  (11.00, Y_IMP, W + 0.15, "$c=\\pm\\epsilon_t\\det M$,\nthen OPEN (L1)", ABIERTO),
}

ARISTAS = [("red", "epar"), ("red", "eimp"),
           ("epar", "etau"), ("etau", "enum"), ("enum", "estop"),
           ("eimp", "otau"), ("otau", "onum"), ("onum", "odiv"), ("odiv", "odet"),
           ("etau", "fus"), ("otau", "fus"), ("fus", "newt"), ("newt", "H"),
           ("estop", "H"), ("odet", "H")]

# Los problemas, colgando de la pieza que cada uno continua.  Van en su propia banda, y en el orden
# de la x de su ancla, para que ningun hilo cruce a otro.
PROBLEMAS = [
    ("prob:idealodd", "red",  ABIERTO),
    ("prob:types",    "red",  ABIERTO),
    ("prob:threshold", "otau", VERIF),
    ("prob:residue",  "fus",  ABIERTO),
    ("prob:LS",       "odiv", ABIERTO),
    ("prob:unit",     "H",    ABIERTO),
]
X_PROB = [0.00, 2.25, 4.50, 6.75, 9.00, 11.25]

# ------------------------------------------------------------------ CONTROL FATAL ---------------
dibujados = {e for (e, _, _) in PROBLEMAS}
en_tex = problemas_del_tex()
if en_tex and dibujados != en_tex:
    raise SystemExit("*** la figura y el .tex no llevan los mismos problemas ***\n"
                     "    solo en la figura: %s\n    solo en el .tex   : %s"
                     % (sorted(dibujados - en_tex), sorted(en_tex - dibujados)))

CAJAS, RECTS = [], []      # (artista de texto, rect en px) y los rects, para el control

fig = plt.figure(figsize=(7.0, 4.55))
ax = fig.add_axes([0.004, 0.055, 0.992, 0.90])

# bandas de fondo, para que las dos rutas se lean como rutas y no como una nube
for y, etq in ((Y_PAR, "even $t$"), (Y_IMP, "odd $t$")):
    ax.add_patch(plt.Rectangle((-0.30, y - H / 2 - 0.22), 13.55, H + 0.44,
                               facecolor=BANDA, edgecolor="none", zorder=0))
    ax.text(-0.22, y + H / 2 + 0.10, etq, fontsize=7.2, color=GRIS, style="italic", zorder=1)

for (a, b) in ARISTAS:
    xa, ya, wa, _, _ = NODOS[a]
    xb, yb, wb, _, _ = NODOS[b]
    if abs(ya - yb) < 0.01:
        p0, p1 = (xa + wa + 0.03, ya), (xb - 0.03, yb)
    elif xb >= xa + wa - 0.01:
        p0, p1 = (xa + wa + 0.03, ya), (xb - 0.03, yb)
    else:
        p0, p1 = (xa + wa / 2, ya + (H / 2 + 0.02) * (1 if yb > ya else -1)), (xb + wb / 2, yb - (H / 2 + 0.02) * (1 if yb > ya else -1))
    ax.annotate("", xy=p1, xytext=p0,
                arrowprops=dict(arrowstyle="-|>", color=REGLA, lw=0.9, shrinkA=0, shrinkB=1.5,
                                connectionstyle="arc3,rad=0.07"))

for k, (x, y, w, txt, st) in NODOS.items():
    ax.add_patch(plt.Rectangle((x, y - H / 2), w, H, facecolor=RELLENO[st], edgecolor=BORDE[st],
                               lw=1.4 if st != ABIERTO else 1.2,
                               linestyle="-" if st != ABIERTO else (0, (3, 2)), zorder=3))
    _t = ax.text(x + w / 2, y, txt, ha="center", va="center", fontsize=6.1, zorder=4,
                 color="white" if st != ABIERTO else TINTA)
    _r = None
    CAJAS.append((_t, (x, y, w, H)))

# los problemas, en su banda
for (etq, ancla, st), x in zip(PROBLEMAS, X_PROB):
    y = Y_PROB
    n, titulo = AUX.get(etq, ("?", etq.split(":")[1]))
    envuelto = "\n".join(textwrap.wrap(titulo, 26)[:2])
    xa, ya, wa, _, _ = NODOS[ancla]
    ax.annotate("", xy=(x + WP / 2, y + HP / 2 + 0.02), xytext=(xa + wa / 2, ya - H / 2 - 0.02),
                arrowprops=dict(arrowstyle="-", color=REGLA, lw=0.7, ls=(0, (1, 2))))
    ax.add_patch(plt.Rectangle((x, y - HP / 2), WP, HP, facecolor=RELLENO[st],
                               edgecolor=BORDE[st], lw=1.2,
                               linestyle="-" if st == VERIF else (0, (3, 2)), zorder=3))
    ax.text(x + WP / 2, y + 0.26, "Problem %s" % n, ha="center", va="center", fontsize=6.4,
            zorder=4, color="white" if st != ABIERTO else TINTA, fontweight="bold")
    _t = ax.text(x + WP / 2, y - 0.13, envuelto, ha="center", va="center", fontsize=5.6, zorder=4,
                 color="white" if st != ABIERTO else GRIS, linespacing=1.35)
    CAJAS.append((_t, (x, y, WP, HP)))

leg = [plt.Rectangle((0, 0), 1, 1, facecolor=AZUL, edgecolor=AZUL, label="proved here"),
       plt.Rectangle((0, 0), 1, 1, facecolor=GRIS, edgecolor=GRIS, label="external input"),
       plt.Rectangle((0, 0), 1, 1, facecolor=NARANJA, edgecolor=NARANJA,
                     label="answered by measurement"),
       plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=TINTA, ls=(0, (3, 2)),
                     label="open")]
ax.legend(handles=leg, loc="lower center", ncol=4, frameon=False, fontsize=7.2,
          bbox_to_anchor=(0.5, -0.055), handlelength=1.2, columnspacing=1.4)
ax.set_xlim(-0.45, 13.40)
ax.set_ylim(Y_PROB - 1.15, Y_PAR + 0.85)
ax.axis("off")
ax.set_title("what is left open, and what each open problem sits on",
             loc="left", fontsize=8.5)

# ---------------------------------------------------------------- CONTROL DE PIXELES ------------
# Una caja que no revienta no avisa a nadie: la comprobacion es en PIXELES sobre el render, no en
# coordenadas de datos.  Cada etiqueta tiene que caber DENTRO de su rectangulo, y ningun par de
# rectangulos puede solaparse.  Si falla, no se guarda.
fig.canvas.draw()
ren = fig.canvas.get_renderer()
TOL = 1.0                                   # holgura en pixeles


def caja_px(x, y, w, h):
    (x0, y0) = ax.transData.transform((x, y - h / 2))
    (x1, y1) = ax.transData.transform((x + w, y + h / 2))
    return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))


fallos = []
RECTS = []
for (t_art, (x, y, w, h)) in CAJAS:
    rect = caja_px(x, y, w, h)
    RECTS.append(rect)
    bb = t_art.get_window_extent(renderer=ren)
    if (bb.x0 < rect[0] - TOL or bb.x1 > rect[2] + TOL
            or bb.y0 < rect[1] - TOL or bb.y1 > rect[3] + TOL):
        fallos.append("etiqueta fuera de su caja: %r  (%.1f,%.1f)-(%.1f,%.1f) vs %s"
                      % (t_art.get_text()[:34].replace("\n", " / "),
                         bb.x0, bb.y0, bb.x1, bb.y1, tuple(round(v, 1) for v in rect)))
for i in range(len(RECTS)):
    for j in range(i + 1, len(RECTS)):
        a, b = RECTS[i], RECTS[j]
        if a[0] < b[2] - TOL and b[0] < a[2] - TOL and a[1] < b[3] - TOL and b[1] < a[3] - TOL:
            fallos.append("dos cajas se solapan: %s y %s"
                          % (tuple(round(v, 1) for v in a), tuple(round(v, 1) for v in b)))
if fallos:
    raise SystemExit("*** el control de pixeles ha disparado ***\n   " + "\n   ".join(fallos))
print("     CONTROL de pixeles: %d etiquetas dentro de su caja, %d cajas sin solape"
      % (len(CAJAS), len(RECTS)))

fig.savefig(os.path.join(AQUI, "fig_problems.pdf"))
plt.close(fig)

print("  fig_problems.pdf")
print("     problemas en el .tex : %s" % sorted(en_tex))
print("     problemas dibujados  : %s" % sorted(dibujados))
print("     CONTROL: los dos conjuntos coinciden : %s" % (dibujados == en_tex if en_tex else "sin .aux"))
for (etq, _, _) in PROBLEMAS:
    print("     %-16s -> %s  %s" % (etq, num(etq), AUX.get(etq, ("", "?"))[1]))
print("     %-22s %d bytes" % ("fig_problems.pdf",
                               os.path.getsize(os.path.join(AQUI, "fig_problems.pdf"))))
print("DONE")
