# -*- coding: utf-8 -*-
# LA FIBRA Y EL DETERMINANTE SON LA MISMA SUMA.   16 de agosto de 2026.
#
# POR QUE ESTA FIGURA.  El paper enuncia dos cosas que parecen distintas y son la misma:
#
#   Prop. (inversion cerrada)   c(X) = sum_{k impar} nu~(X + t k)      -- una suma sobre una FIBRA
#   Prop. (determinante)        c(X) = eps_t . det M(Lambda, X)        -- una suma sobre PERMUTACIONES
#
# La segunda es la primera reagrupada: cada termino no nulo del desarrollo del determinante ES un
# corte de la progresion con el soporte de nu.  Dibujarlas lado a lado es la unica manera de que eso
# se vea, y de paso ensena la dicotomia: una fibra que corta UNA vez da +-1, y una que corta DOS se
# cancela exactamente -- que es (L1) entera.
#
# Se dibuja el caso mas pequeno donde las dos cosas ocurren: t=3, r=2, Lambda=(2,1,0), o sea
# V = (9,5,1) y matrices 3x3.  Arriba X=(6,2), que corta una vez.  Abajo X=(6,-4), que corta dos.
#
# CONTROLES -- la figura SE NIEGA A DIBUJARSE si alguno falla
#   C1  para los dos X dibujados:  c por la formula de progresion  ==  eps_t . det M.
#   C2  los terminos de permutacion que se marcan son EXACTAMENTE los no nulos del desarrollo, y su
#       suma con signo es det M.  (Si marcara los que me convienen, la figura mentiria.)
#   C3  el numero de cortes dibujados es 1 arriba y 2 abajo, contados sobre el soporte de verdad.
#   C4  el soporte dibujado es el de nu_de, no uno escrito a mano.
#   C5  ninguna etiqueta se solapa con otra, medido en PIXELES sobre el lienzo renderizado.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python fig_determinant.py

import itertools
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gates"))
from divided_differences import nu_de, nu_extendida, enderezar_D, eps_t          # noqa: E402
from unimodularidad_barrido import matriz, det_entero                            # noqa: E402

INK, MUTED, GRID = "#1b1b1b", "#8a8a8a", "#dedede"
PLUS, MINUS, SECOND = "#1f6f8b", "#c1462f", "#5a5a5a"
SURFACE = "#fcfcfb"

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

T, R, MP = 3, 2, 1
RP = MP + R
LAM = [2, 1, 0]
V = [2 * LAM[i] + 2 * (RP - i) - 1 for i in range(RP)]
E = eps_t(T, MP)
NU = nu_de(LAM, T, R)
CASOS = [((6, 2), "the fibre meets the support once"),
         ((6, -4), "it meets it twice, and they cancel")]


def fmt(v):
    """Un cero es 0, no +0."""
    return "0" if v == 0 else "%+d" % v


def terminos_permutacion(M):
    """Los terminos NO NULOS del desarrollo de det M, como (permutacion, signo, valor)."""
    out = []
    n = len(M)
    for perm in itertools.permutations(range(n)):
        val = 1
        for i in range(n):
            val *= M[i][perm[i]]
            if val == 0:
                break
        if val == 0:
            continue
        s, visto = 1, [False] * n
        for i in range(n):
            if visto[i]:
                continue
            j, L = i, 0
            while not visto[j]:
                visto[j] = True
                j = perm[j]
                L += 1
            if L % 2 == 0:
                s = -s
        out.append((perm, s, s * val))
    return out


def cortes(X, tope):
    """Los (k, Y, valor) de la progresion que caen en el soporte extendido."""
    out = []
    for k in itertools.product(range(1, tope + 1, 2), repeat=R):
        Y = tuple(X[j] + T * k[j] for j in range(R))
        v = nu_extendida(NU, Y)
        if v:
            out.append((k, Y, v))
    return out


# ------------------------------------------------------------------ controles
TOPE = 9
fallos = []
datos = []
for (X, titulo) in CASOS:
    M = matriz(V, list(X), T, MP, R)
    d = det_entero(M)
    cs = cortes(X, TOPE)
    c_fibra = sum(v for _, _, v in cs)
    tp = terminos_permutacion(M)
    if c_fibra != E * d:
        fallos.append("C1 X=%s: fibra %d != eps.det %d" % (str(X), c_fibra, E * d))
    if sum(v for _, _, v in tp) != d:
        fallos.append("C2 X=%s: los terminos marcados no suman det" % str(X))
    if len(tp) != len(cs):
        fallos.append("C2 X=%s: %d terminos de permutacion contra %d cortes"
                      % (str(X), len(tp), len(cs)))
    datos.append({"X": X, "titulo": titulo, "M": M, "det": d, "cortes": cs,
                  "terminos": tp, "c": c_fibra})

if len(datos[0]["cortes"]) != 1 or len(datos[1]["cortes"]) != 2:
    fallos.append("C3 el numero de cortes no es 1 y 2: %d y %d"
                  % (len(datos[0]["cortes"]), len(datos[1]["cortes"])))
if not NU:
    fallos.append("C4 soporte vacio")
if fallos:
    print("LA FIGURA NO SE DIBUJA.  Controles fallidos:")
    for f in fallos:
        print("   " + f)
    sys.exit(1)

print("  controles C1 C2 C3 C4 ok")
print("  V = %s   eps_t = %+d   soporte de nu: %s" % (V, E, {k: v for k, v in sorted(NU.items())}))
for D in datos:
    print("  X=%s : %d corte(s), c = %+d, det M = %+d, %d termino(s) de permutacion"
          % (str(D["X"]), len(D["cortes"]), D["c"], D["det"], len(D["terminos"])))

# ------------------------------------------------------------------ dibujo
fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.4),
                         gridspec_kw={"width_ratios": [1.32, 1.0]})

VX0, VX1, VY0, VY1 = -12, 18, -12, 18   # ventana, asimetrica: la progresion va hacia arriba-derecha
PUNTOS = []                              # posiciones dibujadas, para el control de solape en pixeles
PENDIENTES = []                          # etiquetas por colocar, con su ancla


def panel_reticulo(ax, D):
    X = D["X"]
    puntos = []
    # el soporte extendido de nu dentro de la ventana
    for a in range(VX0, VX1 + 1):
        for b in range(VY0, VY1 + 1):
            v = nu_extendida(NU, (a, b))
            if v:
                ax.scatter([a], [b], s=52, marker="o",
                           facecolors=(PLUS if v > 0 else MINUS), edgecolors="none",
                           zorder=3, alpha=0.85)
                puntos.append((a, b))
    # la progresion X + t k, k impar: retículo de paso 2t que arranca en X + t(1,...,1)
    prog, hit = [], []
    for k in itertools.product(range(1, TOPE + 1, 2), repeat=R):
        Y = (X[0] + T * k[0], X[1] + T * k[1])
        if not (VX0 <= Y[0] <= VX1 and VY0 <= Y[1] <= VY1):
            continue
        prog.append(Y)
        if nu_extendida(NU, Y):
            hit.append(Y)
    xs = sorted({y[0] for y in prog})
    ys = sorted({y[1] for y in prog})
    for a in xs:
        ax.plot([a, a], [min(ys), max(ys)], color=GRID, linewidth=0.7, zorder=1)
    for b in ys:
        ax.plot([min(xs), max(xs)], [b, b], color=GRID, linewidth=0.7, zorder=1)
    ax.plot([X[0], xs[0]], [X[1], ys[0]], color=MUTED, linewidth=0.9, linestyle=(0, (2, 2)),
            zorder=1)
    ax.scatter([p[0] for p in prog], [p[1] for p in prog], s=30, marker="s",
               facecolors=SURFACE, edgecolors=MUTED, linewidths=0.9, zorder=2)
    puntos.extend(prog)
    for Y in hit:
        ax.scatter([Y[0]], [Y[1]], s=210, facecolors="none",
                   edgecolors=INK, linewidths=1.2, zorder=4)
    ax.scatter([X[0]], [X[1]], s=64, marker="D", facecolors=SURFACE,
               edgecolors=INK, linewidths=1.2, zorder=4)
    puntos.append(X)
    # Las etiquetas NO se colocan a ojo: se prueban candidatos y se queda el primero que no pisa
    # nada.  Adivinar el desplazamiento fue lo que rompio la primera version de esta figura.
    PENDIENTES.append((ax, [(r"$\tilde\nu=%s$" % fmt(nu_extendida(NU, Y)), Y, 8.0) for Y in hit]
                       + [(r"$X=(%d,%d)$" % X, X, 8.4)]))
    ax.axhline(0, color="#eeece7", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#eeece7", linewidth=0.8, zorder=0)
    ax.set_xlim(VX0 - 1, VX1 + 1)
    ax.set_ylim(VY0 - 1, VY1 + 1)
    ax.set_aspect("equal")
    ax.set_xticks([-9, 0, 9])
    ax.set_yticks([-9, 0, 9])
    ax.tick_params(labelsize=7.5, colors=SECOND)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.set_title(D["titulo"], fontsize=9.5, color=INK, pad=7)
    ax.text(0.5, -0.165, r"$c=\sum_k\tilde\nu(X+tk)=%s$" % fmt(D["c"]),
            transform=ax.transAxes, ha="center", fontsize=9, color=INK)
    PUNTOS.append((ax, puntos))


def panel_matriz(ax, D):
    M = D["M"]
    n = len(M)
    for i in range(n):
        for j in range(n):
            ax.text(j, -i, "%+d" % M[i][j] if M[i][j] else "0",
                    ha="center", va="center", fontsize=11,
                    color=(INK if M[i][j] else "#c4c2bc"))
    # Los terminos de permutacion, marcados.  Dos terminos pueden COMPARTIR una celda -- abajo los
    # dos usan V_1 -> X_1 -- y si se dibujan encima solo se ve el ultimo.  Se anidan.
    for m, (perm, s, val) in enumerate(D["terminos"]):
        col = PLUS if val > 0 else MINUS
        d_ = 0.055 * m
        for i in range(n):
            ax.add_patch(Rectangle((perm[i] - 0.34 + d_, -i - 0.30 + d_),
                                   0.68 - 2 * d_, 0.60 - 2 * d_,
                                   fill=False, edgecolor=col, linewidth=1.35, zorder=3))
    ax.text(-1.05, 0.0, r"$V_1{=}9$", ha="right", va="center", fontsize=7.6, color=SECOND)
    ax.text(-1.05, -1.0, r"$V_2{=}5$", ha="right", va="center", fontsize=7.6, color=SECOND)
    ax.text(-1.05, -2.0, r"$V_3{=}1$", ha="right", va="center", fontsize=7.6, color=SECOND)
    ax.text(0, 0.72, "frozen", ha="center", fontsize=7.6, color=SECOND)
    ax.text(1, 0.72, r"$X_1$", ha="center", fontsize=7.6, color=SECOND)
    ax.text(2, 0.72, r"$X_2$", ha="center", fontsize=7.6, color=SECOND)
    ax.set_xlim(-2.1, 2.7)
    ax.set_ylim(-2.9, 1.25)
    ax.axis("off")
    piezas = "+".join(r"(%s)" % fmt(val) for (_, _, val) in D["terminos"])
    ax.text(0.5, -0.03, r"$\det M=%s=%s$" % (piezas, fmt(D["det"])),
            transform=ax.transAxes, ha="center", fontsize=9, color=INK)
    ax.text(0.5, -0.165, r"$\epsilon_t\det M=%s$" % fmt(E * D["det"]),
            transform=ax.transAxes, ha="center", fontsize=9, color=INK)


for fila, D in enumerate(datos):
    panel_reticulo(axes[fila][0], D)
    panel_matriz(axes[fila][1], D)

fig.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.10, hspace=0.42, wspace=0.02)

# ------------------------------------------------------------------ colocacion guiada por el control
CANDIDATOS = [(2.2, 1.8, "left", "bottom"), (-2.2, 1.8, "right", "bottom"),
              (2.2, -1.8, "left", "top"), (-2.2, -1.8, "right", "top"),
              (0.0, 3.0, "center", "bottom"), (0.0, -3.0, "center", "top"),
              (4.5, 0.0, "left", "center"), (-4.5, 0.0, "right", "center"),
              (5.5, 3.2, "left", "bottom"), (-5.5, 3.2, "right", "bottom")]
RADIO = 5.0     # pixeles, el radio aproximado de un marcador


def choca(ax, caja, puntos, otras):
    for (px, py) in puntos:
        dx, dy = ax.transData.transform((px, py))
        if caja.x0 - RADIO < dx < caja.x1 + RADIO and caja.y0 - RADIO < dy < caja.y1 + RADIO:
            return True
    for b in otras:
        if caja.x0 < b.x1 and b.x0 < caja.x1 and caja.y0 < b.y1 and b.y0 < caja.y1:
            return True
    return False


fig.canvas.draw()
ren = fig.canvas.get_renderer()
for (ax, etiquetas), (_, puntos) in zip(PENDIENTES, PUNTOS):
    puestas = [t_.get_window_extent(renderer=ren) for t_ in ax.texts]
    for (txt, ancla, fs) in etiquetas:
        colocada = None
        for (dx, dy, ha, va) in CANDIDATOS:
            t_ = ax.text(ancla[0] + dx, ancla[1] + dy, txt, fontsize=fs, color=INK,
                         ha=ha, va=va, zorder=6)
            caja = t_.get_window_extent(renderer=ren)
            if not choca(ax, caja, puntos, puestas):
                colocada = caja
                break
            t_.remove()
        if colocada is None:
            print("  NO HAY SITIO para la etiqueta %s en %s -- la figura no se dibuja" % (txt, ancla))
            sys.exit(1)
        puestas.append(colocada)

# ------------------------------------------------------------------ C5  solape en PIXELES
fig.canvas.draw()
cajas = []
for ax in axes.ravel():
    for t_ in ax.texts:
        b = t_.get_window_extent(renderer=fig.canvas.get_renderer())
        cajas.append((t_.get_text(), b, ax))
solapes = []
for i in range(len(cajas)):
    for j in range(i + 1, len(cajas)):
        if cajas[i][2] is not cajas[j][2]:
            continue
        a, b = cajas[i][1], cajas[j][1]
        if a.x0 < b.x1 and b.x0 < a.x1 and a.y0 < b.y1 and b.y0 < a.y1:
            solapes.append((cajas[i][0], cajas[j][0]))
# y texto CONTRA PUNTO DIBUJADO: una etiqueta encima de un dato es igual de ilegible que dos
# etiquetas encima, y la version anterior de esta figura tenia justo eso sin que el control lo viera.
RADIO = 5.0     # pixeles, el radio aproximado de un marcador
for (ax, pts) in PUNTOS:
    for (px, py) in pts:
        dx, dy = ax.transData.transform((px, py))
        for (txt, b, axt) in cajas:
            if axt is not ax:
                continue
            if b.x0 - RADIO < dx < b.x1 + RADIO and b.y0 - RADIO < dy < b.y1 + RADIO:
                solapes.append(("punto (%d,%d)" % (px, py), txt))
if solapes:
    print("  C5 FALLA: solapes -> %s" % solapes[:6])
    sys.exit(1)
print("  C5 sin solapes de etiqueta ni etiqueta-sobre-punto (medido en pixeles)")

fig.savefig("fig_determinant.pdf")
print("  fig_determinant.pdf escrito")
