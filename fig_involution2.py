# -*- coding: utf-8 -*-
# ============================================================================================
#  fig_pairing.pdf -- POR QUE SE ANULA: la involucion que empareja TODOS los terminos.
#
#  POR QUE.  El Teorema 8.34 es el enunciado mas nuevo del paper y llega en ocho enunciados de
#  prosa sin una sola imagen.  Y es un mecanismo visual: S es simetrico respecto de C/2, una
#  transversal elige un valor por clase, la reflexion manda esa eleccion a otra, y los dos terminos
#  de Laplace se cancelan.  Lo unico que hay que ver es POR QUE no hay transversal fija.
#
#  QUE DIBUJA, y nada esta puesto a mano -- todo se calcula del beta:
#    (a) S sobre una recta, coloreado por clase, con el eje C/2 y los arcos v <-> C-v.
#        Se marcan las DOS clases que la reflexion fija y se ve que una tiene tamano PAR, luego
#        no contiene C/2: eso es exactamente la segunda clausula de (40).
#    (b) una transversal g y su reflejada g-daga, con los valores elegidos resaltados, y los dos
#        signos w(g), w(g-daga) CALCULADOS, no escritos: salen opuestos.
#    (c) el recuento entero: todas las transversales, emparejadas, con 0 puntos fijos.
#
#  CONTROL, dentro de la propia figura: se dibuja tambien un beta que cumple la simetria pero NO la
#  segunda clausula.  Alli las dos clases fijas son de tamano impar, ambas contienen C/2, y existe
#  una transversal FIJA -- que es la que impide la cancelacion.  Sin ese panel la figura ensenaria
#  que la simetria basta, y no basta.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python fig_involution2.py
# ============================================================================================

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(os.path.dirname(OUT), "gates"))
from criterion_base import setup, all_transversals            # noqa: E402

AZUL, NARANJA, TEAL, MAGENTA = "#2A78D6", "#EB6834", "#00A19A", "#B14BC8"
BANDA, REGLA, GRIS, TINTA = "#F2F1EC", "#B9B7AE", "#6B6560", "#2B2B2B"
CLASE = [AZUL, NARANJA, TEAL, MAGENTA, "#8C6BB1", "#D4A017"]

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "text.color": TINTA, "figure.facecolor": "white", "savefig.facecolor": "white",
})


def anat(beta, t, r):
    cl, E, Cd = setup(beta, t)
    S = sorted(x for k in E for x in Cd[k])
    C = S[0] + S[-1]
    inc = []
    for k in E:
        c = sorted(Cd[k], reverse=True)
        inc += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    fijas = [k for k in E if (C - k) % t == k]
    return cl, E, Cd, S, C, inc, fijas


def parejas(beta, t, r, C, E):
    """empareja las transversales por la reflexion y CUENTA los puntos fijos."""
    cl, _, _ = setup(beta, t)
    tr = all_transversals(beta, cl, r, t)
    idx = {tuple(sorted(x[0].items())): x for x in tr}
    fijas = 0
    opuestos = 0
    total = 0
    ejemplo = None
    for (sel, T, w, deg) in tr:
        ssel = {((C - k) % t if k in E else k): ((C - v) if k in E else v) for k, v in sel.items()}
        key = tuple(sorted(ssel.items()))
        if key not in idx:
            continue
        y = idx[key]
        total += 1
        if key == tuple(sorted(sel.items())):
            fijas += 1
        else:
            opuestos += (y[2] == -w)
            if ejemplo is None:
                ejemplo = (sorted(sel.values(), reverse=True), w,
                           sorted(ssel.values(), reverse=True), y[2])
    return total, fijas, opuestos, ejemplo


def busca(t, r, W, quiere_incremento):
    """El caso de control NO se pone a mano: se busca el beta mas pequeno con S simetrico y con o
    sin la segunda clausula, segun se pida.  Poner a mano un beta de control es como elegir el
    resultado; ademas el primero que escribi tenia una clase vacia y ni siquiera era del dominio."""
    import itertools
    N = t + 2 * r
    for W0 in range(N, W + 1):
        for mid in itertools.combinations(range(1, W0 + 1), N - 1):
            beta = tuple(sorted(mid, reverse=True)) + (0,)
            if beta[0] != W0:
                continue
            st = setup(beta, t)
            if st is None:
                continue
            cl, E, Cd = st
            if not E:
                continue
            S = sorted(x for k in E for x in Cd[k])
            C = S[0] + S[-1]
            if sorted(C - x for x in S) != S:
                continue
            inc = []
            for k in E:
                c = sorted(Cd[k], reverse=True)
                inc += [c[i] + c[i + 1] for i in range(len(c) - 1)]
            if (C in inc) == quiere_incremento:
                return beta
    return None


#   panel (a): cumple las dos clausulas    |    panel (b): simetrico pero SIN incremento = C
CASOS = [((12, 11, 10, 9, 3, 2, 1, 0), 4, 2, "both clauses hold"),
         (busca(4, 2, 14, False), 4, 2, "symmetric, but no increment $=C$")]
assert CASOS[1][0] is not None, "no se encontro caso de control"

fig = plt.figure(figsize=(7.0, 3.15))
for panel, (beta, t, r, titulo) in enumerate(CASOS):
    ax = fig.add_axes([0.045 + panel * 0.505, 0.10, 0.44, 0.80])
    cl, E, Cd, S, C, inc, fijas_cl = anat(beta, t, r)
    tot, fijas, opu, ej = parejas(beta, t, r, C, E)
    lo, hi = min(S), max(S)
    ax.plot([lo - 0.6, hi + 0.6], [0, 0], color=REGLA, lw=0.9, zorder=1)
    # arcos v <-> C-v
    for v in S:
        u = C - v
        if u in S and u > v:
            m = (v + u) / 2.0
            ax.plot([v, m, u], [0.16, 0.52, 0.16], color=GRIS, lw=0.8, alpha=0.75, zorder=2)
    for v in S:
        c = CLASE[(v % t) % len(CLASE)]
        ax.scatter([v], [0], s=170, color=c, edgecolor="white", lw=1.0, zorder=3)
        ax.text(v, 0, str(v), ha="center", va="center", fontsize=5.9, color="white", zorder=4)
    ax.axvline(C / 2.0, color=TINTA, lw=1.0, ls=(0, (3, 2)), zorder=1)
    ax.text(C / 2.0, -0.92, r"$C/2=%g$" % (C / 2.0), ha="center", fontsize=7, color=TINTA)
    # las clases que la reflexion fija, y su tamano
    et = []
    for k in fijas_cl:
        et.append(r"class $%d$: %d elements%s" % (k, len(Cd[k]),
                                                  " (even)" if len(Cd[k]) % 2 == 0 else " (odd)"))
    ax.text(lo - 0.6, 1.16, "\n".join(et), fontsize=6.5, color=TINTA, va="top")
    ax.text(hi + 0.6, 1.16,
            ("%d transversals, %d fixed" % (tot, fijas)) +
            ("\nall pair up with opposite sign" if fijas == 0
             else "\na fixed one cannot cancel"),
            fontsize=6.5, color=(TINTA if fijas == 0 else NARANJA), va="top", ha="right")
    if ej is not None and fijas == 0:
        g, w, gd, wd = ej
        ax.text(lo - 0.6, -1.35,
                r"$g=\{%s\}$, $w=%+d$" % (",".join(map(str, g)), w) + "\n" +
                r"$g^{\dagger}=\{%s\}$, $w=%+d$" % (",".join(map(str, gd)), wd),
                fontsize=6.4, color=TINTA, va="top")
    ax.set_xlim(lo - 1.4, hi + 1.4)
    ax.set_ylim(-2.15, 1.35)
    ax.axis("off")
    ax.set_title("(%s) %s" % ("ab"[panel], titulo), loc="left", fontsize=8.2)

fig.savefig(os.path.join(OUT, "fig_pairing.pdf"))
plt.close(fig)

print("  fig_pairing.pdf -- la involucion, con su control")
for beta, t, r, titulo in CASOS:
    cl, E, Cd, S, C, inc, fijas_cl = anat(beta, t, r)
    tot, fijas, opu, _ = parejas(beta, t, r, C, E)
    print("     t=%d r=%d beta=%s" % (t, r, list(beta)))
    print("        C=%d  clases fijas %s  tamanos %s  incremento=C: %s"
          % (C, fijas_cl, [len(Cd[k]) for k in fijas_cl], C in inc))
    print("        %d transversales, %d fijas, %d parejas con signo opuesto" % (tot, fijas, opu))
print("  DONE")
