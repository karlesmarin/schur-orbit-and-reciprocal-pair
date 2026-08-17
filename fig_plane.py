# -*- coding: utf-8 -*-
# ============================================================================================
#  fig_plane.pdf -- EL PLANO (t, r), que la introduccion promete dibujar y no dibujaba.
#
#  POR QUE.  El parrafo "It is worth drawing the map" enumera en prosa tres lineas resueltas y una
#  region abierta, y el lector tiene que sostener cuatro enunciados a la vez para verlo.  Es una
#  imagen, y cabe en un cuarto de pagina.
#
#  QUE DIBUJA, y NADA se coloca a mano: para cada (t, r) del rango se decide el estatus con las
#  mismas reglas que el paper enuncia, escritas aqui una sola vez:
#     r = 1        -> Teorema 3.1: la forma cerrada, o sea el valor y no solo su anulacion
#     t = 2        -> Teorema 8.6: criterio completo, con el input externo (PvW)
#     t impar      -> Corolario 8.22: criterio completo, SIN input externo
#     el resto     -> Corolario 8.34 da una implicacion; la conversa es Conjetura 8.43
#  El (2,1) se marca aparte porque es la esquina donde las dos mitades del paper se tocan.
#
#  CONTROL.  El recuento de celdas de cada estatus se imprime, y la suma tiene que ser el total de
#  la rejilla.  Una figura que colorea celdas sin contarlas puede perder o duplicar una en silencio.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python fig_plane.py
# ============================================================================================

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def num_del_aux(etiqueta, por_defecto):
    """El numero de un enunciado, leido del .aux y no escrito a mano.

    Un numero copiado dentro de una figura no se entera de nada. Al insertar el teorema de
    suficiencia y la conjetura, LaTeX renumero lo que iba detras y esta leyenda se quedo diciendo
    Cor 8.20 y Conj 8.39 cuando ya eran 8.21 y 8.42 -- junto a un pie que daba los numeros buenos,
    porque el pie sale del .tex. Ninguna guardia lo vio: en el .tex esos numeros no existen.
    El por_defecto es solo para poder dibujar sin haber compilado; si sale ese, es que falta el aux.
    """
    for nombre in ("_v2_trial.aux", "orbit_pair.aux"):
        p = os.path.join(OUT, nombre)
        if os.path.exists(p):
            for l in open(p, encoding="utf-8", errors="replace"):
                if l.startswith(r"\newlabel{%s}" % etiqueta):
                    try:
                        return l.split("{{")[1].split("}")[0]
                    except IndexError:
                        pass
    return por_defecto


AZUL, NARANJA, TEAL = "#2A78D6", "#EB6834", "#00A19A"
BANDA, REGLA, GRIS, TINTA = "#F2F1EC", "#B9B7AE", "#6B6560", "#2B2B2B"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "text.color": TINTA, "figure.facecolor": "white", "savefig.facecolor": "white",
})

TS = list(range(2, 11))
RS = list(range(1, 6))


def estatus(t, r):
    """las reglas del paper, escritas una sola vez.  El orden importa: r=1 da el VALOR."""
    if r == 1:
        return "valor"
    if t % 2 == 1:
        return "criterio_libre"
    if t == 2:
        return "criterio_pvw"
    return "media"


COL = {"valor": AZUL, "criterio_libre": TEAL, "criterio_pvw": NARANJA, "media": "white"}
NOMBRE = {
    "valor": "the value itself (Thm %s)" % num_del_aux("thm:main", "3.1"),
    "criterio_libre": "criterion, no external input (Cor %s)" % num_del_aux("cor:oddgen", "8.22"),
    "criterio_pvw": "criterion, modulo one input (Thm %s)" % num_del_aux("conj:crit", "8.6"),
    "media": "one implication proved (Thm %s);\nconverse conjectural (Conj %s)"
             % (num_del_aux("thm:suff", "8.35"), num_del_aux("conj:general", "8.43")),
}

fig = plt.figure(figsize=(7.0, 2.45))
ax = fig.add_axes([0.055, 0.17, 0.55, 0.74])

cuenta = {k: 0 for k in COL}
for t in TS:
    for r in RS:
        st = estatus(t, r)
        cuenta[st] += 1
        relleno = COL[st]
        ax.add_patch(plt.Rectangle((t - 0.44, r - 0.40), 0.88, 0.80,
                                   facecolor=relleno, edgecolor=REGLA if st == "media" else "white",
                                   lw=0.9, linestyle="-" if st != "media" else (0, (2.2, 1.6)),
                                   zorder=2))
# la esquina donde las dos mitades se tocan
ax.scatter([2], [1], s=54, marker="o", facecolor="white", edgecolor=TINTA, lw=1.3, zorder=4)
ax.annotate(r"$(t,r)=(2,1)$: $\Phi_2=\Psi_1$", xy=(2, 1), xytext=(3.15, 1.95),
            fontsize=6.8, color=TINTA,
            arrowprops=dict(arrowstyle="-", color=TINTA, lw=0.7))

ax.set_xticks(TS)
ax.set_yticks(RS)
ax.set_xlabel(r"$t$  (size of the frozen orbit $\mu_t$)", fontsize=8)
ax.set_ylabel(r"$r$  (free pairs)", fontsize=8)
ax.set_xlim(1.4, max(TS) + 0.6)
ax.set_ylim(0.4, max(RS) + 0.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("bottom", "left"):
    ax.spines[s].set_color(REGLA)
ax.tick_params(length=2.5, labelsize=7.5, color=REGLA)
ax.set_title("where the zero locus is known, and how", loc="left", fontsize=8.5)

# leyenda a la derecha, en texto: los nombres son largos y una leyenda de cajas los aplasta
lx = 0.635
ly = 0.80
for k in ("valor", "criterio_libre", "criterio_pvw", "media"):
    fig.patches.append(plt.Rectangle((lx, ly - 0.035), 0.022, 0.055,
                                     facecolor=COL[k], edgecolor=REGLA if k == "media" else "none",
                                     linestyle="-" if k != "media" else (0, (2.2, 1.6)),
                                     transform=fig.transFigure, zorder=3))
    fig.text(lx + 0.032, ly - 0.008, NOMBRE[k], fontsize=6.6, va="center", color=TINTA)
    ly -= 0.135 if "\n" not in NOMBRE[k] else 0.175

# Tres lineas y no dos: a dos, la segunda medía 3 pt mas que la figura y matplotlib le cortaba la
# ultima letra, de modo que en el PDF ponia "is know".  Recortado no se ve, y el pie de figura si
# tiene sitio para decirlo entero.
fig.text(lx, 0.10, "columns of odd $t$ are settled outright;\n"
                   "the row $r=1$ is the only place the value,\n"
                   "and not merely its vanishing, is known",
         fontsize=6.3, color=GRIS, linespacing=1.45)

fig.savefig(os.path.join(OUT, "fig_plane.pdf"))
plt.close(fig)

total = sum(cuenta.values())
print("  fig_plane.pdf -- rejilla t=%d..%d, r=%d..%d" % (TS[0], TS[-1], RS[0], RS[-1]))
for k, v in cuenta.items():
    print("     %-16s %3d celdas" % (k, v))
print("     %-16s %3d celdas   (rejilla %d x %d = %d)"
      % ("TOTAL", total, len(TS), len(RS), len(TS) * len(RS)))
assert total == len(TS) * len(RS), "el recuento no cuadra con la rejilla"
print("  DONE")
