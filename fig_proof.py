# -*- coding: utf-8 -*-
# ============================================================================================
#  fig_laplace.pdf -- la ARQUITECTURA de la prueba del Teorema 3.1, que son 142 lineas de prosa
#  sin un solo elemento visual.  14 de agosto de 2026.
#
#  QUE DIBUJA.  Lo que dice el parrafo de apertura de la seccion: "expandir por las t filas
#  congeladas deja un termino por cada eleccion de las DOS columnas que no se congelan, y el conteo
#  de exceso dos dice que toda eleccion cae en exactamente uno de tres tipos combinatorios".
#  Izquierda: el determinante N x N partido por Laplace, con las t filas congeladas, el menor
#  M_S = +- V (Lema 4.1) y el bloque 2 x 2 libre f(beta_{j1} - beta_{j2}) (Lema 4.3).
#  Derecha: los TRES TIPOS, uno por fila, con un beta REAL cada uno y el numero de pares U que
#  sobreviven CONTADO, no escrito a mano: 0, 4 y 3.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python fig_proof.py
# ============================================================================================

import os
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AZUL, NARANJA, TEAL, MAGENTA = "#2A78D6", "#EB6834", "#00A19A", "#B14BC8"
BANDA, REGLA, GRIS, TINTA = "#F2F1EC", "#B9B7AE", "#6B6560", "#2B2B2B"
CLASE = [AZUL, NARANJA, TEAL, MAGENTA]

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "text.color": TINTA, "figure.facecolor": "white", "savefig.facecolor": "white",
})


def sobreviven(beta, t):
    """pares U = {j1<j2} cuyo complemento lleva los t residuos.  Se CUENTA, no se afirma."""
    N = len(beta)
    out = []
    for U in combinations(range(N), 2):
        S = [j for j in range(N) if j not in U]
        if len(set(beta[j] % t for j in S)) == t:
            out.append(U)
    return out


def perfil(beta, t):
    d = {}
    for v in beta:
        d[v % t] = d.get(v % t, 0) + 1
    return tuple(d.get(i, 0) for i in range(t))


EJEMPLOS = [((8, 6, 5, 4, 1, 0), "a residue class is empty"),
            ((9, 8, 7, 6, 5, 4), "two classes of size two"),
            ((12, 8, 7, 6, 5, 4), "one class of size three")]

fig = plt.figure(figsize=(7.0, 3.15))
ax = fig.add_axes([0.015, 0.10, 0.44, 0.84])
bx = fig.add_axes([0.53, 0.10, 0.46, 0.84])

# ---------------------------------------------------------------- izquierda: el Laplace ---------
t, N = 4, 6
U = (1, 4)                                   # las dos columnas NO congeladas del dibujo
filas = [r"$1$", r"$\zeta$", r"$\zeta^{2}$", r"$\zeta^{3}$", r"$z$", r"$z^{-1}$"]
for i in range(N):
    for j in range(N):
        congelada = i < t
        enU = j in U
        if congelada and not enU:
            fc, ec = BANDA, REGLA
        elif not congelada and enU:
            fc, ec = "#DCE9FA", AZUL
        else:
            fc, ec = "white", REGLA
        ax.add_patch(plt.Rectangle((j, -i), 0.92, 0.92, facecolor=fc, edgecolor=ec, lw=0.9))
        ax.text(j + 0.46, -i + 0.46, r"$x_{%d}^{\beta_{%d}}$" % (i + 1, j + 1),
                ha="center", va="center", fontsize=6.2,
                color=TINTA if (congelada != enU) else "#9A9A9A")
for i, f in enumerate(filas):
    ax.text(-0.28, -i + 0.46, f, ha="right", va="center", fontsize=8,
            color=GRIS if i < t else AZUL)
for j in range(N):
    ax.text(j + 0.46, 1.12, r"$\beta_{%d}$" % (j + 1), ha="center", fontsize=7,
            color=AZUL if j in U else GRIS)
ax.plot([-0.05, N - 0.02], [-t + 0.98, -t + 0.98], color=TINTA, lw=1.1)
ax.annotate(r"$M_S=(-1)^{\mathrm{inv}(b_S)}V$", xy=(1.4, -1.6), xytext=(1.4, -6.35),
            ha="center", fontsize=8, color=GRIS,
            arrowprops=dict(arrowstyle="-", color=GRIS, lw=0.8))
ax.text(1.4, -6.95, r"Lemma 4.1: zero unless the $t$", ha="center", fontsize=6.6, color=GRIS)
ax.text(1.4, -7.45, r"residues on $S$ are distinct", ha="center", fontsize=6.6, color=GRIS)
ax.annotate(r"$f(\beta_{j_1}-\beta_{j_2})$", xy=(4.6, -4.6), xytext=(5.1, -6.35),
            ha="center", fontsize=8, color=AZUL,
            arrowprops=dict(arrowstyle="-", color=AZUL, lw=0.8))
ax.text(5.1, -6.95, r"the free $2\times2$ block", ha="center", fontsize=6.6, color=AZUL)
ax.set_xlim(-1.5, N + 0.3)
ax.set_ylim(-7.9, 1.7)
ax.axis("off")
ax.set_title(r"one term per choice of the two free columns $U=\{j_1<j_2\}$",
             loc="left", fontsize=8.5)

# ---------------------------------------------------------------- derecha: los tres tipos -------
bx.set_title("and every choice falls into one of three types", loc="left", fontsize=8.5)
info = []
for fila, (beta, nombre) in enumerate(EJEMPLOS):
    surv = sobreviven(beta, t)
    p = perfil(beta, t)
    info.append((beta, p, len(surv)))
    y = -fila * 1.62
    for k, v in enumerate(beta):
        c = CLASE[v % t]
        bx.scatter([k * 0.62], [y], s=190, marker="o", color=c, edgecolor="white", linewidth=1.0, zorder=3)
        bx.text(k * 0.62, y, str(v), ha="center", va="center", fontsize=6.2, color="white", zorder=4)
    bx.text(-0.55, y, r"$\beta$", ha="right", va="center", fontsize=8, color=GRIS)
    bx.text(3.95, y + 0.30, nombre, fontsize=7.6, color=TINTA, va="center")
    bx.text(3.95, y - 0.26, r"profile $(%s)$" % ",".join(map(str, p)), fontsize=7, color=GRIS,
            va="center")
    bx.text(8.35, y, r"$\mathbf{%d}$" % len(surv), fontsize=13, color=TINTA, ha="center",
            va="center", fontweight="bold")
    bx.text(8.35, y - 0.44, "terms", fontsize=6.4, color=GRIS, ha="center", va="center")
leg = [plt.Line2D([], [], marker="o", ls="", color=CLASE[i], mec="white", ms=8,
                  label=r"residue $%d$" % i) for i in range(t)]
bx.legend(handles=leg, loc="lower center", ncol=4, frameon=False, fontsize=7,
          bbox_to_anchor=(0.45, -0.10), handletextpad=0.2, columnspacing=1.0)
bx.set_xlim(-1.1, 9.1)
bx.set_ylim(-4.35, 0.75)
bx.axis("off")

fig.savefig(os.path.join(OUT, "fig_laplace.pdf"))
plt.close(fig)

# ============================================================================================
#  fig_map.pdf -- que se apoya en que, y con que estatus.  La seccion de problemas abiertos son
#  392 lineas de prosa sin un solo elemento visual, y es justo donde el lector necesita ver de un
#  vistazo lo que esta probado, lo que solo esta medido, lo que es de otros y lo que sigue abierto.
#
#  LOS NOMBRES, NO LOS NUMEROS.  Los nodos llevan el enunciado por su nombre y no por su numero,
#  porque los numeros se mueven al editar y una figura que los lleva dentro miente en silencio.
#  El unico numero que aparece, el del teorema principal de la seccion, se LEE del .aux si existe.
# ============================================================================================
PROBADO, VERIF, EXTERNO, ABIERTO = "proved", "verified", "external", "open"
COL = {PROBADO: AZUL, VERIF: NARANJA, EXTERNO: GRIS, ABIERTO: "white"}
BORDE = {PROBADO: AZUL, VERIF: NARANJA, EXTERNO: GRIS, ABIERTO: TINTA}


def num_del_aux(etiqueta, por_defecto):
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


NUM = num_del_aux("conj:crit", "8.6")

# (x, y, ancho, texto, estatus)
NODOS = {
    "pvw":   (0.00, 2.05, 1.62, "external:\nrigidity of $s_\\lambda s_\\mu$", EXTERNO),
    "lap":   (0.00, 0.75, 1.62, "Laplace along the $t$\nfrozen rows", PROBADO),
    "inc":   (0.00, -0.35, 1.62, "increments and the\n2-torsion of $\\mathbb{Z}/t$", PROBADO),
    "dic":   (2.05, 2.05, 1.62, "the dictionary:\ntwo Schur factors", PROBADO),
    "G2":    (2.05, -0.35, 1.62, "$|G|\\leq2$, tie in\nclasses $i,\\ i{+}t/2$", PROBADO),
    "refl":  (4.10, 1.35, 1.62, "$T_{\\mathrm{b}}=\\tau-T_{\\mathrm{a}}$:\nthe two reflect", PROBADO),
    "a1":    (6.15, 1.95, 1.62, "the extremes of $\\mathcal{S}$\navoid $g_{\\mathrm{com}}$", PROBADO),
    "ctau":  (6.15, 0.75, 1.62, "$C=\\tau$", PROBADO),
    "crit":  (8.20, 1.35, 1.78, "the criterion, every $r$\n(Theorem %s)" % NUM, PROBADO),
    "litt":  (4.10, -1.55, 1.62, "Littlewood's reduction\nto the $C_\\mu$", PROBADO),
    "cert":  (6.15, -1.55, 1.62, "a certificate for it\n(refuted)", PROBADO),
    "p105":  (8.20, -1.55, 1.78, "OPEN: explain the\n$C_\\mu$ (Problem %s)"
              % num_del_aux("prob:instrument", "10.5"), ABIERTO),
    # ---- la tercera via, anadida el 14 de agosto: la direccion SUFICIENTE para todo t y r.
    #      Empieza en la hipotesis, no en el Laplace, y por eso arranca su propia fila: no es un
    #      refinamiento del argumento extremal sino otro mecanismo, una involucion GLOBAL.
    "hyp":   (0.00, 3.25, 1.62, "$C-\\mathcal{S}=\\mathcal{S}$ and\nsome $\\Delta_i(k)=C$", PROBADO),
    "refla": (2.05, 3.25, 1.62, "the reflection acts on\nevery transversal", PROBADO),
    "sgn":   (4.10, 3.25, 1.62, "$w(g^\\dagger)=-w(g)$\nand $A(C-T)=A(T)$", PROBADO),
    "suff":  (6.15, 3.25, 1.62, "vanishing, every $t$\nand $r$ (Thm %s)"
              % num_del_aux("thm:suff", "8.35"), PROBADO),
    "conj":  (8.20, 3.25, 1.78, "OPEN: the converse\n(Conjecture %s)"
              % num_del_aux("conj:general", "8.43"), ABIERTO),
}
ARISTAS = [("lap", "dic"), ("lap", "G2"), ("inc", "G2"), ("pvw", "dic"), ("dic", "refl"),
           ("G2", "refl"), ("refl", "a1"), ("refl", "ctau"), ("a1", "ctau"), ("ctau", "crit"),
           ("G2", "crit"), ("litt", "cert"), ("cert", "p105"), ("litt", "p105"),
           ("hyp", "refla"), ("refla", "sgn"), ("sgn", "suff"), ("suff", "conj")]

figm = plt.figure(figsize=(7.0, 4.55))
mx = figm.add_axes([0.005, 0.07, 0.99, 0.88])
H = 0.72
for (a, b) in ARISTAS:
    xa, ya, wa, _, _ = NODOS[a]
    xb, yb, wb, _, _ = NODOS[b]
    mx.annotate("", xy=(xb - 0.04, yb), xytext=(xa + wa + 0.04, ya),
                arrowprops=dict(arrowstyle="-|>", color=REGLA, lw=0.9,
                                connectionstyle="arc3,rad=0.06", shrinkA=0, shrinkB=2))
for k, (x, y, w, txt, st) in NODOS.items():
    c = COL[st]
    relleno = "white" if st == ABIERTO else c
    mx.add_patch(plt.Rectangle((x, y - H / 2), w, H, facecolor=relleno, edgecolor=BORDE[st],
                               lw=1.4 if st != ABIERTO else 1.2,
                               linestyle="-" if st != ABIERTO else (0, (3, 2)), zorder=3))
    mx.text(x + w / 2, y, txt, ha="center", va="center", fontsize=6.4, zorder=4,
            color="white" if st in (PROBADO, VERIF, EXTERNO) else TINTA)
mx.text(NODOS["cert"][0] + NODOS["cert"][2] / 2, NODOS["cert"][1] - 0.55,
        "a counterexample, so proved", ha="center", fontsize=6.2, color=GRIS)
leg = [plt.Rectangle((0, 0), 1, 1, facecolor=AZUL, edgecolor=AZUL, label="proved here"),
       plt.Rectangle((0, 0), 1, 1, facecolor=GRIS, edgecolor=GRIS, label="external input"),
       plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=TINTA, ls=(0, (3, 2)),
                     label="open")]
mx.legend(handles=leg, loc="lower center", ncol=3, frameon=False, fontsize=7.5,
          bbox_to_anchor=(0.5, -0.07), handlelength=1.3)
mx.set_xlim(-0.15, 10.15)
mx.set_ylim(-2.55, 3.95)
mx.axis("off")
mx.set_title("three routes to the zero locus, and what each one reaches",
             loc="left", fontsize=8.5)
figm.savefig(os.path.join(OUT, "fig_map.pdf"))
plt.close(figm)

print("  fig_laplace.pdf  --  los tres tipos, con el conteo HECHO y no afirmado:")
for beta, p, n in info:
    print("     beta=%-24s perfil=%-12s pares U que sobreviven: %d" % (str(beta), str(p), n))
print("")
print("     (el texto de la seccion dice 'none survives / four survive / three do';")
print("      los tres numeros de arriba tienen que ser 0, 4 y 3)")
assert [n for _, _, n in info] == [0, 4, 3], "*** el conteo NO coincide con el texto del paper ***"
print("     COINCIDE.")
print("     %-22s %d bytes" % ("fig_laplace.pdf", os.path.getsize(os.path.join(OUT, "fig_laplace.pdf"))))
print("DONE")
