# -*- coding: utf-8 -*-
r"""figuras_P3b.py -- las figuras de la Nota II.  La FORMA la elige el trabajo del dato.

Nada de dos curvas y barras apiladas: cada figura usa la forma que corresponde a lo que el
dato tiene que hacer, que es la unica razon legitima para elegir un tipo de grafico.

  fig_umbral.pdf   El dato es un HUECO entre dos cantidades por cada n, no dos series.
                   Forma: GRAFICO DE MANCUERNA, una fila por n, eje logaritmico, el segmento
                   ES el hueco, y la estrella es el primo medido.  El ojo lee la distancia,
                   que es el mensaje, en vez de comparar dos lineas.

  fig_ley.pdf      El dato es un RESTO que debe caber dentro de una cota demostrada.
                   Forma: EMBUDO.  La banda sombreada es el O(p^-1) del teorema y los puntos
                   son la medida; la figura ensena teorema y dato a la vez, y se ve que el
                   dato entra.

  fig_censo.pdf    La afirmacion es "estas dos mitades COINCIDEN".  Forma: PLEGADO Y
                   SUPERPUESTO.  Se dibuja la mitad n <= p/2 rellena y la mitad espejada
                   encima como contorno: si coinciden, la figura DEMUESTRA la simetria en
                   vez de afirmarla.  Al lado, las crestas en 3D, que es donde la tercera
                   dimension si paga porque ensena la subida hacia Delta/b = 1.

PALETA validada con el comprobador de la guia de la casa: las cinco pasan las seis pruebas,
incluida la separacion en deuteranopia y protanopia.  La primera que probe (azul/rojo/verde)
FALLABA con Delta E 5.7 entre verde y rojo: la trampa clasica.

LAS DOS EDICIONES SALEN DE AQUI.  Mismo codigo y mismos datos; lo unico que cambia son las
cadenas, que es la unica forma de que la figura castellana no mienta respecto de la inglesa.
Mismo patron que note_P3a/figuras_P3a.py y note_segmento/figuras.py.

    python figuras_P3b.py [en|es]      ("es" -> fig_x_es.pdf, para P3b_es.tex)

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import io
import os
import re
import sys
from math import factorial

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402
import numpy as np                                                     # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
GATES = os.path.join(os.path.dirname(AQUI), "gates")

IDIOMA = sys.argv[1] if len(sys.argv) > 1 else "en"
if IDIOMA not in ("en", "es"):
    raise SystemExit("idioma desconocido: %s  (usa en o es)" % IDIOMA)


def T(en, es):
    return es if IDIOMA == "es" else en


def num_sec(etiqueta):
    """El numero que LaTeX imprime para esa seccion, leido del .aux.

    Un numero de seccion escrito a mano dentro de una figura es TINTA: no se renumera cuando la
    seccion se mueve, ningun \\ref lo alcanza y el pie de la figura, que si usa \\ref, acaba
    diciendo una cosa distinta de lo que el dibujo pone al lado.  Nos paso con esta misma figura
    al reordenar la nota: el dibujo seguia rotulando la ley como la 5, a Gauss como la 4 y a
    Bring como la 3, que era el orden anterior, y los tres numeros estaban ademas permutados.
    """
    aux = os.path.join(AQUI, "P3b.aux")
    if not os.path.exists(aux):
        raise SystemExit("falta %s: compila la nota una vez antes de generar las figuras" % aux)
    for linea in io.open(aux, encoding="utf-8", errors="replace"):
        m = re.match(r"\\newlabel\{%s\}\{\{(\d+)\}" % re.escape(etiqueta), linea)
        if m:
            return m.group(1)
    raise SystemExit("la etiqueta %s no esta en %s" % (etiqueta, aux))


def guardar(fig, nombre):
    """Anade el sufijo del idioma al nombre, y nada mas."""
    ruta = os.path.join(AQUI, nombre + ("" if IDIOMA == "en" else "_es") + ".pdf")
    fig.savefig(ruta)
    assert os.path.getsize(ruta) > 2000, ruta
    return os.path.basename(ruta)


plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "axes.linewidth": 0.6, "axes.edgecolor": "#8a8a8a",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.color": "#5a5a5a", "ytick.color": "#5a5a5a",
    "legend.frameon": False, "figure.dpi": 160,
    "axes.spines.top": False, "axes.spines.right": False,
})
# paleta validada: node scripts/validate_palette.js -> ALL CHECKS PASS
AZUL, NARANJA, MORADO, VERDE, OCRE = "#1F6FB2", "#D55E00", "#7D3C98", "#117733", "#B8860B"
TINTA, SUAVE = "#333333", "#8a8a8a"


def genero(n):
    return 1 + factorial(n - 2) * (n * n - 5 * n + 2) // 4


# =============================================================== 1. LA MANCUERNA
ns = list(range(5, 10))
inicio = [factorial(n) for n in ns]
fin = [4 * genero(n) ** 2 for n in ns]
# PRIMERAS APARICIONES MEDIDAS, de gates/umbral_real.py, BARRIDO DESDE p = n+1.
# Antes ponia {5: 233, 6: 607} y las dos eran falsas: el 607 era el primer primo de una ventana
# que empezaba en 601, y el 233 se leyo de una tabla con filtro de impresion.  Los de verdad son
# 67, 163 y 601, y caen MUY por debajo de n!, no encima.
medido = {5: 67, 6: 163, 7: 601}

fig, ax = plt.subplots(figsize=(5.8, 2.9))
for i, n in enumerate(ns):
    ax.plot([inicio[i], fin[i]], [i, i], "-", color=SUAVE, lw=5, alpha=0.35,
            solid_capstyle="round", zorder=1)
    ax.plot([inicio[i]], [i], "o", color=AZUL, ms=8, zorder=3)
    ax.plot([fin[i]], [i], "o", color=NARANJA, ms=8, zorder=3)
    # la etiqueta va SIEMPRE a la derecha del punto mas a la derecha: en n = 5 la cota es
    # mas barata que la existencia y el orden de los dos puntos se invierte.
    r = fin[i] / inicio[i]
    ax.text(max(inicio[i], fin[i]) * 2.2, i,
            r"$\times%s$" % ("%.0f" % r if r >= 1 else "%.1f" % r),
            va="center", fontsize=8, color=TINTA)
    if n in medido:
        ax.plot([medido[n]], [i], "*", color=VERDE, ms=14, zorder=4)
        ax.text(medido[n], i + 0.33, r"$p=%d$" % medido[n], ha="center",
                fontsize=7.6, color=VERDE)
ax.set_xscale("log")
ax.set_yticks(range(len(ns)))
ax.set_yticklabels([r"$n=%d$" % n for n in ns])
ax.set_xlim(40, 4e11)
ax.set_ylim(-0.6, len(ns) - 0.25)
ax.set_xlabel(T("prime", "primo"))
ax.plot([], [], "o", color=AZUL, ms=7,
        label=T(r"$n!$: one expected on average", r"$n!$: una esperada en promedio"))
ax.plot([], [], "o", color=NARANJA, ms=7,
        label=T(r"Hasse--Weil becomes usable, $(2g_n)^2$",
                r"Hasse--Weil sirve ya, $(2g_n)^2$"))
ax.plot([], [], "*", color=VERDE, ms=11, label=T("first one measured", "el primero medido"))
ax.legend(loc="lower right", fontsize=7.6, handletextpad=0.3)
ax.grid(axis="x", alpha=0.16, lw=0.5)
fig.tight_layout()
print("%-18s mancuerna, %d filas" % (guardar(fig, "fig_umbral"), len(ns)))
plt.close(fig)

# =============================================================== 2. LAS DOS VARIABLES
# ESTA FIGURA ESTABA MAL, y el fallo era el de la tabla: dibujaba el LOCUS de momentos
# llamandolo F_8.  Ahora dibuja las DOS, porque la diferencia entre ellas es justo lo que el
# cuarto paso del teorema controla, y ademas la curva -105/p que la explica.  Una figura que
# ensena el residuo y la cantidad que lo predice dice mas que una envolvente ajustada al dato.
SAL_EST = os.path.join(GATES, "ley_2p_estabilizadores_OUT.txt")
texto = io.open(SAL_EST, encoding="utf-8", errors="replace").read()
trozo = texto[texto.index("p     U (locus)"):]
P8, LOC, COR = [], [], []
for linea in trozo.splitlines()[1:]:
    f = linea.split()
    if len(f) < 6 or not f[0].isdigit():
        continue
    P8.append(int(f[0]))
    LOC.append(float(f[4]))
    COR.append(float(f[5]))
assert len(P8) == len(LOC) == len(COR) and len(P8) >= 10, (len(P8), len(LOC), len(COR))

fig, ax = plt.subplots(figsize=(6.2, 3.2))
pg = np.linspace(min(P8) - 1, max(P8) + 4, 300)
K105 = 40320.0 / 384.0                              # 8!/384, el coeficiente de la correccion
ax.plot(pg, -K105 / pg, "--", color=MORADO, lw=1.2,
        label=T(r"$-8!/(384\,p)$, the stabiliser term",
                r"$-8!/(384\,p)$, el termino de estabilizadores"))
ax.axhline(0, color=TINTA, lw=0.9)
ax.plot(P8, LOC, "o", color=AZUL, ms=5.5, zorder=3,
        label=T(r"$\widetilde F_8$: the moment locus", r"$\widetilde F_8$: el locus de momentos"))
ax.plot(P8, COR, "o", color=NARANJA, ms=5.5, zorder=3,
        label=T(r"$F_8$: defect $\geq 2$", r"$F_8$: defecto $\geq 2$"))
ax.set_xlabel(r"$p$")
ax.set_ylabel(T(r"$(\cdot-2/p)\,p^2+1$", r"$(\cdot-2/p)\,p^2+1$"))
ax.set_ylim(-5.0, 2.2)
ax.legend(loc="lower right", fontsize=7.4, handletextpad=0.4)
ax.grid(axis="y", alpha=0.16, lw=0.5)
fig.tight_layout()
print("%-18s dos series, %d primos, locus de %.2f a %.2f y corregido de %.2f a %.2f"
      % (guardar(fig, "fig_ley"), len(P8), LOC[0], LOC[-1], COR[0], COR[-1]))
plt.close(fig)

# =============================================================== 3. PLEGADO + CRESTAS
SALIDA = os.path.join(GATES, "nivel_rango_OUT.txt")
repartos = {}
for linea in io.open(SALIDA, encoding="utf-8", errors="replace"):
    f = linea.split()
    if len(f) < 8 or not f[0].isdigit():
        continue
    d = {}
    for par in f[7:]:
        if ":" in par:
            k, v = par.split(":")
            d[int(k)] = int(v)
    repartos[(int(f[0]), int(f[1]))] = d

P = max(p for p, _ in repartos)
primos = sorted({p for p, _ in repartos})

fig = plt.figure(figsize=(8.6, 3.2))

# --- izquierda: PLEGADO.  Si las dos mitades coinciden, la figura lo demuestra.
ax = fig.add_subplot(121)
mitad = [n for n in sorted(n for q, n in repartos if q == P) if n <= P / 2]
deltas = sorted({k for (q, n), d in repartos.items() if q == P for k in d})
COLOR = {0: "#C9CDD2", 1: AZUL, 2: NARANJA, 3: MORADO, 4: VERDE}
abajo = np.zeros(len(mitad))
for k in deltas:
    alt = np.array([repartos[(P, n)].get(k, 0) for n in mitad], dtype=float)
    ax.bar(mitad, alt, bottom=abajo, width=0.74, color=COLOR.get(k, OCRE),
           edgecolor="white", linewidth=1.0, label=r"$\delta=%d$" % k, zorder=2)
    abajo += alt
espejo = np.array([sum(repartos[(P, P - n)].values()) for n in mitad], dtype=float)
ax.plot(mitad, espejo, "_", color=TINTA, ms=17, mew=1.6, zorder=4,
        label=T(r"total at $p-n$", r"total en $p-n$"))
ax.set_xlabel(T(r"$n$   (folded at $n=p/2$)", r"$n$   (plegado en $n=p/2$)"))
ax.set_ylabel(T("orbits", "orbitas"))
ax.set_xticks(mitad)
ax.legend(fontsize=7.4, ncol=2, loc="upper left")
ax.grid(axis="y", alpha=0.16, lw=0.5)
ax.set_title(T(r"the mirrored totals land on the bars",
               r"los totales espejados caen sobre las barras"), fontsize=9, pad=6)

# --- derecha: las crestas en 3D, donde la tercera dimension si paga
ax3 = fig.add_subplot(122, projection="3d")
TONOS = plt.cm.viridis(np.linspace(0.15, 0.9, len(primos)))
for c, p in zip(TONOS, primos):
    nn = sorted(n for q, n in repartos if q == p)
    dd = [sum(k * v for k, v in repartos[(p, n)].items()) / sum(repartos[(p, n)].values())
          for n in nn]
    ax3.plot(np.array(nn) / p, [p] * len(nn), dd, "-o", color=c, ms=2.2, lw=1.2)
ax3.set_xlabel(r"$n/p$", labelpad=-3, fontsize=8)
ax3.set_ylabel(r"$p$", labelpad=-3, fontsize=8)
ax3.set_zlabel(r"$\Delta/b$", labelpad=-5, fontsize=8)
ax3.set_zlim(0, 1.05)
ax3.tick_params(labelsize=6.5, pad=-1)
ax3.view_init(elev=20, azim=-60)
ax3.set_title(T(r"and the plateau rises towards $1$",
                r"y la meseta sube hacia $1$"), fontsize=9, pad=-2)

fig.tight_layout()
nombre_censo = guardar(fig, "fig_censo")
plt.close(fig)

# control: el plegado solo demuestra algo si de verdad coinciden
fallos = sum(1 for n in mitad
             if sum(repartos[(P, n)].values()) != sum(repartos[(P, P - n)].values()))
print("%-18s plegado a p = %d (%d parejas, %d discrepancias) + crestas de %d primos"
      % (nombre_censo, P, len(mitad), fallos, len(primos)))

# =============================================================== 4. EL HILO DE HISTORIA
# No es adorno: es el MAPA de la nota.  Cada ancla esta atada a la seccion que justifica y
# lleva su color.  Lo que se ve y la prosa no da: el vano de casi un siglo entre los dos
# objetos clasicos que indexan los extremos, y el racimo de los anos sesenta y setenta,
# cuando el recuento del centro ya estaba hecho y catalogado sin que lo supieramos.
#
# Una linea temporal es un documento sin auditar hasta que se le mide el solape de etiquetas.
# El primer intento puso las etiquetas a ojo y la auditoria lo tumbo: 1965/1970/1973 se
# pisaban.  Aqui los carriles se ASIGNAN, y luego se vuelve a medir.
# EL HILO CAMBIO DE TESIS, y la figura con el.  Antes decia "dos objetos clasicos indexan los
# dos extremos del censo", que es FALSO: los semisistemas son una de cuatro familias de defecto
# tres.  Lo que los datos SI sostienen, y se midio en gates/P3b_hilo_historia.py, es que todo
# OBJETO que la nota encuentra es anterior a lo que mide (Verlinde 1988) y todo INSTRUMENTO con
# que lo alcanza es posterior.  La particion es limpia: ultimo objeto 1978, primer instrumento
# 2006.  Asi que la figura tiene ahora DOS CARRILES y un eje vertical que los separa.
HILO = [
    (1786, "Bring", "objeto", r"$X^n+aX+b$"),
    (1801, "Gauss", "objeto", T("half-systems", "semisistemas")),
    (1894, "Dedekind", "objeto", T("the conductor", "el conductor")),
    (1965, "Rosa–Znám", "objeto", T("the count", "el recuento")),
    (1970, "Kunz", "objeto", T("symmetric", "simetrico")),
    (1978, "Odlyzko–Stanley", "objeto", T("power sums", "sumas de potencias")),
    (1988, "Verlinde", "medido", T("the column", "la columna")),
    (2006, "Gepner", "instrumento", T("the discriminant", "el discriminante")),
    (2012, "Li–Wan", "instrumento", T("the sieve", "el tamiz")),
    (2021, "Korchmáros", "instrumento", T("the genus", "el genero")),
    (2026, T("this note", "esta nota"), "aqui", ""),
]
CSEC = {"objeto": MORADO, "instrumento": AZUL, "medido": NARANJA, "aqui": VERDE}
ANCHO = 26.0          # anchura tipica de una etiqueta, en anos del eje

# asignacion de carriles: el primero libre en el que no se solape horizontalmente
carril, ocupados = {}, {}
for ano, quien, sec, que in HILO:
    k = 0
    while any(abs(ano - o) < ANCHO for o in ocupados.get(k, [])):
        k += 1
    carril[ano] = k
    ocupados.setdefault(k, []).append(ano)

fig, ax = plt.subplots(figsize=(7.4, 2.7))
ax.axhline(0, color=SUAVE, lw=1.0, zorder=1)
for ano, quien, sec, que in HILO:
    y = 0.30 + 0.60 * carril[ano]
    ax.plot([ano, ano], [0.04, y - 0.04], "-", color=CSEC[sec], lw=0.8, alpha=0.55, zorder=2)
    marca = "*" if sec == "aqui" else "o"
    ax.plot([ano], [0], marca, color=CSEC[sec], ms=15 if sec == "aqui" else 6.5, zorder=3)
    etiqueta = "%s %d" % (quien, ano) if not que else "%s %d\n%s" % (quien, ano, que)
    ax.text(ano, y, etiqueta, ha="center", va="bottom", fontsize=6.9,
            color=VERDE if sec == "aqui" else TINTA, linespacing=1.3)

# el eje que parte la historia: lo que se mide cae entre los objetos y los instrumentos
# la linea sube SOLO por encima del eje: si cruza abajo, se come el rotulo de los 28 anos
ax.vlines(1988, 0.0, 0.30 + 0.60 * (max(carril.values()) + 0.9),
          color=NARANJA, lw=0.9, ls=(0, (3, 2)), alpha=0.75, zorder=1)
ax.annotate("", xy=(2006, -0.30), xytext=(1978, -0.30),
            arrowprops=dict(arrowstyle="<->", color=SUAVE, lw=0.8))
ax.text(1992, -0.38, T("28 years, and nothing in between",
                       "28 anos, y nada en medio"),
        ha="center", va="top", fontsize=7, color=SUAVE)
for k in ("objeto", "medido", "instrumento"):
    ax.plot([], [], "o", color=CSEC[k], ms=6,
            label={"objeto": T("objects it runs into", "objetos que encuentra"),
                   "medido": T("what is measured", "lo que se mide"),
                   "instrumento": T("instruments that reach them",
                                    "instrumentos que los alcanzan")}[k])
ax.legend(loc="upper left", fontsize=7.2, ncol=3, handletextpad=0.3, columnspacing=1.2)
ax.set_xlim(1762, 2050)
ax.set_ylim(-0.62, 0.30 + 0.60 * (max(carril.values()) + 1.35))
ax.axis("off")
fig.tight_layout()
nombre_hilo = guardar(fig, "fig_hilo")

fig.canvas.draw()
ren = fig.canvas.get_renderer()
cajas = [(t.get_text().split("\n")[0], t.get_window_extent(renderer=ren)) for t in ax.texts]
choques = [(cajas[i][0], cajas[j][0])
           for i in range(len(cajas)) for j in range(i + 1, len(cajas))
           if cajas[i][1].x1 > cajas[j][1].x0 and cajas[j][1].x1 > cajas[i][1].x0
           and cajas[i][1].y1 > cajas[j][1].y0 and cajas[j][1].y1 > cajas[i][1].y0]
plt.close(fig)
if choques:
    raise SystemExit("ETIQUETAS QUE SE PISAN (%s): %s" % (IDIOMA, choques))
print("%-18s %d anclas en %d carriles, 0 solapes de etiqueta"
      % (nombre_hilo, len(HILO), max(carril.values()) + 1))

# =============================================================== 5. EL MAPA, SOBRE EL CENSO
# Un indice en forma de tabla es lo de siempre.  El censo tiene un eje propio --- el defecto
# delta --- y cada seccion de la nota es DUENA DE UN TRAMO de ese eje.  Asi que el mapa se
# dibuja encima del objeto: el lector ve de una vez el reparto real y donde vive cada seccion.
# El tramo del extremo se dibuja VACIO, porque a este primo lo esta: no aparece hasta p ~ n!.
agreg = {}
for (p, n), d in repartos.items():
    if p != P:
        continue
    for k, v in d.items():
        agreg[k] = agreg.get(k, 0) + v
dmax = max(agreg) + 1                       # el estrato extremo, vacio a este primo
ks = list(range(dmax + 1))
alturas = [agreg.get(k, 0) for k in ks]

fig, ax = plt.subplots(figsize=(7.2, 2.8))
for k, h in zip(ks, alturas):
    if h:
        ax.bar([k], [h], width=0.62, color=AZUL if k == 1 else (
            "#C9CDD2" if k == 0 else MORADO), edgecolor="white", lw=1.0, zorder=3)
        # el numero va DENTRO de la barra si es alta: fuera chocaba con la etiqueta del tramo
        if h >= 100:
            ax.text(k, h * 0.55, "%d" % h, ha="center", va="center", fontsize=7.4,
                    color="white", zorder=4)
        else:
            ax.text(k, h * 1.35, "%d" % h, ha="center", fontsize=7.4, color=TINTA)
    else:
        ax.bar([k], [1], width=0.62, color="none", edgecolor=NARANJA, lw=1.1,
               linestyle=(0, (2, 1.6)), zorder=3)
        # NO es "el" defecto maximo: el maximo es C(n-1,2) y DEPENDE DE n, asi que agregando
        # sobre n no hay una sola casilla que lo represente.  La columna se dibuja separada y
        # rotulada como lo que es: un estrato que a este primo esta vacio para todo n.
        ax.text(k, 1.6, T("$\\delta=\\binom{n-1}{2}$\nvaries with $n$;\nempty at this prime",
                          "$\\delta=\\binom{n-1}{2}$\ndepende de $n$;\nvacio a este primo"),
                ha="center", fontsize=6.2, color=NARANJA, style="italic", linespacing=1.25)
ax.set_yscale("log")
ax.set_ylim(0.7, 1.2e5)   # techo alto: los corchetes de tramo deben DESPEJAR la barra mas alta
ax.set_xlim(-0.75, dmax + 0.75)
ax.set_xticks(ks)
ax.set_xlabel(T(r"defect $\delta$", r"defecto $\delta$"))
ax.set_ylabel(T("orbits", "orbitas"))
ax.grid(axis="y", alpha=0.14, lw=0.5)

# los tramos: cada seccion es duena de un trozo del eje
TRAMOS = [(0, 0, "#7a7a7a", T("maximal orders\n(Part I)", "ordenes maximales\n(Parte I)"), 0),
          (1, 1, AZUL, T("almost everything\n§%s: $2/p-1/p^2$" % num_sec("sec:ley"),
                         "casi todo\n§%s: $2/p-1/p^2$" % num_sec("sec:ley")), 0),
          (2, 3, MORADO, T("the centre\n§%s: Gauss's half-systems" % num_sec("sec:centro"),
                           "el centro\n§%s: semisistemas de Gauss" % num_sec("sec:centro")), 0),
          (dmax, dmax, NARANJA,
           T("the extreme\n§%s: Bring, from $p\\gtrsim n!$" % num_sec("sec:bring"),
             "el extremo\n§%s: Bring, desde $p\\gtrsim n!$" % num_sec("sec:bring")), 0)]
for a, b, c, txt, _ in TRAMOS:
    ax.annotate("", xy=(a - 0.34, 0.80), xytext=(b + 0.34, 0.80),
                xycoords=("data", "axes fraction"), textcoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="|-|,widthA=0.28,widthB=0.28", color=c, lw=1.0))
    ax.text((a + b) / 2, 0.86, txt, transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=7.0, color=c, linespacing=1.3)
ax.set_title(T(r"the whole census at $p=%d$, and which section owns which part" % P,
               r"el censo entero a $p=%d$, y de que seccion es cada tramo" % P),
             fontsize=9.5, pad=22)
fig.tight_layout()
nombre_mapa = guardar(fig, "fig_mapa")

# misma auditoria de solapes que la linea temporal: una figura con etiquetas es un
# documento sin auditar hasta que se le mide el choque.
fig.canvas.draw()
ren = fig.canvas.get_renderer()
cajas = [(t.get_text().split("\n")[0], t.get_window_extent(renderer=ren))
         for t in ax.texts if t.get_text().strip()]
choques = [(cajas[i][0], cajas[j][0])
           for i in range(len(cajas)) for j in range(i + 1, len(cajas))
           if cajas[i][1].x1 > cajas[j][1].x0 and cajas[j][1].x1 > cajas[i][1].x0
           and cajas[i][1].y1 > cajas[j][1].y0 and cajas[j][1].y1 > cajas[i][1].y0]
plt.close(fig)
if choques:
    raise SystemExit("MAPA, etiquetas que se pisan (%s): %s" % (IDIOMA, choques))
print("%-18s eje delta 0..%d, %d orbitas, el extremo vacio a p = %d"
      % (nombre_mapa, dmax, sum(alturas), P))


# =============================================================== 6. LA ESCALERA DE LOS DEFECTOS
# POR QUE ESTA FIGURA.  Es el corazon nuevo de la nota y no tenia dibujo.  Y tiene una propiedad
# que una tabla no puede ensenar: EL MISMO PEINE SE LEE DOS VECES.  Las casillas llenas son a la
# vez los HUECOS del semigrupo y los INDICES de los momentos que se anulan --- que es exactamente
# lo que dice el lema.  Ademas el numero de filas de cada banda ES el coeficiente de la
# expansion: una, dos, cuatro.  El lector ve de donde salen el 2 de 2/p y el 4 de 4/p^2, y de ahi
# el 3 = -1 + 4 del defecto medio.  Los recuentos salen de gates/P3b_escalera_defectos.py.
ESCALERA = [
    (1, AZUL, [("2,3", (1,), 12228)]),
    (2, MORADO, [("3,4,5", (1, 2), 451), ("2,5", (1, 3), 403)]),
    (3, NARANJA, [("4,5,6,7", (1, 2, 3), 16), ("3,5,7", (1, 2, 4), 19),
                  ("3,4", (1, 2, 5), 27), ("2,7", (1, 3, 5), 2)]),
]
TERMINO = {1: T(r"$1-O(p^{-1})$", r"$1-O(p^{-1})$"),
           2: r"$2/p$", 3: r"$4/p^{2}$"}
NMAX = 8

fig, ax = plt.subplots(figsize=(7.6, 4.3))
y = 0.0
filas = []
for delta, color, familias in ESCALERA:
    y0 = y
    for nombre, huecos, medidas in familias:
        for k in range(NMAX + 1):
            hueco = k in huecos
            ax.add_patch(plt.Rectangle((k - 0.40, y - 0.30), 0.80, 0.60,
                                       facecolor=color if hueco else "none",
                                       edgecolor=color if hueco else "#D8DCE0",
                                       linewidth=0.8, zorder=2))
            if hueco:
                # la casilla llena, leida por segunda vez: el momento que se anula
                ax.text(k, y - 0.52, r"$M_{%d}$" % k, ha="center", va="top",
                        fontsize=6.6, color=color)
        ax.text(-1.5, y, r"$\langle%s\rangle$" % nombre, ha="right", va="center",
                fontsize=8.6, color=TINTA)
        ax.text(NMAX + 1.9, y, "%d" % medidas, ha="right", va="center",
                fontsize=7.4, color=SUAVE)
        filas.append(y)
        y -= 1.12
    # el margen derecho: cuantas familias, que es el coeficiente
    ym = (y0 + y + 1.12) / 2
    ax.text(NMAX + 2.6, ym, r"$\times%d$" % len(familias), ha="center", va="center",
            fontsize=15, color=color)
    ax.text(NMAX + 2.6, ym - 0.42, TERMINO[delta], ha="center", va="top",
            fontsize=8, color=color)
    ax.text(-3.5, ym, r"$\delta=%d$" % delta, ha="center", va="center",
            fontsize=11, color=color)
    y -= 0.55

ax.text(NMAX + 1.9, filas[0] + 0.92, T("orbits", "orbitas"), ha="right", va="bottom",
        fontsize=6.8, color=SUAVE, style="italic")
ax.text((NMAX) / 2.0, filas[0] + 0.92,
        T("filled = a gap of $S_T$ = a vanishing moment",
          "llena = un hueco de $S_T$ = un momento nulo"),
        ha="center", va="bottom", fontsize=8, color=TINTA)
# el remate: de donde sale el coeficiente del defecto medio.  El -1 es de la expansion de
# F_n y el 4 es el numero de filas de la banda de abajo, que se esta viendo.
ax.text((NMAX) / 2.0, y + 0.10,
        T(r"the mean defect per orbit is $1+2/p+3/p^{2}$,   with   $3=-1+4$",
          r"el defecto medio por orbita es $1+2/p+3/p^{2}$,   con   $3=-1+4$"),
        ha="center", va="top", fontsize=9.2, color=VERDE)
ax.set_xlim(-4.6, NMAX + 4.4)
ax.set_ylim(y - 0.75, filas[0] + 1.5)
ax.axis("off")
fig.tight_layout()
nombre_esc = guardar(fig, "fig_escalera")

# la misma auditoria de solapes que las otras dos figuras con etiquetas
fig.canvas.draw()
ren = fig.canvas.get_renderer()
cajas = [(t.get_text(), t.get_window_extent(renderer=ren)) for t in ax.texts if t.get_text().strip()]
choques = [(cajas[i][0], cajas[j][0])
           for i in range(len(cajas)) for j in range(i + 1, len(cajas))
           if cajas[i][1].x1 > cajas[j][1].x0 and cajas[j][1].x1 > cajas[i][1].x0
           and cajas[i][1].y1 > cajas[j][1].y0 and cajas[j][1].y1 > cajas[i][1].y0]
plt.close(fig)
if choques:
    raise SystemExit("ESCALERA, etiquetas que se pisan (%s): %s" % (IDIOMA, choques[:4]))
print("%-18s escalera 1-2-4, %d familias, %d orbitas medidas"
      % (nombre_esc, sum(len(f) for _, _, f in ESCALERA),
         sum(m for _, _, f in ESCALERA for _, _, m in f)))
