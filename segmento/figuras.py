# -*- coding: utf-8 -*-
"""Figuras de la nota del segmento.

fig_plano.pdf -- LA LEY ES UN PLANO.  Cada casilla medida es un punto (W_1, dim de capa, v_p), y
la ley  v_p = 2 dim - 1 - W_1  es un plano.  El caso REAL vive con dim = phi(q')/2 y la SOMBRA con
dim = phi(q') entera: si el eje es la dimension de la capa, las dos poblaciones caen en el MISMO
plano, que es la tesis hecha imagen -- el grupo no importa, importa el segmento.  Lo que se despega
no se oculta, y se separa en sus dos causas: la zona de la escalera (W_1 = 0) y el resto de las
casillas fuera de la hipotesis.

MUCHAS casillas comparten posicion: se agregan y el tamano del punto dice cuantas hay debajo, que
es lo honesto -- dispersarlas con ruido seria inventar datos.

Los datos salen de gates/capas_tipo/plano_denso_OUT.txt, no se escriben a mano.

Authors: Carles Marin, Claude (AI assistant).
Run:  python figuras.py
"""
import io
import math
import sys
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import numpy as np

AQUI = pathlib.Path(__file__).resolve().parent
# Los datos de las figuras viven en gates/capas_tipo/ del arbol de trabajo; en el paquete auxiliar
# viajan AL LADO de este guion, que es plano y no puede leer fuera de si mismo.
G = AQUI if (AQUI / "plano_denso_OUT.txt").exists() else AQUI.parent / "gates" / "capas_tipo"

AZUL = "#1F4E79"
NARANJA = "#C2610A"
ROJO = "#B03030"
GRIS = "#7A7A7A"

# "en" -> fig_x.pdf ;  "es" -> fig_x_es.pdf.  Mismo patron que note_conductor/figuras.py.
#
# LAS FIGURAS SE TRADUCEN PORQUE SON TEXTO.  Un pie de figura en castellano sobre un dibujo cuyos
# rotulos estan en ingles deja al lector de la version castellana leyendo media figura, y es la
# mitad que lleva los numeros.  Las dos ediciones se generan del MISMO codigo, con los mismos
# datos: lo unico que cambia son las cadenas, asi que una figura no puede decir cosas distintas en
# un idioma y en otro.
IDIOMA = "en"


def T(en, es):
    return es if IDIOMA == "es" else en


def guardar(fig, nombre):
    """Anade el sufijo del idioma al nombre, y nada mas."""
    out = AQUI / (nombre + ("" if IDIOMA == "en" else "_" + IDIOMA) + ".pdf")
    fig.savefig(out, bbox_inches="tight")
    return out.name


def lee():
    """(tipo, dim, W_1, v_p) del barrido ancho."""
    out = []
    for ln in io.open(G / "plano_denso_OUT.txt", encoding="utf-8", errors="replace"):
        if not ln.strip() or ln.startswith("#"):
            continue
        c = ln.split()
        if len(c) != 8:
            continue
        out.append((c[0], int(c[5]), int(c[6]), int(c[7])))
    return out


def fig_plano():
    d = lee()
    enp = [x for x in d if x[3] == 2 * x[1] - 1 - x[2]]
    fue = [x for x in d if x[3] != 2 * x[1] - 1 - x[2]]
    escal = [x for x in fue if x[2] == 0]
    otros = [x for x in fue if x[2] != 0]

    fig = plt.figure(figsize=(9.6, 4.6))
    ax = fig.add_subplot(121, projection="3d")
    # vista DE CANTO: la normal del plano es (1,-2,1), y mirar en la direccion (1,1,1) -- que es
    # perpendicular a ella -- lo proyecta sobre una RECTA.  Lo que esta en el plano se alinea; lo
    # que no, se sale.  Es la prueba visual, no el adorno.
    ax2 = fig.add_subplot(122, projection="3d")

    wmax = max(x[2] for x in d) + 1
    dmax = max(x[1] for x in d) + 1
    W, D = np.meshgrid(np.linspace(0, wmax, 30), np.linspace(1, dmax, 30))
    V = 2 * D - 1 - W
    V[V < 0] = np.nan
    for a in (ax, ax2):
        a.plot_surface(W, D, V, alpha=0.15, color=AZUL, linewidth=0, antialiased=True,
                       rstride=1, cstride=1)

    def agrega(lista):
        c = {}
        for x in lista:
            c[(x[2], x[1], x[3])] = c.get((x[2], x[1], x[3]), 0) + 1
        return c

    grupos = [("real", AZUL, "o", T("real $\\mathrm{SU}(n)$", "real $\\mathrm{SU}(n)$"), enp),
              ("sombra", NARANJA, "^", T("shadow $\\mathrm{GL}(n)$", "sombra $\\mathrm{GL}(n)$"), enp),
              (None, ROJO, "X", T("ladder zone, $W_1=0$", "zona de la escalera, $W_1=0$"), escal),
              (None, GRIS, "s", T("other cells off the hypothesis",
                                 "otras casillas fuera de la hip\u00f3tesis"), otros)]

    for eje in (ax, ax2):
        for tipo, color, marca, nombre, fuente in grupos:
            lista = [x for x in fuente if (tipo is None or x[0] == tipo)]
            if not lista:
                continue
            a = agrega(lista)
            et = None
            if eje is ax:
                if tipo is not None:
                    et = T("%s: %d cells, %d sites", "%s: %d casillas, %d sitios") % (
                        nombre, sum(a.values()), len(a))
                else:
                    et = "%s (%d)" % (nombre, sum(a.values()))
            base = 16 if tipo is not None else 40
            eje.scatter([k[0] for k in a], [k[1] for k in a], [k[2] for k in a],
                        c=color, marker=marca,
                        s=[base + 20 * (v ** 0.5) for v in a.values()],
                        depthshade=False, edgecolors="white", linewidths=0.5, alpha=0.95,
                        label=et)

    # la elevacion de canto, calculada y no estimada: ver el docstring del parche.
    AZIM2 = 45.0
    def elev_de_canto(eje, azim, aspecto):
        (x0, x1), (y0, y1), (z0, z1) = eje.get_xlim3d(), eje.get_ylim3d(), eje.get_zlim3d()
        bx, by, bz = aspecto
        sx, sy, sz = (x1 - x0) / bx, (y1 - y0) / by, (z1 - z0) / bz
        nx, ny, nz = sx, -2.0 * sy, sz
        a = np.radians(azim)
        return np.degrees(np.arctan2(-(nx * np.cos(a) + ny * np.sin(a)), nz))

    ASPECTO = (1.2, 1.0, 0.85)
    ax2.set_box_aspect(ASPECTO)
    ELEV2 = elev_de_canto(ax2, AZIM2, ASPECTO)
    print("              | elevacion de canto calculada: %.2f grados" % ELEV2)

    for eje, titulo, vista in (
            (ax, T("(a) the measured cells and the plane",
                    "(a) las casillas medidas y el plano"), (18, -58)),
            (ax2, T("(b) edge-on: the plane becomes a line",
                    "(b) de canto: el plano se vuelve una recta"), (ELEV2, AZIM2))):
        eje.set_xlabel("$W_1$", labelpad=3)
        eje.set_ylabel(T("$d$ or $D$", "$d$ o $D$"), labelpad=5)
        # la etiqueta z solo en el panel (a): en (b) cae justo donde la de (a) y se pisan
        if eje is ax:
            eje.set_zlabel("$v_p$", labelpad=6)
        eje.set_title(titulo, pad=2, fontsize=9.5)
        eje.view_init(elev=vista[0], azim=vista[1])
        eje.set_box_aspect((1.2, 1.0, 0.85))
        eje.tick_params(labelsize=7.5, pad=1)
        for a in (eje.xaxis, eje.yaxis, eje.zaxis):
            a.pane.set_alpha(0.06)

    fig.suptitle(T("$v_p=2\\,\\dim-1-W_1$: one plane, two rings",
                 "$v_p=2\\,\\dim-1-W_1$: un plano, dos anillos"), y=0.99, fontsize=11)
    ax.legend(loc="upper left", fontsize=7.2, framealpha=0.93, borderpad=0.4,
              bbox_to_anchor=(-0.10, 1.00))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.02, wspace=0.10)
    guardar(fig, "fig_plano")
    plt.close(fig)
    bajo = [x for x in d if 2 * x[2] > x[1]]
    bajo_ok = [x for x in bajo if x[3] == 2 * x[1] - 1 - x[2]]
    print("fig_plano.pdf | %d filas, %d en el plano, %d fuera (%d escalera, %d otras)" % (
        len(d), len(enp), len(fue), len(escal), len(otros)))
    print("              | bajo W_1 > dim/2: %d de %d" % (len(bajo_ok), len(bajo)))




def fig_formas():
    """Las tres formas de segmento son EL MISMO BLOQUE en tres posiciones.

    Lo que cambia no es el conjunto sino donde cae respecto del eje de simetria c -> -c:
    centrado en 0 (real, diente de sierra), centrado en q'/2 (real, diente desplazado), y sin
    centrar (no real, la sombra).  Esa es la tesis del paper en dos centimetros de papel.
    """
    import matplotlib.patches as mpatches
    QP = 12          # clases modulo q'
    LARGO = 9        # el bloque tiene el mismo tamano en los tres

    casos = [
        (0.0, T("centred at $0$", "centrado en $0$"),
         "$\\mathrm{Sp}(2m)$, $\\mathrm{SU}(2m{+}1)$",
         T("Stickelberger's sawtooth", "el diente de sierra de Stickelberger"), AZUL),
        (QP / 2.0, T("centred at $q'/2$", "centrado en $q'/2$"), "$\\mathrm{SU}(2m)$",
         T("the half-period shift", "el desplazamiento de medio per\u00edodo"), "#4E8C4A"),
        ((LARGO - 1) / 2.0, T("not centred", "sin centrar"), "$\\mathrm{GL}(n)$",
         T("no symmetry: the shadow", "sin simetr\u00eda: la sombra"), NARANJA),
    ]

    # La altura sube de 3.25 a 3.96 para dejar sitio a la segunda flecha SIN encoger los paneles:
    # con bottom=0.24 la caja de ejes vale 0.64 de la figura, y 2.535/0.64 = 3.96 devuelve a los
    # circulos su tamano absoluto de antes.  Si se deja 3.25, los circulos encogen, los rotulos no,
    # y "complement, q odd" se come el "c -> -c" del panel izquierdo.
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 3.96))
    for eje, (centro, titulo, grupo, gobierna, color) in zip(axes, casos):
        eje.set_aspect("equal")
        eje.axis("off")
        ang = lambda c: 2 * np.pi * c / QP
        th = np.linspace(0, 2 * np.pi, 200)
        eje.plot(np.cos(th), np.sin(th), "-", lw=0.7, color="#DCE2E8", zorder=0)
        # todas las clases, en gris
        for c in range(QP):
            eje.plot(np.cos(ang(c)), np.sin(ang(c)), "o", ms=4.5, color="#C3CBD4", zorder=2)
        # el eje de simetria c -> -c
        eje.plot([-1.22, 1.22], [0, 0], "-", lw=0.8, color="#9AA7B4", zorder=1)
        eje.text(1.26, 0.03, "$c\\mapsto-c$", fontsize=7, color="#6B7885", va="bottom")
        # el bloque
        bloque_cl = [(round(centro - (LARGO - 1) / 2.0) + k) % QP for k in range(LARGO)]
        xs = [np.cos(ang(c)) for c in bloque_cl]
        ys = [np.sin(ang(c)) for c in bloque_cl]
        eje.plot(xs, ys, "o", ms=8, color=color, zorder=4, markeredgecolor="white", mew=0.8)
        # EL REFLEJO c -> -c, en aros huecos: si el segmento es simetrico los aros caen justo
        # encima de los puntos; si no lo es, se despegan.  Eso es lo que separa los tres casos.
        refl = [(-c) % QP for c in bloque_cl]
        eje.plot([np.cos(ang(c)) for c in refl], [np.sin(ang(c)) for c in refl], "o",
                 ms=13, mfc="none", mec=color, mew=1.1, zorder=3, alpha=0.75)
        # EL EJE PROPIO DEL BLOQUE: la recta por su punto medio.  Un marcador puntual sobre el
        # circulo caia encima de un dato y no se entendia; el eje, en cambio, dice lo que importa:
        # en los dos casos centrados coincide con el eje gris c -> -c, y en el tercero no.
        ejex = [-0.82 * np.cos(ang(centro)), 0.82 * np.cos(ang(centro))]
        ejey = [-0.82 * np.sin(ang(centro)), 0.82 * np.sin(ang(centro))]
        eje.plot(ejex, ejey, "--", lw=1.3, color=color, zorder=2, dashes=(5, 3))
        eje.set_xlim(-1.45, 1.45)
        # las etiquetas VAN FUERA del circulo: a -1.16 se comian los puntos de abajo
        eje.set_ylim(-1.72, 1.58)
        eje.text(0, 1.42, titulo, ha="center", fontsize=10, color=color)
        if centro not in (0.0, QP / 2.0):
            eje.text(0, -1.17, T("its axis misses $c\\mapsto-c$",
                                 "su eje no es $c\\mapsto-c$"), ha="center", fontsize=7.5,
                     color=color)
        eje.text(0, -1.40, grupo, ha="center", fontsize=9)
        eje.text(0, -1.60, gobierna, ha="center", fontsize=8, color="#55606C", style="italic")
    # LA FLECHA DEL COMPLEMENTO.  Para q impar el complemento n -> q-n lleva el caso par al caso
    # centrado: el panel del medio deja de ser un regimen aparte y cae dentro del primero.  Va
    # entre los dos paneles, en coordenadas de figura, y por debajo de los circulos para no tapar
    # ningun dato.
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.24, wspace=0.05)
    izq = axes[0].get_position()
    med = axes[1].get_position()
    der = axes[2].get_position()
    y0 = izq.y0 + 0.30 * izq.height
    fig.patches.append(matplotlib.patches.FancyArrowPatch(
        (med.x0 + 0.16 * med.width, y0), (izq.x1 - 0.16 * izq.width, y0),
        transform=fig.transFigure, arrowstyle="->", mutation_scale=11, lw=1.0,
        color="#4E8C4A", connectionstyle="arc3,rad=0.28", zorder=10))
    fig.text((izq.x1 + med.x0) / 2.0, y0 + 0.085,
             T("complement,\n$q$ odd", "complemento,\n$q$ impar"), ha="center",
             va="bottom", fontsize=7.5, color="#4E8C4A")
    # LA SEGUNDA FLECHA.  El complemento hace del panel del medio un caso del primero: es una
    # IGUALDAD de anillos.  El tercer panel no se reduce por igualdad --- el anillo de GL(n) no es
    # real --- pero tampoco es un regimen aparte: es una extension libre de GRADO 2 del anillo del
    # primer panel, el de Sp(2(r-1)).  Va por DEBAJO de los tres ejes, no dentro, para no cruzarse
    # con la verde ni pisar los rotulos; por eso el bottom sube de 0.10 a 0.24.
    yb = izq.y0 - 0.085
    fig.patches.append(matplotlib.patches.FancyArrowPatch(
        (der.x0 + 0.30 * der.width, yb), (izq.x1 - 0.30 * izq.width, yb),
        transform=fig.transFigure, arrowstyle="->", mutation_scale=11, lw=1.0,
        color=NARANJA, connectionstyle="arc3,rad=0.06", zorder=10))
    fig.text((izq.x1 + der.x0) / 2.0, yb - 0.095,
             "free of rank $2$ over it: $B=R\\oplus\\theta R$,\n"
             "$p$ odd", ha="center", va="bottom", fontsize=7.5,
             color=NARANJA)
    fig.suptitle(T("one block, three positions: the conductor depends on where the segment sits",
                   "un bloque, tres posiciones: el conductor depende de d\u00f3nde cae el segmento"),
                 y=1.00, fontsize=10.5)
    guardar(fig, "fig_formas")
    plt.close(fig)
    print("fig_formas.pdf | tres posiciones del mismo bloque de %d clases modulo %d" % (LARGO, QP))




def fig_operaciones():
    """LAS TRES OPERACIONES QUE NO CAMBIAN EL ANILLO, en un solo circulo.

    El espectro de SU(n) en g_q son n de los 2q exponentes de UNA clase de paridad modulo 2q:
    n-1, n-3, ..., 1-n.  Sobre ese circulo se leen las tres operaciones del Lema de division, y
    las tres son el mismo hecho -- multiplicar por un factor monico de Z[T] no cambia el anillo de
    coeficientes:

      (a) anadir el autovalor 1        el factor es  T - 1        SU(2m+1) = Sp(2m)
      (b) tomar el complemento         el producto es T^q -+ 1    SU(n) = SU(q-n)
      (c) dar una vuelta entera        el factor es  T^q -+ 1     SU(n+q) = SU(n)

    Nada aqui esta dibujado a ojo: las posiciones son las de los exponentes, calculadas.
    """
    Q = 9
    DOS = 2 * Q

    def pos(e):
        a = 2 * np.pi * (e % DOS) / DOS
        return np.cos(a), np.sin(a)

    def exps(n):
        return [(n - 1 - 2 * j) % DOS for j in range(n)]

    fig, axes = plt.subplots(1, 3, figsize=(8.8, 3.5))
    for eje in axes:
        eje.set_aspect("equal")
        eje.axis("off")
        th = np.linspace(0, 2 * np.pi, 240)
        eje.plot(np.cos(th), np.sin(th), "-", lw=0.7, color="#DCE2E8", zorder=0)
        for e in range(DOS):
            x, y = pos(e)
            eje.plot([x], [y], "o", ms=3.2, color="#C3CBD4", zorder=1)
        eje.set_xlim(-1.5, 1.5)
        eje.set_ylim(-1.62, 1.62)

    # (a) anadir el autovalor 1:  SU(9) = Sp(8) + el punto 0
    ea = exps(Q)                       # n = q = 9 impar: los exponentes PARES
    for e in ea:
        x, y = pos(e)
        col = GRIS if e == 0 else AZUL
        axes[0].plot([x], [y], "o", ms=9, color=col, mec="white", mew=0.9, zorder=3)
    x, y = pos(0)
    axes[0].annotate("$x=1$", (x, y), textcoords="offset points", xytext=(16, 10),
                     fontsize=8, color=GRIS,
                     arrowprops=dict(arrowstyle="-", lw=0.7, color=GRIS))
    axes[0].set_title(T("(a) adjoin the eigenvalue $1$", "(a) a\u00f1adir el autovalor $1$"),
                      fontsize=9.5, pad=6)
    axes[0].text(0, -1.34, "$\\mathrm{SU}(9)=\\mathrm{Sp}(8)$", ha="center", fontsize=9)
    axes[0].text(0, -1.56, T("the factor is $T-1$", "el factor es $T-1$"),
                 ha="center", fontsize=8, color="#55606C",
                 style="italic")

    # (b) complemento:  SU(4) y lo que queda, que es -(espectro de SU(5))
    n = 4
    eb = exps(n)
    resto = [e for e in range(DOS) if e % 2 == eb[0] % 2 and e not in eb]
    VERDE = "#4E8C4A"                  # el mismo verde que marca SU(2m) en fig_formas
    for e in eb:
        x, y = pos(e)
        axes[1].plot([x], [y], "o", ms=9, color=AZUL, mec="white", mew=0.9, zorder=4)
    for e in resto:
        x, y = pos(e)
        axes[1].plot([x], [y], "o", ms=9, color=VERDE, mec="white", mew=0.9, zorder=3)
        xx, yy = pos(e + Q)            # x -> -x  es girar media vuelta: multiplicar por t^q
        axes[1].plot([xx], [yy], "o", ms=13, mfc="none", mec=VERDE, mew=1.2, zorder=2,
                     alpha=0.85)
    axes[1].annotate("", xy=pos(resto[0] + Q), xytext=pos(resto[0]),
                     arrowprops=dict(arrowstyle="->", lw=0.9, color=VERDE,
                                     connectionstyle="arc3,rad=0.25"), zorder=5)
    axes[1].text(-0.06, 0.02, "$\\times\\,t^{q}$", ha="right", va="center", fontsize=8,
                 color=VERDE)
    # sin estas dos etiquetas el lector ve tres colores y no sabe cual es cual
    xb, yb = pos(eb[0])
    axes[1].annotate("$\\mathrm{SU}(4)$", (xb, yb), textcoords="offset points", xytext=(13, 7),
                     fontsize=8, color=AZUL)
    xr, yr = pos(resto[len(resto) // 2] + Q)
    axes[1].annotate("$\\mathrm{SU}(5)$", (xr, yr), textcoords="offset points", xytext=(14, -4),
                     fontsize=8, color=VERDE)
    # El pie del paper ya dice que este panel es level-rank; el rotulo lo dice tambien, para que
    # la figura no dependa del pie para ser honesta sobre la atribucion.
    axes[1].set_title(T("(b) take the complement: level-rank",
                        "(b) tomar el complemento: nivel-rango"), fontsize=9.5, pad=6)
    axes[1].text(0, -1.34, T("$\\mathrm{SU}(4)=\\mathrm{SU}(5)$ at $q=9$",
                             "$\\mathrm{SU}(4)=\\mathrm{SU}(5)$ en $q=9$"), ha="center", fontsize=9)
    axes[1].text(0, -1.56, T("the product is $T^{q}+1$", "el producto es $T^{q}+1$"),
                 ha="center", fontsize=8, color="#55606C",
                 style="italic")

    # (c) una vuelta entera:  n = 4  y  n = 13  ocupan lo mismo, con 4 dobles
    ec = exps(n + Q)
    mult = {}
    for e in ec:
        mult[e] = mult.get(e, 0) + 1
    for e, k in mult.items():
        x, y = pos(e)
        axes[2].plot([x], [y], "o", ms=9, color=GRIS, mec="white", mew=0.9, zorder=3)
        if k > 1:
            axes[2].plot([x], [y], "o", ms=9, color=AZUL, mec="white", mew=0.9, zorder=4)
            axes[2].plot([x], [y], "o", ms=14, mfc="none", mec=AZUL, mew=1.2, zorder=2)
    # el sobrante no ES SU(4): es SU(4) girada media vuelta, t^q = -1.  Dibujarlo sin ese paso
    # seria una figura que miente, asi que se dibuja el paso.
    dob = sorted(e for e, k in mult.items() if k > 1)
    for e in dob:
        xx, yy = pos(e + Q)
        axes[2].plot([xx], [yy], "o", ms=13, mfc="none", mec=AZUL, mew=1.2, zorder=2, alpha=0.85)
    axes[2].annotate("", xy=pos(dob[0] + Q), xytext=pos(dob[0]),
                     arrowprops=dict(arrowstyle="->", lw=0.9, color=AZUL,
                                     connectionstyle="arc3,rad=0.25"), zorder=5)
    axes[2].text(-0.06, 0.02, "$\\times\\,t^{q}$", ha="right", va="center", fontsize=8,
                 color=AZUL)
    xd, yd = pos(dob[-1])
    axes[2].annotate(T("twice", "dos veces"), (xd, yd), textcoords="offset points", xytext=(-14, 12),
                     fontsize=8, color=AZUL, ha="right")
    xs4, ys4 = pos(dob[len(dob) // 2] + Q)
    axes[2].annotate("$\\mathrm{SU}(4)$", (xs4, ys4), textcoords="offset points", xytext=(14, -4),
                     fontsize=8, color=AZUL)
    axes[2].set_title(T("(c) one full turn", "(c) una vuelta entera"), fontsize=9.5, pad=6)
    axes[2].text(0, -1.34, T("$\\mathrm{SU}(13)=\\mathrm{SU}(4)$ at $q=9$",
                             "$\\mathrm{SU}(13)=\\mathrm{SU}(4)$ en $q=9$"), ha="center", fontsize=9)
    axes[2].text(0, -1.56, T("the factor is $T^{q}-1$", "el factor es $T^{q}-1$"),
                 ha="center", fontsize=8, color="#55606C",
                 style="italic")
    dobles = len(dob)

    fig.suptitle(T("three operations on the segment, one division lemma: the ring does not change",
                   "tres operaciones sobre el segmento, un lema de divisi\u00f3n: el anillo no cambia"),
                 y=1.02, fontsize=10.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.06, wspace=0.05)
    guardar(fig, "fig_operaciones")
    plt.close(fig)
    print("fig_operaciones.pdf | q=%d: (a) %d exponentes, (b) %d+%d, (c) %d dobles" % (
        Q, len(ea), len(eb), len(resto), dobles))


def kron(a, n):
    """simbolo de Kronecker (a|n), suficiente para n impar libre de cuadrados."""
    if n == 1:
        return 1
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def orden_mayor():
    """Las casillas (f, m) donde muere un caracter impar de ORDEN > 2.

    Salen de gates/capas_tipo/caracteres_orden_mayor_OUT.txt, que las calcula en aritmetica
    entera exacta y las controla contra el simbolo de Kronecker en el caso de orden 2.  No se
    escriben a mano: si la compuerta cambia, la figura cambia.
    """
    import ast
    marca = "casillas (f, m, orden) con orden > 2:"
    for ln in io.open(G / "caracteres_orden_mayor_OUT.txt", encoding="utf-8",
                      errors="replace"):
        if ln.startswith(marca):
            return ast.literal_eval(ln[len(marca):].strip())
    return []


def fig_ceros():
    """Donde muere una componente impar: una familia en diagonal y el resto disperso.

    T(m,f) = sum_{c<=m} c chi(c).  La recta m = (f-1)/2 esta CERRADA para todo caracter impar
    primitivo: alli T = (f/2) B_{1,chi} (chi(2)^{-1} - 1), que se anula exactamente cuando
    chi(2) = 1 -- para el cuadratico, f = 7 (mod 8).  Los puntos fuera de ella son la parte
    abierta.  Los rombos vacios marcan las casillas donde ademas muere un caracter de ORDEN > 2:
    no estan en el rango que barrio el paper, y son las que impiden decir "solo el cuadratico".
    """
    fs = [f for f in range(3, 160, 2)]
    diag, triv, otros = [], [], []
    for f in fs:
        if any(f % (d * d) == 0 for d in range(2, int(f ** 0.5) + 1)):
            continue
        for m in range(1, f):
            s = sum(c * kron(c, f) for c in range(1, m + 1))
            if s != 0:
                continue
            if m == f - 1 and f % 4 == 1:
                triv.append((f, m))        # caracter PAR: c <-> f-c anula la suma sola
            elif 2 * m + 1 == f:
                diag.append((f, m))        # la familia cerrada por Dirichlet
            else:
                otros.append((f, m))

    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    ax.plot([3, 160], [1, 79.5], "-", lw=0.9, color="#C3CBD4", zorder=1)
    ax.text(152, 88, "$m=(f-1)/2$", fontsize=8, color="#7A8794", ha="right", va="bottom")
    if triv:
        ax.plot([x[0] for x in triv], [x[1] for x in triv], "s", ms=4.5, color=GRIS,
                mec="white", mew=0.5, zorder=2,
                label=T("$m=f-1$, $f\\equiv1\\ (4)$: even character, trivial (%d)",
                         "$m=f-1$, $f\\equiv1\\ (4)$: car\u00e1cter par, trivial (%d)") % len(triv))
    if otros:
        ax.plot([x[0] for x in otros], [x[1] for x in otros], "o", ms=5, color=NARANJA,
                mec="white", mew=0.5, zorder=3,
                label=T("quadratic $\\chi$, the rest: congruence families, open (%d)",
                         "$\\chi$ cuadr\u00e1tico, el resto: familias de congruencias, abierto (%d)")
                       % len(otros))
    if diag:
        ax.plot([x[0] for x in diag], [x[1] for x in diag], "D", ms=5.5, color=AZUL,
                mec="white", mew=0.5, zorder=4,
                label=T("$m=(f-1)/2$: closed by the criterion $\\chi(2)=1$ (%d)",
                         "$m=(f-1)/2$: cerrado por el criterio $\\chi(2)=1$ (%d)") % len(diag))
    mayores = orden_mayor()
    if mayores:
        # el mismo (f,m) puede aparecer con dos ORDENES distintos (f=127 con 6 y con 18): las
        # casillas son menos que las entradas, y la leyenda cuenta CASILLAS, que es lo que se ve
        celdas = sorted(set((x[0], x[1]) for x in mayores))
        ax.plot([x[0] for x in celdas], [x[1] for x in celdas], "D", ms=13, mfc="none",
                mec=ROJO, mew=1.3, zorder=5,
                label=T("there an odd $\\chi$ of order $>2$ vanishes too (%d cells, $f$ prime)",
                         "all\u00ed se anula tambi\u00e9n un $\\chi$ impar de orden $>2$ (%d casillas, $f$ primo)")
                      % len(celdas))
        # NADA de etiqueta dentro: en esta zona el texto se come la diagonal y dos rombos.  Cual es
        # la casilla de la izquierda y que defecto tiene lo dice el PIE, que es donde se mantiene.
        print("              | orden > 2 en %d casillas: %s" % (len(celdas), celdas))
    ax.set_xlabel(T("conductor $f$ of the character", "conductor $f$ del car\u00e1cter"))
    ax.set_ylabel("$m$")
    ax.set_title(T("where an odd component vanishes: $\\sum_{c\\leq m}c\\,\\chi(c)=0$",
                 "d\u00f3nde se anula una componente impar: $\\sum_{c\\leq m}c\\,\\chi(c)=0$"),
                 fontsize=10.5)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.94)
    ax.grid(alpha=0.18, lw=0.5)
    ax.set_xlim(0, 162)
    fig.tight_layout()
    guardar(fig, "fig_ceros")
    plt.close(fig)
    print("fig_ceros.pdf | ceros: %d en la diagonal, %d fuera" % (len(diag), len(otros)))




def fig_historia():
    """TRES carriles, y la nota esta en el cruce.

    Hasta la compuerta del complemento la figura tenia dos: la pregunta del VALOR, contestada tres
    veces en cincuenta anos, y la del ANILLO, planteada una vez para grupos finitos y una para un
    grupo compacto.  Faltaba el tercero, y es el que enlaza: la pregunta del EMPAREJAMIENTO --
    cuando dos grupos distintos dan los mismos valores -- contestada desde 1982 en otra comunidad,
    la de los grupos de lazos y las algebras de fusion, sin tocar nunca la palabra "anillo".

    La nota vive en el cruce: nuestra operacion (b) ES ese emparejamiento, redemostrada con el
    lema de division y dicha como igualdad de ordenes.  Por eso el carril de abajo no es una
    anomalia aislada: es la esquina vacia de un cuadro que ya tenia dos lados.
    """
    hitos_valor = [
        (1976, T("Kostant\nCoxeter class:\n$0,\\pm1$",
                  "Kostant\nclase de Coxeter:\n$0,\\pm1$")),
        # el hito del medio son DOS trabajos y el criterio de anulacion es de KOSTANT (2004);
        # atribuirselo a Cellini-Moseneder Frajria-Papi, que dan la descripcion por grupos de Weyl
        # afines en los tipos clasicos (2006), era un error nuestro.  Punto entre las dos fechas.
        (2006, T("Kostant 2004; Cellini\u2013M\u00f6seneder\nFrajria\u2013Papi 2006\nwhich weights survive",
                  "Kostant 2004; Cellini\u2013M\u00f6seneder\nFrajria\u2013Papi 2006\nqu\u00e9 pesos sobreviven")),
        (2025, T("Nadimpalli\u2013Pattanayak\nPrasad\nvalues as dimensions",
                  "Nadimpalli\u2013Pattanayak\nPrasad\nvalores como dimensiones")),
    ]
    # El carril del emparejamiento.  Frenkel 1982 es el origen; Pressley-Segal 1986 lo tiene como
    # simetria k <-> N-k de representaciones de grupos de lazos; Nakanishi-Tsuchiya 1992 lo prueba
    # para modelos WZW; Witten 1993 lo deriva de G(k,N) = G(N-k,N).
    hitos_par = [
        (1982, T("Frenkel\nrank\u2013level", "Frenkel\nrango\u2013nivel")),
        (1986, "Pressley\u2013Segal\n$k\\leftrightarrow N{-}k$"),
        (1992, "Nakanishi\u2013Tsuchiya\nWZW"),
        (1993, "Witten\n$G(k,N)\\cong G(N{-}k,N)$"),
    ]
    hitos_anillo = [
        (2019, T("B\u00e4chle\u2013Sambale\norders, finite groups",
                  "B\u00e4chle\u2013Sambale\n\u00f3rdenes, grupos finitos")),
        (2026, T("the conductor for $\\mathrm{Sp}(2m)$\nand this note",
                  "el conductor de $\\mathrm{Sp}(2m)$\ny esta nota")),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 4.15))
    VERDE = "#4E8C4A"
    for y, nombre, color in ((2.0, T("the VALUE", "el VALOR"), AZUL),
                             (1.0, T("the PAIRING", "el EMPAREJAMIENTO"), VERDE),
                             (0.0, T("the RING", "el ANILLO"), NARANJA)):
        ax.axhline(y, color="#C3CBD4", lw=1.0, zorder=1)
        ax.text(1970.5, y, nombre, fontsize=9, color=color, va="center", ha="right", weight="bold")
    for (x, txt) in hitos_valor:
        ax.plot([x], [2.0], "o", ms=8, color=AZUL, mec="white", mew=1.0, zorder=3)
        ax.annotate(txt, (x, 2.0), textcoords="offset points", xytext=(0, 12), ha="center",
                    fontsize=7.0, color="#33404D")
    # las cuatro del carril del medio caen en once anos: se escalonan arriba/abajo para no pisarse
    for i, (x, txt) in enumerate(hitos_par):
        ax.plot([x], [1.0], "o", ms=7, color=VERDE, mec="white", mew=1.0, zorder=3)
        arriba = (i % 2 == 0)
        ax.annotate(txt, (x, 1.0), textcoords="offset points",
                    xytext=(0, 11 if arriba else -26), ha="center", fontsize=6.8,
                    color="#33404D")
    for i, (x, txt) in enumerate(hitos_anillo):
        ax.plot([x], [0.0], "o", ms=8, color=NARANJA, mec="white", mew=1.0, zorder=3)
        ax.annotate(txt, (x, 0.0), textcoords="offset points", xytext=(0, -30 if i == 0 else -56),
                    ha="center", fontsize=7.0, color="#33404D")
    # EL CRUCE: la nota toca los tres carriles.  Se dibuja como una llave vertical en 2026, no
    # como una flecha, porque no es una implicacion: es el mismo trabajo visto desde tres sitios.
    ax.annotate("", xy=(2026, 2.0), xytext=(2026, 0.0),
                arrowprops=dict(arrowstyle="-", lw=1.2, color=ROJO, ls=(0, (4, 3))), zorder=2)
    ax.plot([2026], [1.0], "o", ms=6, color=ROJO, mec="white", mew=1.0, zorder=4)
    # El rotulo NO va centrado en y=1.0: ahi la linea horizontal del carril le cruza el renglon
    # del medio y parece tachado.  Va apoyado encima de la linea, y empieza en 2028.8 para que
    # tampoco lo toque la vertical de puntos.
    ax.text(2028.8, 1.10, T("this note:\nthe same ring,\nseen from all three",
                            "esta nota:\nel mismo anillo,\nvisto desde los tres"),
            fontsize=7.2, color=ROJO, va="bottom", ha="left")
    ax.set_xlim(1967, 2041)
    ax.set_ylim(-1.35, 2.72)
    ax.set_yticks([])
    ax.set_xticks([1976, 1982, 1993, 2006, 2019, 2026])
    ax.set_xticklabels(["1976", "1982", "1993", "2006", "2019", "2026"], fontsize=8)
    for lado in ("left", "right", "top"):
        ax.spines[lado].set_visible(False)
    ax.spines["bottom"].set_color("#C3CBD4")
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    guardar(fig, "fig_historia")
    plt.close(fig)
    print("fig_historia.pdf | %d hitos en el valor, %d en el emparejamiento, %d en el anillo" % (
        len(hitos_valor), len(hitos_par), len(hitos_anillo)))


def fig_puente():
    """LA SOMBRA ES EL ANILLO REAL MAS EL MISMO, DESPLAZADO.

    El enunciado central de la seccion de la sombra es una convolucion:
    W_i(B) = W_i(R) + W_{i-s}(R).  Es una frase sobre sumas de dimensiones y se lee mucho mejor
    apilada: la barra clara es W(R), la oscura es la MISMA barra corrida s posiciones, y la altura
    total es W(B).  La linea de puntos es la dimension llena 2d: donde la suma la alcanza, el
    conductor para -- y por eso e_B = e_R + s se VE, no hay que creerselo.

    Tres paneles, y los tres son casillas medidas, no dibujos:
      (a) delta = 0, s = 1   -- el caso generico,      (r,p) = (11,3)
      (b) delta = 1, s = 1   -- el defecto de h^-,     (r,p) = (23,3), q = 207
      (c) s = 3              -- p | N, el desplazamiento general, (r,p,N) = (5,3,3), q = 45
    Los perfiles de (a) y (b) salen de la ley; el de (c) es la medida de puente_general_s.sage.
    """
    casos = [
        (T("(a) $\\delta=0$, shift $1$", "(a) $\\delta=0$, desplazamiento $1$"), [1, 5, 5, 5, 5], 1, 5,
         "$r=11$, $p=3$"),
        (T("(b) $\\delta=1$, shift $1$", "(b) $\\delta=1$, desplazamiento $1$"), [1, 10, 11, 11, 11], 1, 11,
         "$r=23$, $p=3$, $q=207$"),
        (T("(c) $p\\mid N$, shift $3$", "(c) $p\\mid N$, desplazamiento $3$"), [1, 0, 1, 2, 2, 2], 3, 2,
         "$r=5$, $p=3$, $N=3$, $q=45$"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.05))
    for eje, (titulo, WR, s, d, pie) in zip(axes, casos):
        L = len(WR)
        xs = np.arange(L)
        base = np.array(WR, dtype=float)
        desp = np.array([WR[i - s] if i - s >= 0 else 0 for i in range(L)], dtype=float)
        eje.bar(xs, base, width=0.66, color="#BFD3E6", edgecolor=AZUL, linewidth=0.8,
                label="$W_i(R)$", zorder=3)
        eje.bar(xs, desp, width=0.66, bottom=base, color=NARANJA, edgecolor="white",
                linewidth=0.8, alpha=0.9, label="$W_{i-s}(R)$", zorder=3)
        eje.axhline(2 * d, ls=":", lw=1.0, color=GRIS, zorder=2)
        eje.text(L - 0.45, 2 * d, "$2d$", fontsize=7.5, color=GRIS, va="bottom", ha="right")
        # El conductor: primer grado a partir del cual la suma llena.  OJO: si la ventana medida
        # se acaba antes de llenar, el "primer grado" seria el borde de la ventana y NO el
        # conductor.  Eso pasa en el panel (c), donde K = 6 y el perfil no llega a 2d: ahi se
        # dibuja una cota, no un valor.  Un capped sweep tiene que decir que esta capped.
        tot = base + desp
        lleno = (tot[L - 1] == 2 * d)
        if lleno:
            e = next(i for i in range(L) if all(tot[j] == 2 * d for j in range(i, L)))
            et, x = "$e_B=%d$" % e, e - 0.5
        else:
            et, x = "$e_B>%d$" % (L - 1), L - 0.55
        eje.annotate("", xy=(x, 2 * d * 1.12), xytext=(x, 0),
                     arrowprops=dict(arrowstyle="-", lw=1.1, color=ROJO, ls="--"), zorder=4)
        eje.text(x - (0 if lleno else 1.15) + 0.08, 2 * d * 1.13, et, fontsize=8, color=ROJO,
                 va="bottom")
        eje.set_xticks(xs)
        eje.set_xticklabels([str(i) for i in xs], fontsize=8)
        eje.set_xlabel(T("layer $i$", "capa $i$"), fontsize=8.5)
        eje.set_title(titulo, fontsize=9.5, color=AZUL)
        eje.text(0.5, -0.30, pie, transform=eje.transAxes, ha="center", fontsize=7.5,
                 color="#55606C", style="italic")
        eje.set_ylim(0, 2 * d * 1.30)
        eje.tick_params(axis="y", labelsize=8)
        for lado in ("top", "right"):
            eje.spines[lado].set_visible(False)
        eje.spines["left"].set_color("#C3CBD4")
        eje.spines["bottom"].set_color("#C3CBD4")
    axes[0].legend(loc="upper left", fontsize=7.5, frameon=False, handlelength=1.1)
    fig.subplots_adjust(bottom=0.30, wspace=0.28, top=0.86)
    fig.suptitle(T("the shadow is the real ring plus itself, shifted",
                   "la sombra es el anillo real m\u00e1s \u00e9l mismo, desplazado"),
                 y=1.01, fontsize=10.5)
    guardar(fig, "fig_puente")
    plt.close(fig)
    print("fig_puente.pdf | 3 paneles: shift 1 sin defecto, shift 1 con defecto, shift 3")


def fig_intervalos():
    """POR QUE UN INTERVALO NO PUEDE SER UNION DE ESTRATOS, salvo casi todo o casi nada.

    Es la figura de lem:inter, que es lo que cierra prop:sop -- la unica frontera que bloqueaba
    el envio.  El argumento es geometrico y se ve entero en tres circulos:

      (a) el arco I que contiene las unidades tiene que contener 0, porque contiene 1 y q'-1;
      (b) su complemento K no tiene ninguna unidad, y si c y c+1 estan en K sus gcd D y D' son
          COPRIMOS; cada estrato X_D contiene D y q'-D, CASI ANTIPODALES, asi que K tiene que ser
          largo -- de ahi ((B+1)/2)^2 <= D D' <= q';
      (c) q' = 6 es el unico sitio donde cabe, y ahi phi(q') = 2.

    Los colores son los estratos gcd(c,q'); las unidades van en azul.
    """
    QPS = [30, 30, 6]
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.6))
    PAL = {1: AZUL, 2: "#7FA8C9", 3: "#4E8C4A", 5: NARANJA, 6: "#B08A3E",
           10: "#8E6FB0", 15: "#C08A8A", 30: GRIS}
    for idx, (eje, qp) in enumerate(zip(axes, QPS)):
        eje.set_aspect("equal")
        eje.axis("off")
        ang = lambda c: np.pi / 2 - 2 * np.pi * c / qp
        th = np.linspace(0, 2 * np.pi, 300)
        eje.plot(np.cos(th), np.sin(th), "-", lw=0.7, color="#DCE2E8", zorder=0)
        # el arco destacado
        if idx == 0:
            arco = [c % qp for c in range(-7, 8)]        # I, alrededor de 0
            et = T("$I\\ni 0,\\pm1$: must contain every unit",
                   "$I\\ni 0,\\pm1$: tiene que contener toda unidad")
        elif idx == 1:
            arco = [c % qp for c in range(8, 23)]        # K, el complementario
            # OJO al leer el panel: los puntos marcados (D y q'-D) caen FUERA del arco rojo, y eso
            # ES la contradiccion -- X_D tendria que estar entero dentro de K.  El rotulo lo dice.
            et = T("$c,c{+}1\\in K$ give coprime $D,D'$, but $X_D$ spills out of $K$",
                   "$c,c{+}1\\in K$ dan $D,D'$ coprimos, pero $X_D$ se sale de $K$")
        else:
            arco = [2, 3, 4]
            et = T("$q'=6$, $B=3$: the only one", "$q'=6$, $B=3$: el \u00fanico")
        A = set(arco)
        for c in range(qp):
            D = math.gcd(c, qp)
            col = PAL.get(D, GRIS)
            dentro = c in A
            eje.plot(np.cos(ang(c)), np.sin(ang(c)), "o", ms=9 if dentro else 5.5,
                     color=col if dentro else "white", mec=col,
                     mew=1.4 if dentro else 1.0, zorder=3)
        # el arco, dibujado como sector
        if arco:
            a0, a1 = ang(arco[0]), ang(arco[-1])
            tt = np.linspace(a0, a1 if a1 < a0 else a1 - 2 * np.pi, 120)
            eje.plot(1.13 * np.cos(tt), 1.13 * np.sin(tt), "-", lw=2.4,
                     color=ROJO if idx == 1 else AZUL, alpha=0.75, zorder=2)
        # en el panel (b), marcar D y q'-D casi antipodales
        if idx == 1:
            for D, col in ((2, "#7FA8C9"), (3, "#4E8C4A")):
                for cc in (D, qp - D):
                    eje.plot(np.cos(ang(cc)), np.sin(ang(cc)), "o", ms=13, mfc="none",
                             mec=col, mew=1.6, zorder=4)
                eje.annotate("", xy=(np.cos(ang(qp - D)), np.sin(ang(qp - D))),
                             xytext=(np.cos(ang(D)), np.sin(ang(D))),
                             arrowprops=dict(arrowstyle="<->", lw=1.0, color=col,
                                             connectionstyle="arc3,rad=0.18"), zorder=2)
            eje.text(0, -0.08, T("$D$ and $q'-D$ nearly antipodal,\nand both outside $K$",
                                 "$D$ y $q'-D$ casi ant\u00edpodas,\ny los dos fuera de $K$"), ha="center",
                     fontsize=7.5, color="#55606C", style="italic")
        eje.set_xlim(-1.45, 1.45)
        eje.set_ylim(-1.62, 1.45)
        eje.text(0, 1.30, "(%s) $q'=%d$" % ("abc"[idx], qp), ha="center", fontsize=9.5,
                 color=AZUL)
        eje.text(0, -1.46, et, ha="center", fontsize=7.8, color="#55606C")
    fig.suptitle(T("an interval that is a union of $\\gcd$ strata is almost all or almost none",
                   "un intervalo que es uni\u00f3n de estratos $\\gcd$ es casi todo o casi nada"),
                 y=1.00, fontsize=10.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.06, wspace=0.05)
    guardar(fig, "fig_intervalos")
    plt.close(fig)
    print("fig_intervalos.pdf | 3 circulos: el arco de las unidades, el complementario, y q'=6")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] not in ("en", "es"):
            sys.exit("uso: python figuras.py [en|es]")
        IDIOMA = sys.argv[1]
    print("=== figuras en '%s'" % IDIOMA)
    fig_plano()
    fig_operaciones()
    fig_formas()
    fig_ceros()
    fig_historia()
    fig_puente()
    # fig_intervalos() ilustraba la prueba de seis casos de lem:inter; el 19-sep-2026 la prueba
    # pasa a ser por autocorrelacion y la figura sale del texto.  Se conserva la funcion.
