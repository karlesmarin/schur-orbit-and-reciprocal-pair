# -*- coding: utf-8 -*-
"""Figuras de la nota II.  Los datos salen de las salidas archivadas en esta misma carpeta, no se escriben
a mano.

fig_mezcla.pdf -- EL RECUENTO a(m) CONTRA LAS DOS PREDICCIONES.  Cociente a(m)/prediccion para la
gaussiana global (media y varianza de todas las elecciones de signos) y para la mezcla (gaussiana
condicionada a los signos de los primos <= sqrt(m)).  La mezcla se asienta en 1; la gaussiana global
se queda corta por un factor casi constante.  Fuente: II_mezcla_OUT.txt.

fig_margen.pdf -- POR QUE LA LINEA A CIERRA.  (a) K(m)/m, el ancho que las sumas con signo de los
primos de (m/2, m] dejan sin cubrir en cada extremo, medido exacto: tiende a ~4.  (b) El resto |R|
tras fijar los signos inferiores, frente a la ventana P - 4m que los primos superiores cubren,
medido hasta m = 10^6 (gris); la prueba: B*/D en cada uno de los 104 bloques certificados de
[5000, 10^8) (azul) y la cota analitica (B+K)/P_Q desde 10^8 (naranja).  Todo por debajo de 1.
Fuentes: II_lema_signos_OUT.txt, II_ruta_A_OUT.txt, II_certificados_A_bloques.csv,
II_cola_analitica.py.

Authors: Carles Marin, Claude (AI assistant).
Run:  python figuras_II.py [en|es]      (es -> fig_x_es.pdf, para segmento_II_es.tex)
"""
import math
import pathlib
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = pathlib.Path(__file__).resolve().parent
G = AQUI  # en el paquete publicado, las salidas viven junto a este script
AZUL = "#1F4E79"
NARANJA = "#C2610A"
GRIS = "#7A7A7A"

# "en" -> fig_x.pdf ;  "es" -> fig_x_es.pdf.  Mismo patron que note_segmento/figuras.py.
# Las dos ediciones salen del MISMO codigo con los MISMOS datos: solo cambian las cadenas.
IDIOMA = "en"


def T(en, es):
    return es if IDIOMA == "es" else en


def guardar(fig, nombre):
    """Anade el sufijo del idioma al nombre, y nada mas."""
    out = AQUI / (nombre + ("" if IDIOMA == "en" else "_" + IDIOMA) + ".pdf")
    fig.savefig(out, bbox_inches="tight")
    return out.name


def fig_mezcla():
    ms, rg, rm = [], [], []
    for ln in (G / "II_mezcla_OUT.txt").read_text(encoding="utf-8", errors="replace").splitlines():
        c = ln.split()
        if len(c) == 4 and c[0].isdigit():
            m = int(c[0])
            if m >= 20:
                ms.append(m); rg.append(float(c[2])); rm.append(float(c[3]))
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.axhline(1, color=GRIS, lw=0.8, ls="--")
    ax.plot(ms, rg, "o", ms=3.5, color=NARANJA, label=T("global Gaussian", "gaussiana global"))
    ax.plot(ms, rm, "o", ms=3.5, color=AZUL,
            label=T("mixture (conditioned on the primes $\\leq\\sqrt{m}$)",
                    "mezcla (condicionada a los primos $\\leq\\sqrt{m}$)"))
    ax.set_xlabel(T("$m$ (admissible rows, $m\\equiv0,3$ mod $4$)",
                    "$m$ (filas admisibles, $m\\equiv0,3$ mod $4$)"))
    ax.set_ylabel(T("$a(m)$ / prediction", "$a(m)$ / predicción"))
    ax.set_ylim(0.6, 3.0)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(T("The count of solutions against two predictions",
                   "El número de soluciones frente a dos predicciones"), fontsize=10)
    guardar(fig, "fig_mezcla")
    return len(ms)


def fig_margen():
    km, kk = [], []
    for ln in (G / "II_lema_signos_OUT.txt").read_text(encoding="utf-8", errors="replace").splitlines():
        mo = re.match(r"m=\s*(\d+) .* K=\s*(\d+)", ln)
        if mo and int(mo.group(1)) >= 100:
            km.append(int(mo.group(1))); kk.append(int(mo.group(2)) / int(mo.group(1)))
    rm, rr = [], []
    for ln in (G / "II_ruta_A_OUT.txt").read_text(encoding="utf-8", errors="replace").splitlines():
        mo = re.match(r"m=\s*(\d+) .*\|R\|/\(P-4m\)=([0-9.e+-]+)", ln)
        if mo:
            rm.append(int(mo.group(1))); rr.append(float(mo.group(2)))
    # (b) la prueba: bloques certificados (B*/D por bloque, II_certificados_A_bloques.csv) y,
    # desde 10^8, la cota analitica (B+K)/P_Q <= (4/3) L sqrt(L+13)/sqrt(M) + 3808 L^2/M
    # (II_cola_analitica.py)
    bx, by = [], []
    for ln in (G / "II_certificados_A_bloques.csv").read_text(encoding="utf-8").splitlines()[1:]:
        c = ln.split(",")
        bx += [int(c[0]), int(c[1])]
        r = math.sqrt(float(c[6]))
        by += [r, r]
    xs = [10 ** (8 + i * 0.05) for i in range(0, 41)]
    ys = [(4 / 3) * math.log(x) * math.sqrt(math.log(x) + 13) / math.sqrt(x)
          + 3808 * math.log(x) ** 2 / x for x in xs]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.6))
    a1.plot(km, kk, "o-", ms=3.5, color=AZUL)
    a1.axhline(4, color=GRIS, lw=0.8, ls="--")
    a1.set_xlabel("$m$"); a1.set_ylabel("$K(m)/m$")
    a1.set_title(T("(a) uncovered width of the signed sums",
                   "(a) anchura sin cubrir de las sumas con signo"), fontsize=10)
    a2.loglog(rm, rr, "o", ms=4, color=GRIS, label=T("measured $|R|/(P-4m)$", "$|R|/(P-4m)$ medido"))
    a2.loglog(bx, by, "-", color=AZUL, lw=1.2,
              label=T("certified blocks: $B_*/D$", "bloques certificados: $B_*/D$"))
    a2.loglog(xs, ys, "-", color=NARANJA, lw=1.4,
              label=T("analytic bound $(B+K)/P_Q$", "cota analítica $(B+K)/P_Q$"))
    a2.axvline(1e8, color=NARANJA, lw=0.8, ls="--")
    a2.axhline(1, color=GRIS, lw=0.8, ls="--")
    a2.set_xlabel("$m$"); a2.set_ylabel(T("remainder / window", "resto / ventana"))
    a2.set_title(T("(b) the margin of the proof", "(b) el margen de la demostración"), fontsize=10)
    a2.legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.tight_layout()
    guardar(fig, "fig_margen")
    return len(km), len(rm)


def fig_historia():
    """TRES carriles y la nota en el cruce, como en la nota I.  Fechas verificadas contra la fuente
    (informe de verificacion del 19-sep): Jacobi 1832 enuncia la semisuma, Dirichlet 1838 la
    demuestra y se la atribuye; Maillet 1913 pregunta, Carlitz-Olson 1955 evaluan (y dicen que
    Chowla y Weil lo tenian); Polya 1919 OBSERVA, Haselgrove 1958 refuta; Erdos hacia 1932, impreso
    1957, Tao 2016; Sarkozy 1994 es 'Finite addition theorems II'."""
    ROJO = "#B03030"
    VERDE = "#4E8C4A"
    semis = [(1832, T("Jacobi\nhalf sums, $h(-p)$", "Jacobi\nsemisumas, $h(-p)$")),
             (1838, T("Dirichlet\nproof", "Dirichlet\ndemostración")),
             (1913, T("Maillet\ndeterminant", "Maillet\ndeterminante")),
             (1955, "Carlitz–Olson\n$p^{(p-3)/2}h^-$"),
             (1975, T("Jager–Lenstra\ncosecants", "Jager–Lenstra\ncosecantes")), (1987, "Girstmair"),
             (2001, "Kučera,\nKanemitsu–Kuzumaki")]
    signos = [(1919, "Pólya\n$\\sum\\lambda(n)\\leq0$?"),
              (1932, T("Erdős\ndiscrepancy", "Erdős\ndiscrepancia")),
              (1958, T("Haselgrove\ndisproof", "Haselgrove\nrefutación")), (2010, "Borwein–Choi–Coons"),
              (2016, "Tao\nEDP"), (2024, T("Conrey\nRH", "Conrey\nHR")),
              (2026, "Klurman–Munsch–Sun")]
    sumas = [(1994, T("Sárközy\nprogressions in\nsubset sums", "Sárközy\nprogresiones en\nsumas de subconjuntos")),
             (2003, T("Lev\nconstants", "Lev\nconstantes")),
             (2006, "Szemerédi–Vu"), (2021, "Bringmann–Wellnitz;\nConlon–Fox–Pham")]
    fig, ax = plt.subplots(figsize=(9.6, 4.4))
    for y, nombre, color in ((2.0, T("HALF SUMS", "SEMISUMAS"), AZUL), (1.0, T("SIGNS", "SIGNOS"), VERDE),
                             (0.0, T("SUBSET SUMS", "SUMAS DE\nSUBCONJUNTOS"), NARANJA)):
        ax.axhline(y, color="#C3CBD4", lw=1.0, zorder=1)
        ax.text(1822, y, nombre, fontsize=9, color=color, va="center", ha="right", weight="bold")
    # posiciones a mano donde los hitos caen juntos (2010/2016/2026 y 1994/2003/2006)
    a_mano = {("SIGNS", 2010): (-6, -30, "right"), ("SIGNS", 2016): (0, -30, "center"),
              ("SIGNS", 2024): (-3, 10, "center"), ("SIGNS", 2026): (6, -30, "left"),
              ("SUBSET", 1994): (-4, 10, "right"), ("SUBSET", 2003): (0, -30, "center"),
              ("SUBSET", 2006): (0, 30, "center"), ("SUBSET", 2021): (4, -48, "center")}
    for lista, y, color, clave in ((semis, 2.0, AZUL, "HALF"), (signos, 1.0, VERDE, "SIGNS"),
                                   (sumas, 0.0, NARANJA, "SUBSET")):
        for i, (x, txt) in enumerate(lista):
            ax.plot([x], [y], "o", ms=7, color=color, mec="white", mew=1.0, zorder=3)
            dx, dy, ha = a_mano.get((clave, x), (0, 10 if i % 2 == 0 else -30, "center"))
            ax.annotate(txt, (x, y), textcoords="offset points", xytext=(dx, dy),
                        ha=ha, fontsize=6.6, color="#33404D")
    ax.annotate("", xy=(2026, 2.0), xytext=(2026, 0.0),
                arrowprops=dict(arrowstyle="-", lw=1.2, color=ROJO, ls=(0, (4, 3))), zorder=2)
    ax.plot([2026], [1.0], "o", ms=6, color=ROJO, mec="white", mew=1.0, zorder=4)
    ax.text(2032, 1.10, T("this note:\nexact zeros of\nweighted sums",
                          "esta nota:\nceros exactos de\nsumas ponderadas"), fontsize=7.2, color=ROJO,
            va="bottom", ha="left")
    ax.set_xlim(1815, 2050)
    ax.set_ylim(-0.9, 2.75)
    ax.set_yticks([])
    ax.set_xticks([1832, 1913, 1958, 1994, 2026])
    ax.set_xticklabels(["1832", "1913", "1958", "1994", "2026"], fontsize=8)
    for lado in ("left", "right", "top"):
        ax.spines[lado].set_visible(False)
    ax.spines["bottom"].set_color("#C3CBD4")
    fig.tight_layout()
    guardar(fig, "fig_historia_II")
    return len(semis) + len(signos) + len(sumas)


def _jacobi(a, n):
    """simbolo de Jacobi (a/n), n impar positivo."""
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


def fig_ceros():
    """Donde se anula T(m, chi) = sum_{c <= m} c chi(c), chi cuadratico de conductor f impar libre de
    cuadrados, f < 160 (antes en la nota I; desde el 19-sep vive aqui).  La diagonal m = (f-1)/2 la
    cierra el teorema de la diagonal; la familia gris m = f-1 es la de los caracteres pares; el resto
    son familias de congruencias por fila, y toda fila admisible esta ocupada (Corolario de las filas).
    Los rombos rojos: casillas donde se anula tambien un caracter impar de orden > 2
    (caracteres_orden_mayor_OUT.txt)."""
    import ast
    ROJO = "#B03030"
    diag, triv, otros = [], [], []
    for f in range(3, 160, 2):
        if any(f % (d * d) == 0 for d in range(2, int(f ** 0.5) + 1)):
            continue
        for m in range(1, f):
            if sum(c * _jacobi(c, f) for c in range(1, m + 1)) != 0:
                continue
            if m == f - 1 and f % 4 == 1:
                triv.append((f, m))
            elif 2 * m + 1 == f:
                diag.append((f, m))
            else:
                otros.append((f, m))
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    ax.plot([3, 160], [1, 79.5], "-", lw=0.9, color="#C3CBD4", zorder=1)
    ax.text(152, 88, "$m=(f-1)/2$", fontsize=8, color="#7A8794", ha="right", va="bottom")
    ax.plot([x[0] for x in triv], [x[1] for x in triv], "s", ms=4.5, color=GRIS, mec="white", mew=0.5,
            zorder=2, label=T("$m=f-1$, $f\\equiv1\\ (4)$: even character, trivial (%d)",
                              "$m=f-1$, $f\\equiv1\\ (4)$: carácter par, trivial (%d)") % len(triv))
    ax.plot([x[0] for x in otros], [x[1] for x in otros], "o", ms=5, color=NARANJA, mec="white",
            mew=0.5, zorder=3, label=T("quadratic $\\chi$, the rest: row congruence families (%d)",
                                       "$\\chi$ cuadrático, el resto: familias de congruencias por fila (%d)")
            % len(otros))
    ax.plot([x[0] for x in diag], [x[1] for x in diag], "D", ms=5.5, color=AZUL, mec="white", mew=0.5,
            zorder=4, label=T("$m=(f-1)/2$: closed by the criterion $\\chi(2)=1$ (%d)",
                              "$m=(f-1)/2$: cerrada por el criterio $\\chi(2)=1$ (%d)") % len(diag))
    fuente = G / "capas_tipo" / "caracteres_orden_mayor_OUT.txt"
    if not fuente.exists():
        fuente = G / "caracteres_orden_mayor_OUT.txt"
    marca = "casillas (f, m, orden) con orden > 2:"
    mayores = []
    for ln in fuente.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith(marca):
            mayores = ast.literal_eval(ln[len(marca):].strip())
    celdas = sorted(set((x[0], x[1]) for x in mayores))
    ax.plot([x[0] for x in celdas], [x[1] for x in celdas], "D", ms=13, mfc="none", mec=ROJO, mew=1.3,
            zorder=5, label=T("there an odd $\\chi$ of order $>2$ vanishes too (%d cells, $f$ prime)",
                              "allí se anula también un $\\chi$ impar de orden $>2$ (%d casillas, $f$ primo)")
            % len(celdas))
    ax.set_xlabel(T("conductor $f$ of the character", "conductor $f$ del carácter"))
    ax.set_ylabel("$m$")
    ax.set_title(T("where $T(m,\\chi)=\\sum_{c\\leq m}c\\,\\chi(c)$ vanishes",
                   "dónde se anula $T(m,\\chi)=\\sum_{c\\leq m}c\\,\\chi(c)$"), fontsize=10.5)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.94)
    ax.grid(alpha=0.18, lw=0.5)
    ax.set_xlim(0, 162)
    fig.tight_layout()
    guardar(fig, "fig_ceros")
    plt.close(fig)
    return len(diag), len(triv), len(otros), len(celdas)


def _tablas(q):
    L = [0] + [-1] * (q - 1)
    for r in range(1, q):
        L[r * r % q] = 1
    C = [0] * q; T = [0] * q; c = t = 0
    for r in range(q):
        c += L[r]; t += r * L[r]; C[r] = c; T[r] = t
    return L, C, T


def _S_rec(q, s, m, C, T):
    if m == 0:
        return 0
    k, r = divmod(m, q)
    return s * q * _S_rec(q, s, k, C, T) + k * (T[q - 1] + q * C[r]) + T[r]


def fig_caracteres():
    """(a) q = 103: F_q(y) = T(r) + y (h - C(r)) en [0, q], con las celdas nulas (R_103 = {47, 51, 55});
    (b) q = 59: F_q en [0, q/2], sin celda nula;  (c) ceros de S_chi (s = +1) en escala logaritmica:
    el arbol de digitos de q = 5, la cadena de q = 11 (h = 1) y los bucles de q = 59 y 131.
    Todo se calcula aqui con la recursion por digitos (Teorema de la nota), en enteros."""
    ROJO = "#B03030"; VERDE = "#4E8C4A"
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(11.2, 3.4), gridspec_kw={"width_ratios": [1.2, 1, 1.2]})
    for ax, q, titulo in ((a1, 103, T("(a) $q=103$: null cells", "(a) $q=103$: celdas nulas")),
                          (a2, 59, T("(b) $q=59$: none", "(b) $q=59$: ninguna"))):
        L, C, Tr = _tablas(q)
        h = -Tr[q - 1] // q
        ys, Fs = [], []
        tope = q if q == 103 else (q + 1) // 2
        for r in range(tope):
            for y in (r, r + 1):
                if q != 103 and y > q / 2:
                    y = q / 2          # (b): hasta x = 1/2, donde la serie se anula por imparidad
                ys.append(y); Fs.append((Tr[r] + y * (h - C[r])) / q)
        ax.plot(ys, Fs, "-", color=AZUL, lw=1.0)
        ax.axhline(0, color=GRIS, lw=0.7, ls="--")
        R = [r for r in range(q) if Tr[q - 1] + q * C[r] == 0 and Tr[r] == 0]
        for r in R:
            ax.axvspan(r, r + 1, color=ROJO, alpha=0.35, lw=0)
            ax.plot([r, r + 1], [0, 0], "-", color=ROJO, lw=2.2)
        ax.set_xlabel("$y=qx$"); ax.set_title(titulo, fontsize=10)
        if R:
            # recuadro ampliado alrededor de las celdas nulas
            ins = ax.inset_axes([0.60, 0.58, 0.37, 0.37])
            ins.plot(ys, Fs, "-", color=AZUL, lw=1.0)
            ins.axhline(0, color=GRIS, lw=0.6, ls="--")
            for r in R:
                ins.axvspan(r, r + 1, color=ROJO, alpha=0.35, lw=0)
                ins.plot([r, r + 1], [0, 0], "-", color=ROJO, lw=2.2)
            ins.set_xlim(44, 59); ins.set_ylim(-0.06, 0.06)
            ins.set_xticks(R); ins.tick_params(labelsize=6.5); ins.set_yticks([])
            ins.set_title(T("cells %s", "celdas %s") % ", ".join(map(str, R)), fontsize=7.5, color=ROJO)
    a1.set_ylabel(T("$F_q(y)/q$  (Conrey's $f_q$, rescaled)", "$F_q(y)/q$  ($f_q$ de Conrey, reescalada)"))
    # (c) ceros
    series = []
    L, C, Tr = _tablas(5)
    arbol = sorted({int("".join(d), 5) for j in range(1, 9) for d in __import__("itertools").product("04", repeat=j)} - {0})
    series.append((T("$q=5$: digit tree", "$q=5$: árbol de dígitos"), arbol, AZUL))
    cad11 = [11]
    while len(cad11) < 7:
        cad11.append(11 * (cad11[-1] + 1))
    series.append((T("$q=11$: class number one", "$q=11$: número de clases uno"), cad11, VERDE))
    b59 = 283200; Q = 59 ** 4
    series.append((T("$q=59$: a block", "$q=59$: un bloque"),
                   [59 * b59 * sum(Q ** i for i in range(t)) for t in range(1, 4)], NARANJA))
    series.append((T("$q=131$: a loop", "$q=131$: un bucle"),
                   [27 + 43 * sum(131 ** j for j in range(1, t + 1)) for t in range(0, 6)], ROJO))
    for i, (nombre, zs, color) in enumerate(series):
        a3.plot([math.log10(z) for z in zs], [i] * len(zs), "o", ms=4, color=color)
    a3.set_yticks(range(len(series))); a3.set_yticklabels([s[0] for s in series], fontsize=8)
    a3.set_xlabel(T("$\\log_{10} m$ of the zeros of $S_\\chi$", "$\\log_{10} m$ de los ceros de $S_\\chi$"))
    a3.set_title(T("(c) infinite families", "(c) familias infinitas"), fontsize=10)
    a3.set_ylim(-0.6, len(series) - 0.4)
    # control: todos son ceros (por la recursion)
    for q, zs in ((5, arbol), (11, cad11), (59, series[2][1]), (131, series[3][1])):
        Lq, Cq, Tq = _tablas(q)
        assert all(_S_rec(q, 1, z, Cq, Tq) == 0 for z in zs), q
    fig.tight_layout()
    guardar(fig, "fig_caracteres")
    return sum(len(s[1]) for s in series)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] not in ("en", "es"):
            sys.exit("uso: python figuras_II.py [en|es]")
        IDIOMA = sys.argv[1]
    print("=== figuras en '%s'" % IDIOMA)
    print("fig_ceros.pdf | (diagonal, pares, resto, casillas de orden > 2):", fig_ceros())
    print("fig_caracteres.pdf | ceros dibujados:", fig_caracteres())
    print("fig_historia_II.pdf | hitos:", fig_historia())
    print("fig_mezcla.pdf | filas:", fig_mezcla())
    print("fig_margen.pdf | puntos K, puntos R:", fig_margen())
