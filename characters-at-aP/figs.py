# -*- coding: utf-8 -*-
# LAS FIGURAS DE LA NOTA 2 PARA PRASAD.   25 de agosto de 2026.
#
# ⚠ ESTE FICHERO GENERA SEIS, PERO LA NOTA SOLO USA CUATRO.  La cabecera decia "Cinco figuras"
#   y ninguno de los tres numeros coincidia con otro; una version anterior anuncio seis.
#
# LAS CUATRO QUE VIAJAN (comprobar con: grep -c includegraphics characters_at_aP.tex):
#   fig_map        el mapa: una raiz (Littlewood 1940), dos ramas, un punto de encuentro (a_P)
#   fig_union3d    por que el contenido esta en las POTENCIAS: 100 puntos +-1 contra 504 hasta 72
#   fig_capacity   el PORQUE: en tipo A las cajas saturan; en tipo C sobra una unidad
#   fig_types      donde R_P(q) y R(q) difieren, por tipo y por q
#
# LAS DOS QUE NO ENTRARON (se conservan porque cuestan poco y pueden volver a hacer falta):
#   fig_folding    la involucion r -> q-r sobre Z/q, con sus dos puntos fijos
#   fig_values     el espectro de chi(a_P^2)
#
# Los datos de fig_values, fig_types y fig_union3d se RECALCULAN aqui, no se copian de nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python figs.py > figs_OUT.txt 2>&1     (desde gates/; escribe en ../note_prasad2/)

import sys
from collections import Counter
from fractions import Fraction
from math import gcd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Arc

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

AZUL, NARANJA, GRIS, VERDE = "#1f4e79", "#c55a11", "#7f7f7f", "#2e6f40"
plt.rcParams.update({"font.size": 9, "font.family": "serif", "axes.linewidth": 0.6})


# ---------------------------------------------------------------- 1. EL MAPA
def fig_map():
    # ⚠ Las posiciones van por INDICE, no por anyo.  Con el eje temporal lineal las cinco
    #   referencias de 2016-2026 se solapaban y la figura era ilegible: comprobado mirandola.
    # ⚠ 9-sep-2026: con 6 nodos arriba y 7 abajo hay que ensanchar; con la anchura vieja
    #   las etiquetas se pisaban.  Comprobado mirando el PNG, no el codigo.
    fig, ax = plt.subplots(figsize=(10.6, 4.0))
    ax.set_xlim(-1.5, 10.0)
    ax.set_ylim(-2.35, 3.00)   # -2.95: la fila de abajo pasó a tres lineas y rozaba el rotulo
    ax.axis("off")

    # ⚠ El Forum Math. de Prasad es de 2016 (arXiv:1402.5504), no de 1993 --- la primera version de
    #   esta figura llevaba una fecha inventada de memoria.  Leida de la bibliografia de la nota 1.
    # ⚠ 9-sep-2026.  FALTABA EL NODO QUE MAS IMPORTA.  Cellini-Moseneder Frajria-Papi 2006 da,
    #   en su Cor. 3.2, el criterio Y EL SIGNO en a_P mismo, en todo tipo y en el estadistico de
    #   RAICES --- esta misma recta, no la de corraices de NPP25.  La figura contaba una historia
    #   de dos literaturas que apenas se cruzaron y se dejaba fuera al que ya estaba en la nuestra,
    #   dieciocho anyos antes.  Y abajo faltaba Lecouvey, a quien el propio Albion (2501.18520)
    #   acredita junto a Ayyer-Kumari por la combinatoria de cocientes en tipos clasicos.
    # ⚠ Solo NOMBRE y ANYO.  Lo que cada uno aporto va en el PIE: medido, con las etiquetas
    #   de tres lineas la fuente caia a 3-4 pt EN PAGINA y no se leia.
    arriba = [
        (1976, "Kostant"),
        (2006, "Cellini–Möseneder\nFrajria–Papi"),
        (2009, "Kumar–Lusztig–Prasad"),
        (2016, "Prasad"),
        (2021, "Lübeck–Prasad"),
        (2025, "NPP  ·  Polo"),
    ]
    abajo = [
        (2002, "Adin–Frumkin"),
        # ⚠ 2007, y VA AQUI: la linea de abajo esta ordenada por anyo y el orden es lo que la
        #   figura afirma.  Poner a Lecouvey detras de Petreolle (2016) seria decir que viene
        #   despues, y es al reves: Albion (2501.18520) lo acredita como antecedente junto a
        #   Ayyer-Kumari.
        (2007, "Lecouvey"),
        # ⚠ 2016, no 2015: la bibliografia dice Adv. Appl. Math. 79 (2016).  Los demas nodos de
        #   esta rama van por anyo de PUBLICACION (Ayyer-Kumari 2022, Albion 2023, Kumari II 2024);
        #   dejar este por el anyo del arXiv era mezclar dos convenios en la misma linea temporal.
        (2016, "Pétréolle"),
        # ⚠ Aqui habia dos fechas mal y en el ORDEN equivocado: Ayyer-Kumari figuraba en 2025
        #   y DESPUES de Albion.  Es de 2022 (arXiv:2109.11310v3, J. Algebra), o sea ANTES; y
        #   Albion (arXiv:2212.07343) es de 2023, y Kumari II de 2024.  Leidas de los propios
        #   PDF en _papers/.  Una linea temporal con el orden cambiado no es una errata: dice
        #   quien viene de quien, que es justo lo que esta figura afirma.
        # ⚠ y con seis nodos las etiquetas se pisaban: comprobado mirando el PNG renderizado.
        #   Van partidas en tres lineas cortas para que cada bloque quepa en su hueco.
        (2022, "Ayyer–Kumari"),
        (2023, "Albion"),
        (2024, "Kumari II"),
        (2026, "Marín"),
    ]

    X0, X1 = 0.0, 8.10                      # extremos de las dos lineas horizontales
    ax.plot([X0, X1], [1.35, 1.35], color=AZUL, lw=1.1, zorder=1)
    ax.plot([X0, X1], [-1.35, -1.35], color=NARANJA, lw=1.1, zorder=1)

    # la raiz
    ax.add_patch(Rectangle((-1.30, -0.44), 1.15, 0.88, fc="white", ec="black", lw=1.0, zorder=3))
    ax.text(-0.725, 0.0, "Littlewood\n1940", ha="center", va="center", fontsize=8.6,
            fontweight="bold", zorder=4)
    for y in (1.35, -1.35):
        ax.add_patch(FancyArrowPatch((-0.15, 0.26 * (1 if y > 0 else -1)), (X0, y),
                                     arrowstyle="-", lw=1.1, color=AZUL if y > 0 else NARANJA,
                                     connectionstyle="arc3,rad=%.2f" % (-0.30 if y > 0 else 0.30),
                                     zorder=2))

    def pon(items, y, escalonar=False):
        col = AZUL if y > 0 else NARANJA
        n = len(items)
        for j, (anyo, txt) in enumerate(items):
            x = X0 + (j + 0.62) * (X1 - X0) / (n + 0.1)
            ax.plot([x], [y], "o", ms=6, color=col, zorder=4)
            s = 1 if y > 0 else -1
            # ⚠ escalonado: un nombre de dos lineas y ancho pisaba al vecino.  Alternando la
            #   altura cada etiqueta dispone del doble de anchura.  Visto en el PNG.
            alto = 0.46 + (0.42 if (escalonar and j % 2 == 1) else 0.0)
            if escalonar:
                ax.plot([x, x], [y + s * 0.10, y + s * (alto - 0.04)], lw=0.6, color=col,
                        alpha=0.45, zorder=3)
            ax.text(x, y + s * alto, txt, ha="center", va="bottom" if y > 0 else "top",
                    fontsize=9.2, linespacing=1.25, zorder=5)
            ax.text(x, y - s * 0.17, str(anyo), ha="center", va="top" if y > 0 else "bottom",
                    fontsize=9.0, color=col, fontweight="bold", zorder=5)

    pon(arriba, 1.35, escalonar=True)
    pon(abajo, -1.35)

    # el punto de encuentro
    ax.add_patch(Rectangle((8.55, -0.66), 1.25, 1.32, fc="#eef3f8", ec=AZUL, lw=1.3, zorder=3))
    ax.text(9.175, 0.0, "$a_P^{\\,k}$\nthis note", ha="center", va="center", fontsize=10.5,
            fontweight="bold", color=AZUL, zorder=4)
    for y in (1.35, -1.35):
        ax.add_patch(FancyArrowPatch((X1, y), (8.55, 0.32 * (1 if y > 0 else -1)),
                                     arrowstyle="-|>", mutation_scale=11, lw=1.2,
                                     color=AZUL if y > 0 else NARANJA,
                                     connectionstyle="arc3,rad=%.2f" % (0.30 if y > 0 else -0.30),
                                     zorder=2))

    # EL ROCE.  Adin-Frumkin 2002 se titula "Rim hook tableaux and KOSTANT'S eta-function
    # coefficients": la rama de abajo estaba calculando los coeficientes de la identidad de la rama
    # de arriba, y las dos siguieron sin verse otros veinticuatro anyos.
    xk = X0 + (0 + 0.62) * (X1 - X0) / (len(arriba) + 0.1)      # Kostant 1976
    xa = X0 + (0 + 0.62) * (X1 - X0) / (len(abajo) + 0.1)       # Adin-Frumkin 2002
    ax.plot([xk, xa], [0.92, -0.92], ls=(0, (2, 2)), lw=1.0, color=GRIS, zorder=1)
    ax.text((xk + xa) / 2 - 0.05, 0.02, "their title names\nhis $\\eta$-function", ha="center",
            va="center", fontsize=7.2, color=GRIS, style="italic",
            bbox=dict(fc="white", ec="none", pad=1.6), zorder=3)

    ax.text(4.05, 2.76, "the representation-theoretic line — character values at special elements",
            ha="center", fontsize=8.6, color=AZUL, style="italic")
    ax.text(4.05, -2.22, "the combinatorial line — cores, quotients, and plethysm in classical types",
            ha="center", fontsize=8.6, color=NARANJA, style="italic")
    fig.tight_layout()
    fig.savefig("../note_prasad2/fig_map.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig_map.pdf")


# ---------------------------------------------------------------- 2. EL PORQUE
def fig_capacity():
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.05))

    def cajas(ax, etiquetas, caps, llenos, titulo, sub):
        x = 0.0
        for et, cap, lleno in zip(etiquetas, caps, llenos):
            for s in range(cap):
                fc = AZUL if s < lleno else "white"
                ax.add_patch(Rectangle((x, s * 0.5), 0.72, 0.42, fc=fc, ec="black", lw=0.7))
            ax.text(x + 0.36, -0.30, et, ha="center", va="top", fontsize=8.2)
            x += 1.0
        ax.set_xlim(-0.35, x - 0.1)
        ax.set_ylim(-0.95, max(caps) * 0.5 + 0.62)
        ax.axis("off")
        ax.set_title(titulo, fontsize=9.6, pad=8)
        ax.text((x - 0.28) / 2 - 0.2, max(caps) * 0.5 + 0.28, sub, ha="center", fontsize=8.6,
                color=VERDE if "saturates" in sub else NARANJA, fontweight="bold")

    # tipo A: m=2, n=3 -> 6 beta-numeros, 3 clases de capacidad 2, saturan
    cajas(axes[0], ["$0$", "$1$", "$2$"], [2, 2, 2], [2, 2, 2],
          "type $A$  (Prasad 2016):  $m=2$, $n=3$",
          "$6$ rows, capacity $6$ — saturates")
    # tipo C: m=5, k=2 -> d=2, q=6.  clases 0*,1,2,3*.  cap 1,2,2,1 = 6, filas 5
    cajas(axes[1], ["$0^{*}$", "$1$", "$2$", "$3^{*}$"], [1, 2, 2, 1], [1, 2, 2, 0],
          "type $C$  ($a_P^2$ in $\\mathrm{Sp}(10)$):  $d=2$, $q=6$",
          "$5$ rows, capacity $6$ — one unit of slack")
    axes[1].text(2.05, -0.62, "$*$ = fixed class of the folding, capacity $\\lfloor d/2\\rfloor$",
                 ha="center", fontsize=7.6, color=NARANJA)
    fig.tight_layout()
    fig.savefig("../note_prasad2/fig_capacity.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig_capacity.pdf")


# ---------------------------------------------------------------- 3. EL PLEGADO
def fig_folding():
    q = 8
    fig, ax = plt.subplots(figsize=(6.4, 2.5))
    for r in range(q):
        col = NARANJA if r in (0, q // 2) else AZUL
        ax.plot([r], [0], "o", ms=8, color=col, zorder=3)
        ax.text(r, -0.22, str(r), ha="center", va="top", fontsize=8.6, color=col)
    for r in range(1, q // 2):
        s = q - r
        ax.add_patch(Arc(((r + s) / 2, 0), s - r, 0.85 * (s - r) / (q / 2), theta1=0, theta2=180,
                         color=AZUL, lw=1.1))
        ax.text((r + s) / 2, 0.46 * (s - r) / (q / 2) + 0.04, f"$\\mathrm{{cls}}={r}$",
                ha="center", fontsize=7.8, color=AZUL)
    ax.set_xlim(-0.7, q - 0.3)
    ax.set_ylim(-0.75, 1.35)
    ax.axis("off")
    ax.text(q / 2 - 0.5, -0.62, "the two orange residues are the fixed points of $r\\mapsto q-r$; "
                               "they carry half capacity", ha="center", fontsize=7.8, color=NARANJA)
    fig.tight_layout()
    fig.savefig("../note_prasad2/fig_folding.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig_folding.pdf")


# ---------------------------------------------------------------- 4. LOS TIPOS
def cartan(tipo, n):
    A = [[2 if i == j else 0 for j in range(n)] for i in range(n)]
    if tipo in ("A", "B", "C"):
        for i in range(n - 1):
            A[i][i + 1] = A[i + 1][i] = -1
    if tipo == "B":
        A[n - 1][n - 2] = -2
    if tipo == "C":
        A[n - 2][n - 1] = -2
    if tipo == "F":
        A = [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -2, 2, -1], [0, 0, -1, 2]]
    if tipo == "G":
        A = [[2, -3], [-1, 2]]
    return A


def longs(tipo, n):
    return {"A": [1] * n, "B": [2] * (n - 1) + [1], "C": [1] * (n - 1) + [2],
            "F": [2, 2, 1, 1], "G": [1, 3]}[tipo]


def positivas(A, n):
    simples = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
    R, front = set(simples), list(simples)
    while front:
        nueva = []
        for a in front:
            for i in range(n):
                p, b = 0, list(a)
                while True:
                    b[i] -= 1
                    if min(b) < 0 or tuple(b) not in R:
                        break
                    p += 1
                if p - sum(a[j] * A[i][j] for j in range(n)) >= 1:
                    c = tuple(a[j] + (1 if j == i else 0) for j in range(n))
                    if c not in R:
                        R.add(c)
                        nueva.append(c)
        front = nueva
    return sorted(R)


def fig_types():
    casos = [("B", 3), ("B", 4), ("B", 5), ("C", 3), ("C", 4), ("C", 5), ("F", 4), ("G", 2)]
    fig, ax = plt.subplots(figsize=(8.2, 3.2))
    qmax = 18
    for row, (tipo, n) in enumerate(casos):
        A, d = cartan(tipo, n), longs(tipo, n)
        R = positivas(A, n)
        r = {"B": 2, "C": 2, "F": 2, "G": 3}[tipo]
        for q in range(2, qmax + 1):
            RP = RQ = 0
            sobra = 0
            for c in R:
                nrm = sum(Fraction(c[i] * c[j] * d[i] * A[i][j], 2) for i in range(n) for j in range(n))
                htP = sum(c[i] * d[i] for i in range(n))
                htd = int(Fraction(htP) / nrm)
                if htP % q == 0:
                    RP += 1
                    if htd % q:
                        sobra += 1
                if htd % q == 0:
                    RQ += 1
            if sobra:
                ax.add_patch(Rectangle((q - 0.42, row - 0.40), 0.84, 0.80, fc=NARANJA, ec="none"))
                ax.text(q, row, str(sobra), ha="center", va="center", fontsize=7.4, color="white")
            elif RP:
                ax.add_patch(Rectangle((q - 0.42, row - 0.40), 0.84, 0.80, fc="#cfd9e4", ec="none"))
        ax.text(1.2, row, f"${tipo}_{{{n}}}$   $r={r}$", ha="right", va="center", fontsize=8.6)
    ax.set_xlim(0.0, qmax + 0.7)
    ax.set_ylim(-0.7, len(casos) - 0.3)
    ax.set_yticks([])
    ax.set_xticks(range(2, qmax + 1))
    ax.tick_params(labelsize=7.6)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("$q$", fontsize=9)
    ax.text(qmax / 2 + 1, len(casos) - 0.15,
            "orange = $R_P(q)$ strictly larger than $R(q)$ (number of extra long roots);   "
            "grey = equal and non-empty", ha="center", fontsize=7.8)
    fig.tight_layout()
    fig.savefig("../note_prasad2/fig_types.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig_types.pdf")


# ---------------------------------------------------------------- 5. LOS VALORES
def det_int(M):
    n = len(M)
    if n == 0:
        return 1
    M = [row[:] for row in M]
    sign, prev = 1, 1
    for c in range(n - 1):
        if M[c][c] == 0:
            piv = next((r for r in range(c + 1, n) if M[r][c] != 0), None)
            if piv is None:
                return 0
            M[c], M[piv] = M[piv], M[c]
            sign = -sign
        for r in range(c + 1, n):
            for k in range(c + 1, n):
                M[r][k] = (M[r][k] * M[c][c] - M[r][c] * M[c][k]) // prev
            M[r][c] = 0
        prev = M[c][c]
    return sign * M[n - 1][n - 1]


def det_c(M):
    """Determinante complejo por eliminacion con pivoteo (para fig_union3d)."""
    n = len(M)
    M = [r[:] for r in M]
    d = 1 + 0j
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-11:
            return 0j
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        inv = 1 / M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] * inv
            if f:
                for k in range(c, n):
                    M[r][k] -= f * M[c][k]
    return d


def fig_values():
    import itertools
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.0))
    for ax, (m, titulo) in zip(axes, [(1, None), (5, None)]):
        pass
    datos = {}
    for m in (5, 6):
        q = m + 1

        def H(k, q=q):
            def c(n):
                return 0 if (n < 0 or n % q) else n // q + 1
            return 0 if k < 0 else c(k) - 2 * c(k - 1) + c(k - 2)

        vals = Counter()
        for lam in itertools.product(range(7), repeat=m):
            if any(lam[i] < lam[i + 1] for i in range(m - 1)):
                continue
            M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
                  for j in range(m)] for i in range(m)]
            dd = det_int(M)
            vals[dd // 2] += 1
        datos[m] = vals
    for ax, m in zip(axes, (5, 6)):
        vals = datos[m]
        ks = sorted(vals)
        ax.bar(ks, [vals[k] for k in ks], width=0.9,
               color=[NARANJA if k == 0 else AZUL for k in ks])
        ax.set_title(f"$\\mathrm{{Sp}}({2*m})$,  $a_P^2$ of order ${m+1}$", fontsize=9.4)
        ax.set_xlabel("$\\chi_\\lambda(a_P^2)$", fontsize=8.8)
        ax.set_ylabel("weights", fontsize=8.8)
        ax.tick_params(labelsize=7.4)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        nz = sorted({abs(k) for k in ks if k})
        ax.text(0.5, 0.93, "$|\\chi|\\in\\{" + ",".join(str(v) for v in nz[:11]) +
                (",\\dots" if len(nz) > 11 else "") + "\\}$",
                transform=ax.transAxes, ha="center", fontsize=7.6, color=AZUL)
    fig.tight_layout()
    fig.savefig("../note_prasad2/fig_values.pdf", bbox_inches="tight")
    plt.close(fig)
    print("fig_values.pdf")



# ---------------------------------------------------------------- 6. LA UNION, EN 3D
def fig_union3d():
    """Sp(6): el mismo retículo de pesos a dos capacidades.  Datos recalculados aqui."""
    import itertools
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    m = 3
    TOP = 15

    def h_de(k):
        """h_n del alfabeto {zeta^{+- j k}}_{j=1..m}, t=2m+2, por serie compleja exacta a doble."""
        import cmath
        t = 2 * m + 2
        z = cmath.exp(2j * cmath.pi / t)
        N = 2 * (TOP + m) + 6
        h = [0j] * (N + 1)
        h[0] = 1 + 0j
        for j in range(1, m + 1):
            for e in ((j * k) % t, (-j * k) % t):
                y = z ** e
                acc, nue = 0j, [0j] * (N + 1)
                for n in range(N + 1):
                    acc = acc * y + h[n]
                    nue[n] = acc
                h = nue
        return h

    fig = plt.figure(figsize=(11.6, 4.6))
    for idx, (k, orden) in enumerate([(1, 8), (2, 4)]):
        h = h_de(k)
        H = lambda n: (h[n] if 0 <= n < len(h) else 0j)
        pts = []
        for lam in itertools.product(range(TOP + 1), repeat=m):
            if any(lam[i] < lam[i + 1] for i in range(m - 1)):
                continue
            M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
                  for j in range(m)] for i in range(m)]
            v = det_c(M) / 2
            if abs(v.imag) > 1e-5 or abs(v.real - round(v.real)) > 1e-5:
                continue
            val = int(round(v.real))
            if val:
                pts.append((lam[0], lam[1], lam[2], abs(val)))
        ax = fig.add_subplot(1, 2, idx + 1, projection="3d")
        if pts:
            xs, ys, zs, vs = zip(*pts)
            mx = max(vs)
            ax.scatter(xs, ys, zs, s=[9 + 42 * (v / mx) for v in vs],
                       c=[AZUL if v == 1 else NARANJA for v in vs],
                       depthshade=True, edgecolors="none", alpha=0.85)
        ax.set_xlabel("$\\lambda_1$", fontsize=8, labelpad=-6)
        ax.set_ylabel("$\\lambda_2$", fontsize=8, labelpad=-6)
        ax.set_zlabel("$\\lambda_3$", fontsize=8, labelpad=-6)
        ax.tick_params(labelsize=6, pad=-2)
        mxv = max(v for *_, v in pts) if pts else 0
        ax.set_title(f"$a_P^{{{k}}}$,  order ${orden}$,  capacity ${k}$"
                     f"\n{len(pts)} surviving weights,  max $|\\chi|$ = {mxv}",
                     fontsize=9, pad=-2)
        ax.view_init(elev=20, azim=-58)
        ax.grid(False)
        for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
            pane.pane.set_edgecolor("#dddddd")
            pane.pane.set_alpha(0.25)
    # ⚠ el rotulo lambda_3 del panel DERECHO salia CORTADO por el borde: tight_layout no cuenta
    #   las etiquetas de un eje 3D.  Comprobado mirando el PNG renderizado, no el codigo.
    # ⚠ DOS PROBLEMAS SEGUIDOS, los dos vistos mirando el PNG y no el codigo:
    #   1. el rotulo lambda_3 del panel DERECHO salia CORTADO --- bbox_inches="tight" no rescata
    #      la etiqueta de un eje z en 3D, asi que hay que dejarle sitio a mano;
    #   2. al quitar bbox_inches para arreglar (1), la figura entro en la nota con margenes
    #      muertos y los paneles salieron pequenos.
    #   La combinacion que funciona: dejar sitio a la derecha Y recortar el blanco sobrante.
    fig.tight_layout()
    fig.subplots_adjust(left=0.02, right=0.90, wspace=0.04)
    fig.savefig("../note_prasad2/fig_union3d.pdf", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("fig_union3d.pdf")

for f in (fig_map, fig_capacity, fig_folding, fig_types, fig_values, fig_union3d):
    f()
print("OK")
