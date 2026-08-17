# -*- coding: utf-8 -*-
# ============================================================================================
#  Las dos figuras de los huecos que quedaban: la introduccion (172 lineas sin nada visual) y
#  Background (80).  14 de agosto de 2026.
#
#    fig_alphabet.pdf  EL OBJETO: la orbita de raices CONGELADA mas r pares reciprocos LIBRES,
#                      en el plano complejo, y por que el par no esta en ninguna orbita.
#    fig_beta.pdf      EL DICCIONARIO que usa el paper entero: lambda -> beta -> clases mod t,
#                      con el conteo de exceso sum_i (n_i - 1) = 2r hecho, no afirmado.
#
#  Todo calculado.  El perfil, el conteo de exceso y el determinante del alfabeto salen de las
#  definiciones, y un assert cierra cada uno contra lo que el paper afirma.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python fig_intro.py
# ============================================================================================

import cmath
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AZUL, NARANJA, TEAL, MAGENTA = "#2A78D6", "#EB6834", "#00A19A", "#B14BC8"
BANDA, REGLA, GRIS, TINTA = "#F2F1EC", "#B9B7AE", "#6B6560", "#2B2B2B"
CLASE = [AZUL, NARANJA, TEAL, MAGENTA, "#8A6D3B", "#3B6E8A"]

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "text.color": TINTA, "figure.facecolor": "white", "savefig.facecolor": "white",
})


# ================================================================== fig_alphabet =================
def panel_alfabeto(ax, t, pares, titulo, nota):
    """pares = lista de z complejos; se dibuja z y 1/z para cada uno."""
    th = [cmath.exp(2j * cmath.pi * k / t) for k in range(t)]
    ang = [i / 400 * 2 * math.pi for i in range(401)]
    ax.plot([math.cos(a) for a in ang], [math.sin(a) for a in ang], color=REGLA, lw=0.9)
    ax.axhline(0, color="#E5E3DE", lw=0.6, zorder=0)
    ax.axvline(0, color="#E5E3DE", lw=0.6, zorder=0)
    for k, w in enumerate(th):
        ax.scatter([w.real], [w.imag], s=95, color=CLASE[k % len(CLASE)], edgecolor="white",
                   linewidth=1.2, zorder=4, marker="o")
    for n, z in enumerate(pares):
        for w, et in ((z, r"$z_{%d}$" % (n + 1)), (1 / z, r"$z_{%d}^{-1}$" % (n + 1))):
            ax.scatter([w.real], [w.imag], s=95, color=TINTA, edgecolor="white", linewidth=1.2,
                       zorder=5, marker="s")
            # la etiqueta del inverso se aparta hacia FUERA del origen: pegada al centro chocaba
            # con los puntos congelados +-1
            dx = 7 if abs(w) > 1 else -8
            ax.annotate(et, (w.real, w.imag), textcoords="offset points",
                        xytext=(dx, 6 if w.imag >= 0 else -11),
                        ha="left" if dx > 0 else "right", fontsize=7.5, color=TINTA)
        ax.plot([z.real, (1 / z).real], [z.imag, (1 / z).imag], color=TINTA, lw=0.8,
                ls=(0, (2, 2)), zorder=3)
    det = 1
    for w in th:
        det *= w
    for z in pares:
        det *= z * (1 / z)
    ax.set_title(titulo, loc="left", fontsize=8.2)
    ax.text(0.5, -1.72, nota, ha="center", fontsize=7, color=GRIS, transform=ax.transData)
    ax.text(0.5, -2.06, r"$\prod(\mathrm{alphabet})=(-1)^{t-1}=%s$"
            % ("+1" if abs(det - 1) < 1e-9 else "-1"),
            ha="center", fontsize=7.5, color=TINTA, transform=ax.transData)
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-2.35, 1.55)
    ax.set_aspect("equal")
    ax.axis("off")
    return det


figA = plt.figure(figsize=(7.0, 2.75))
# Paso 0.328 y no 0.331: con el anterior, el pie del tercer panel se salia 1.4 pt de la figura y
# matplotlib recortaba el ultimo glifo, que es justo el exponente de $O(6)^{-}$.  Se leia "O(6)",
# o sea la componente identidad, que es lo contrario de lo que la frase dice.
axes = [figA.add_axes([0.01 + i * 0.328, 0.06, 0.31, 0.88]) for i in range(3)]
d1 = panel_alfabeto(axes[0], 4, [cmath.exp(0.9j)],
                    r"(a)  $t=4$, one pair, $|z|=1$",
                    # No "off the orbit": at z = zeta^j the pair DOES meet mu_t.  The caption was
                    # corrected for exactly this and the drawn label had kept the absolute claim.
                    r"on the circle, and free of the orbit")
d2 = panel_alfabeto(axes[1], 4, [complex(1.75, 0.0)],
                    r"(b)  $t=4$, one pair, $z$ real",
                    r"$z$ outside, $z^{-1}$ inside: a free direction")
d3 = panel_alfabeto(axes[2], 2, [complex(2.0, 0.8), cmath.exp(2.5j) * 1.35],
                    r"(c)  $t=2$, two pairs",
                    r"an $r$-dimensional free torus, inside $O(6)^{-}$")
leg = [Line2D([], [], marker="o", ls="", color=GRIS, mec="white", ms=8,
              label=r"the frozen orbit $\mu_t$"),
       Line2D([], [], marker="s", ls="", color=TINTA, mec="white", ms=7,
              label=r"a free reciprocal pair $z,\,z^{-1}$")]
figA.legend(handles=leg, loc="lower center", ncol=2, frameon=False, fontsize=7.5,
            bbox_to_anchor=(0.5, -0.015))
figA.savefig(os.path.join(OUT, "fig_alphabet.pdf"))
plt.close(figA)
# El producto del alfabeto es prod(mu_t) * prod(z * 1/z) = zeta^{t(t-1)/2} = (-1)^{t-1}.  Mi primera
# expectativa aqui decia +1 para t=4, y este assert la tumbo: es -1 para TODO t par, no solo para
# t=2.  La figura tenia razon y yo no; el enunciado que sale es el general.
for tt, dd in ((4, d1), (4, d2), (2, d3)):
    assert abs(dd - (-1) ** (tt - 1)) < 1e-9, "*** el determinante no da (-1)^{t-1} ***"


# ================================================================== fig_beta =====================
LAM, T = (5, 3, 2, 1, 1, 0), 4
N = len(LAM)
BETA = [LAM[j] + (N - 1 - j) for j in range(N)]
PERFIL = tuple(sum(1 for b in BETA if b % T == i) for i in range(T))
R = (N - T) // 2
EXCESO = sum(n - 1 for n in PERFIL)

figB = plt.figure(figsize=(7.0, 2.5))
p = figB.add_axes([0.015, 0.10, 0.26, 0.80])
q = figB.add_axes([0.315, 0.10, 0.30, 0.80])
s = figB.add_axes([0.665, 0.10, 0.325, 0.80])

# --- (a) el diagrama de Young
for i, li in enumerate(LAM):
    for k in range(li):
        p.add_patch(plt.Rectangle((k, -i), 0.92, 0.92, facecolor=BANDA, edgecolor=REGLA, lw=0.8))
    p.text(-0.35, -i + 0.46, r"$\lambda_{%d}=%d$" % (i + 1, li), ha="right", va="center",
           fontsize=7, color=GRIS)
p.set_xlim(-2.5, max(LAM) + 0.4)
p.set_ylim(-N + 0.4, 1.5)
p.axis("off")
p.set_title(r"(a) $\lambda=(%s)$" % (",".join(map(str, LAM))), loc="left", fontsize=8.2)

# --- (b) beta_j = lambda_j + N - j
q.set_title(r"(b) $\beta_j=\lambda_j+N-j$, strictly decreasing", loc="left", fontsize=8.2)
for j in range(N):
    y = -j
    q.text(-0.1, y, r"$%d+%d$" % (LAM[j], N - 1 - j), ha="right", va="center", fontsize=7,
           color=GRIS)
    q.annotate("", xy=(0.75, y), xytext=(0.12, y),
               arrowprops=dict(arrowstyle="-|>", color=REGLA, lw=0.8))
    q.scatter([1.05], [y], s=200, color=CLASE[BETA[j] % T], edgecolor="white", linewidth=1.1,
              zorder=3)
    q.text(1.05, y, str(BETA[j]), ha="center", va="center", fontsize=7, color="white", zorder=4)
    q.text(1.45, y, r"$\equiv %d$" % (BETA[j] % T), ha="left", va="center", fontsize=7,
           color=CLASE[BETA[j] % T])
q.set_xlim(-1.5, 2.4)
q.set_ylim(-N + 0.4, 1.5)
q.axis("off")

# --- (c) las clases y el conteo de exceso
s.set_title(r"(c) the classes mod $t=%d$, and the excess" % T, loc="left", fontsize=8.2)
for i in range(T):
    xs = [b for b in BETA if b % T == i]
    for k, b in enumerate(xs):
        s.scatter([i * 1.0], [-k * 0.75], s=200, color=CLASE[i], edgecolor="white", linewidth=1.1,
                  zorder=3)
        s.text(i * 1.0, -k * 0.75, str(b), ha="center", va="center", fontsize=7, color="white",
               zorder=4)
    s.text(i * 1.0, 0.62, r"$n_{%d}=%d$" % (i, PERFIL[i]), ha="center", fontsize=7.5,
           color=CLASE[i])
    if PERFIL[i] >= 2:
        s.add_patch(plt.Rectangle((i - 0.34, -(PERFIL[i] - 1) * 0.75 - 0.34), 0.68,
                                  (PERFIL[i] - 1) * 0.75 + 0.68, facecolor="none",
                                  edgecolor=TINTA, lw=1.0, ls=(0, (3, 2)), zorder=2))
s.text((T - 1) / 2.0, -2.35,
       r"$\sum_i(n_i-1)=%d=2r$,  so $r=%d$" % (EXCESO, R), ha="center", fontsize=8.5, color=TINTA)
s.text((T - 1) / 2.0, -2.85,
       r"the dashed classes are the excess $\mathcal{E}$; $|\mathcal{E}|=e=%d$"
       % sum(1 for n in PERFIL if n >= 2), ha="center", fontsize=7, color=GRIS)
s.set_xlim(-0.8, T - 0.2)
s.set_ylim(-3.2, 1.15)
s.axis("off")

figB.savefig(os.path.join(OUT, "fig_beta.pdf"))
plt.close(figB)
assert EXCESO == 2 * R, "*** el conteo de exceso no cuadra con N = t + 2r ***"

# ================================================================== fig_thread ===================
#  EL HILO: el paper como una cadena de preguntas que se van respondiendo, y cada respuesta abre la
#  siguiente.  No es el mapa de dependencias de la seccion 10 -- aquel dice QUE SE APOYA EN QUE;
#  este dice POR QUE SE PREGUNTA LO SIGUIENTE.  Cada perla lleva su estatus.
#
#  Nada inventado: cada respuesta es un enunciado del paper, citado por su nombre.
#  Y cada perla dice DONDE se responde: asi la figura no resume, GUIA la lectura.
#  Y estan redactadas para que FLUYAN: cada pregunta recoge una palabra de la respuesta anterior
#  --- orbit, factorises, triple, vanishes, one pair, virtual --- de modo que el collar se lea como
#  una sola conversacion y no como siete titulares sueltos.
PERLAS = [
    ("What does a Schur polynomial\ndo on a full orbit $\\mu_t$?",
     "Littlewood: it vanishes exactly\nwhen the $t$-core is nonempty", "ext", "classical"),
    ("And if one free reciprocal\npair joins that orbit?",
     "then it factorises: the value is\na triple $(d_1,d_2,d_3)$ and a sign", "ok", "§3"),
    ("Is that triple all the\nshape leaves behind?",
     "yes — a partition of any size\nis compressed to three integers", "ok", "§3"),
    # "the set where it vanishes does, for every r" read as a claim about every t, which is exactly
    # what Conjecture 8.43 leaves open for even t >= 4.  The pearl that follows already asks the
    # t = 2 question, so the answer says t = 2 too and the chain stops promising the whole plane.
    ("Does the factorising\nsurvive a second pair?",
     "no. But at $t=2$ the set where\nit vanishes does, for every $r$", "ok", "§7–8"),
    ("At $t=2$, which $\\lambda$\nvanish for every $r$?",
     "two families, and nothing else", "ok", "§8"),
    ("And for every $t$\nat once?",
     "one reflection of the excess part\nforces it — for all $t$ and $r$", "ok", "§8"),
    ("Why did one pair\nbehave so differently?",
     "there the folded expansion is\n$\\pm$ genuine; from two on it need not be", "ok", "§8"),
    ("So what is\nstill missing?",
     "one set: is the common omitted\npart $g_{\\mathrm{com}}$ symmetric too?", "abierto", "§8"),
]
CP = {"ok": AZUL, "ext": GRIS, "abierto": TINTA}

# UN SOLO HILO, y las perlas colgando de el: la catenaria hace que se lea de un tiron y quita el
# salto de linea, que en la version de dos filas obligaba a un arco que cruzaba por encima de dos
# perlas.  Los rotulos alternan arriba y abajo, que es lo que da sitio para siete en un ancho.
figT = plt.figure(figsize=(7.0, 3.9))
tx = figT.add_axes([0.005, 0.03, 0.99, 0.94])
X0, X1, B, M = 0.55, 11.45, 5.2, 1.05
CX = (X0 + X1) / 2.0


def hilo(x):
    return M * (math.cosh((x - CX) / B) - 1.0)


xs = [X0 + i * (X1 - X0) / 300 for i in range(301)]
tx.plot(xs, [hilo(x) for x in xs], color=REGLA, lw=1.3, zorder=1, solid_capstyle="round")
for i, (q, a, st, donde) in enumerate(PERLAS):
    x = X0 + i * (X1 - X0) / (len(PERLAS) - 1)
    y = hilo(x)
    arriba = (i % 2 == 0)
    s = 1 if arriba else -1
    tx.scatter([x], [y], s=290, color=CP[st], edgecolor="white", linewidth=1.5, zorder=4,
               marker="o" if st != "abierto" else "D")
    tx.text(x, y, str(i + 1), ha="center", va="center", fontsize=7.5, color="white", zorder=5)
    tx.plot([x, x], [y + s * 0.26, y + s * 0.52], color=REGLA, lw=0.7, zorder=2)
    va_q = "bottom" if arriba else "top"
    tx.text(x, y + s * 0.60, q, ha="center", va=va_q, fontsize=6.3, color=GRIS, style="italic")
    tx.text(x, y + s * (1.32 if arriba else 0.98), a, ha="center", va=va_q, fontsize=6.5,
            color=TINTA)
    tx.text(x, y + s * (2.06 if arriba else 1.72), donde, ha="center", va=va_q, fontsize=6.5,
            color=CP[st])
    if i:
        xp = X0 + (i - 1) * (X1 - X0) / (len(PERLAS) - 1)
        xm = (xp + x) / 2.0
        tx.annotate("", xy=(xm + 0.22, hilo(xm + 0.22)), xytext=(xm - 0.22, hilo(xm - 0.22)),
                    arrowprops=dict(arrowstyle="-|>", color=GRIS, lw=1.0), zorder=3)
legT = [Line2D([], [], marker="o", ls="", color=AZUL, mec="white", ms=8, label="proved here"),
        Line2D([], [], marker="o", ls="", color=GRIS, mec="white", ms=8, label="classical"),
        Line2D([], [], marker="D", ls="", color=TINTA, mec="white", ms=7, label="open")]
tx.legend(handles=legT, loc="lower center", ncol=3, frameon=False, fontsize=7.5,
          bbox_to_anchor=(0.5, -0.055))
tx.set_xlim(-1.65, 13.65)
tx.set_ylim(-2.55, 3.55)
tx.axis("off")
figT.savefig(os.path.join(OUT, "fig_thread.pdf"))
plt.close(figT)

print("  fig_thread.pdf  --  el hilo de %d preguntas, con estatus %s"
      % (len(PERLAS), [st for _, _, st, _ in PERLAS]))
print("")
print("  fig_alphabet.pdf  --  el determinante del alfabeto, CALCULADO:")
print("     (a) t=4, un par : %+d      (b) t=4, un par : %+d      (c) t=2, dos pares : %+d"
      % (round(d1.real), round(d2.real), round(d3.real)))
print("     el paper dice que a t=2 el alfabeto vive en la componente det = -1 de O(2r+2): CUADRA")
print("")
print("  fig_beta.pdf  --  el diccionario, con el conteo hecho:")
print("     lambda = %s   ->   beta = %s" % (str(LAM), str(BETA)))
print("     perfil (n_0..n_{t-1}) = %s ;  sum (n_i - 1) = %d ;  2r = %d ;  e = %d"
      % (str(PERFIL), EXCESO, 2 * R, sum(1 for n in PERFIL if n >= 2)))
print("     N = t + 2r : %d = %d + 2*%d  CUADRA" % (N, T, R))
for f in ("fig_alphabet.pdf", "fig_beta.pdf"):
    print("     %-20s %d bytes" % (f, os.path.getsize(os.path.join(OUT, f))))
print("DONE")
