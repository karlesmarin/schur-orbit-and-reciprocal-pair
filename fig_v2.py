# -*- coding: utf-8 -*-
# ============================================================================================
#  Las tres figuras del material nuevo de la v2.  14 de agosto de 2026.
#
#  POR QUE.  Los tres bloques de la v2 anaden siete paginas de teoria con CERO figuras, y las tres
#  cosas que cuestan de leer son justo las tres que se dibujan bien:
#
#    fig_increments.pdf   el voraz y el EMPATE: por que |G| <= 2 y donde estan los dos maximizadores
#    fig_reflection.pdf   la REFLEXION T_b = tau - T_a y la Proposicion A1, con su contraejemplo
#    fig_virtual.pdf      caracter GENUINO contra propiamente VIRTUAL, en 3D sobre el retículo de Sp(4)
#
#  TODO SE CALCULA, NADA SE DIBUJA A MANO.  Las clases, los incrementos, los maximizadores, K,
#  g_com, el cuadruple del empate y las expansiones simplecticas salen de las definiciones, en
#  aritmetica entera exacta, dentro de este mismo fichero.  Si una figura se contradice con el
#  texto, es que el texto esta mal.
#
#  PALETA.  Azul #2A78D6 y naranja #EB6834 son los de casa (\stproved y \stverif del paper).  Los
#  otros dos, #00A19A y #B14BC8, se eligieron pasando el validador de paletas: las seis
#  comprobaciones en verde, peor par adyacente DeltaE 10.7 en deuteranopia (objetivo 8).  Y la
#  identidad nunca va solo por color: hay etiqueta directa o marcador en cada serie.
#
#  Authors: Carles Marin, Claude (AI assistant).
#  Run: python fig_v2.py   (escribe los tres PDF en el directorio del paper)
# ============================================================================================

import os
from itertools import combinations, permutations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

AZUL, NARANJA, TEAL, MAGENTA = "#2A78D6", "#EB6834", "#00A19A", "#B14BC8"
BANDA, REGLA, GRIS = "#F2F1EC", "#B9B7AE", "#6B6560"
TINTA = "#2B2B2B"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": REGLA, "axes.linewidth": 0.7,
    "xtick.color": GRIS, "ytick.color": GRIS, "text.color": TINTA,
    "axes.labelcolor": TINTA, "figure.facecolor": "white", "savefig.facecolor": "white",
})


# ------------------------------------------------------------------ la combinatoria -------------
def clases(beta, t):
    C = {}
    for v in beta:
        C.setdefault(v % t, []).append(v)
    for k in C:
        C[k].sort(reverse=True)
    return C


def incrementos(C):
    return [(cs[k] + cs[k + 1], i, k + 1) for i, cs in C.items() for k in range(len(cs) - 1)]


def maximizadores(beta, t, r):
    C = clases(beta, t)
    if len(C) < t:
        return None
    inc = sorted(incrementos(C), key=lambda x: -x[0])
    if len(inc) != 2 * r:
        return None
    tau = inc[r - 1][0]
    may = [x for x in inc if x[0] > tau]
    ig = [x for x in inc if x[0] == tau]
    s = r - len(may)
    if s < 1:
        return None
    sel = [may + ig] if s == len(ig) else [may + list(c) for c in combinations(ig, s)]
    js = []
    for S_ in sel:
        j = {i: 0 for i in C}
        for (_, i, k) in S_:
            j[i] = max(j[i], k)
        js.append(j)
    return js, tau, inc, C


def partes(j, C):
    A, g, L = [], [], []
    for i, cs in C.items():
        k = j.get(i, 0)
        A += cs[:k]
        if k < len(cs):
            g.append(cs[k])
        L += cs[k + 1:]
    return sorted(A, reverse=True), sorted(g, reverse=True), sorted(L, reverse=True)


def anatomia(beta, t, r):
    M = maximizadores(beta, t, r)
    if M is None or len(M[0]) != 2:
        return None
    js, tau, inc, C = M
    pa, pb = partes(js[0], C), partes(js[1], C)
    Ta, Tb = tuple(sorted(pa[0] + pa[2], reverse=True)), tuple(sorted(pb[0] + pb[2], reverse=True))
    S = sorted(v for cs in C.values() if len(cs) >= 2 for v in cs)
    K = sorted(set(Ta) & set(Tb))
    dif = sorted((set(Ta) | set(Tb)) - set(K))
    gcom = sorted(set(S) - set(Ta) - set(Tb))
    refl = (tuple(sorted((tau - v for v in Ta), reverse=True)) == Tb)
    return dict(js=js, tau=tau, inc=inc, C=C, Ta=Ta, Tb=Tb, S=S, K=K, dif=dif, gcom=gcom,
                Ha=pa[0], La=pa[2], Hb=pb[0], Lb=pb[2], refl=refl)


# ------------------------------------------------------------------ Laurent + Sp ----------------
def padd(a, b):
    o = dict(a)
    for e, c in b.items():
        v = o.get(e, 0) + c
        o[e] = v
        if not v:
            del o[e]
    return o


def psc(a, k):
    return {} if k == 0 else {e: c * k for e, c in a.items()}


def sgn(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


PERMS = {}


def perms(n):
    if n not in PERMS:
        PERMS[n] = [(p, sgn(p)) for p in permutations(range(n))]
    return PERMS[n]


def ldiv(num, den):
    q, rem = {}, dict(num)
    le = max(den)
    lc = den[le]
    while rem:
        me = max(rem)
        f = rem[me] // lc
        ke = tuple(x - y for x, y in zip(me, le))
        q[ke] = q.get(ke, 0) + f
        for e, c in den.items():
            t = tuple(x + y for x, y in zip(e, ke))
            v = rem.get(t, 0) - f * c
            rem[t] = v
            if not v:
                del rem[t]
    return {e: c for e, c in q.items() if c}


def det_pares(filas, r):
    out = {}
    for p, s in perms(2 * r):
        e = [0] * r
        for i in range(2 * r):
            a, lo = p[i] >> 1, p[i] & 1
            e[a] += -filas[i] if lo else filas[i]
        k = tuple(e)
        v = out.get(k, 0) + s
        out[k] = v
        if not v:
            del out[k]
    return out


def alt(beta, r):
    N = 2 * r + 2
    tot = {}
    for p in range(N):
        for q in range(p + 1, N):
            if (beta[p] & 1) == (beta[q] & 1):
                continue
            w = (1 if (beta[p] & 1) == 0 else -1) * (-1 if (p + q + 1) & 1 else 1)
            tot = padd(tot, psc(det_pares([beta[i] for i in range(N) if i not in (p, q)], r), w))
    return psc(tot, -2)


def det_binom(a, r):
    out = {}
    for p, s in perms(r):
        for ch in range(1 << r):
            e, c = [0] * r, s
            for i in range(r):
                sg = -1 if (ch >> i) & 1 else 1
                e[i] = sg * a[p[i]]
                if sg < 0:
                    c = -c
            k = tuple(e)
            v = out.get(k, 0) + c
            out[k] = v
            if not v:
                del out[k]
    return out


SPC = {}


def sp_char(mu, r):
    key = (tuple(mu), r)
    if key not in SPC:
        SPC[key] = ldiv(det_binom([mu[j] + r - j for j in range(r)], r),
                        det_binom([r - j for j in range(r)], r))
    return SPC[key]


def expansion(lam, r):
    N = 2 * r + 2
    beta = [lam[i] + (N - 1 - i) for i in range(N)]
    A = alt(beta, r)
    if not A:
        return {}
    P = ldiv(A, alt([N - 1 - i for i in range(N)], r))
    coef, rem = {}, dict(P)
    while rem:
        best = max((tuple(sorted((abs(x) for x in e), reverse=True)) for e in rem),
                   key=lambda d: (sum(d), d))
        c = rem[best]
        coef[best] = coef.get(best, 0) + c
        rem = padd(rem, psc(sp_char(list(best), r), -c))
    return {m: c for m, c in coef.items() if c}


# ================================================================== FIGURA 1 =====================
def fig_increments(beta=(8, 7, 6, 5, 3, 2, 1, 0), t=4, r=2):
    a = anatomia(beta, t, r)
    assert a is not None, "esa beta no da |G| = 2"
    inc, tau, C = a["inc"], a["tau"], a["C"]
    col = {i: c for i, c in zip(sorted(C), [AZUL, NARANJA, TEAL, MAGENTA])}

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.7),
                                 gridspec_kw={"width_ratios": [1.35, 1]})

    # --- izquierda: los 2r incrementos ordenados, con el corte en r
    y = list(range(len(inc)))[::-1]
    for yy, (v, i, k) in zip(y, inc):
        ax.barh(yy, v, height=0.6, color=col[i], edgecolor="white", linewidth=1.2)
        # etiqueta DENTRO de la barra: fuera chocaba con la flecha del empate y con sus guias
        ax.text(v - 0.28, yy, r"$\Delta_{%d}(%d)=%d$" % (i, k, v), va="center", ha="right",
                fontsize=7.5, color="white")
    vmax = max(v for v, _, _ in inc)
    corte = y[r - 1] - 0.5
    # la linea del corte se dibuja SOLO sobre los datos: con axhline cruzaba su propio rotulo a la
    # izquierda y la etiqueta del empate a la derecha
    ax.plot([0, vmax * 1.02], [corte, corte], color=GRIS, lw=1.2, ls=(0, (4, 3)), zorder=2)
    ax.text(-0.4, corte, r"the cut at rank $r=%d$" % r, fontsize=7.5, color=GRIS,
            ha="right", va="center")
    tied = [yy for yy, (v, i, k) in zip(y, inc) if v == tau]
    if len(tied) == 2:
        xarr = vmax * 1.16
        ax.annotate("", xy=(xarr, tied[0]), xytext=(xarr, tied[1]),
                    arrowprops=dict(arrowstyle="<->", color=TINTA, lw=1.1))
        ax.text(xarr + 0.3, sum(tied) / 2, r"tie: $\tau=%d$" % tau, fontsize=8.5,
                color=TINTA, va="center", fontweight="bold")
    ax.set_yticks([])
    ax.set_xticks([x for x in ax.get_xticks() if 0 <= x <= vmax])
    ax.set_xlim(-vmax * 1.02, vmax * 1.66)
    ax.set_xlabel(r"increment $\Delta_i(k)=c_{i,k}+c_{i,k+1}$")
    ax.set_title(r"$\beta=(%s)$,  $t=%d$,  $r=%d$" % (",".join(map(str, beta)), t, r), loc="left")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#E5E3DE", lw=0.6)
    ax.set_axisbelow(True)

    # --- derecha: los dos maximizadores, clase a clase
    bx.set_title("the two maximisers, class by class:\nthey differ in the tied classes only",
                 loc="left", fontsize=8.5)
    ncl = len(C)
    for col_j, (etq, j) in enumerate(zip((r"$g_{\mathrm{a}}$", r"$g_{\mathrm{b}}$"), a["js"])):
        for row, i in enumerate(sorted(C)):
            cs = C[i]
            for pos, v in enumerate(cs):
                k = j[i]
                if pos < k:
                    face, edge, lab = col[i], col[i], "H"
                elif pos == k:
                    face, edge, lab = "white", col[i], "g"
                else:
                    face, edge, lab = BANDA, REGLA, "L"
                x = col_j * 4.6 + pos * 1.02
                bx.add_patch(plt.Rectangle((x, -row - 0.34), 0.94, 0.68, facecolor=face,
                                           edgecolor=edge, linewidth=1.1))
                bx.text(x + 0.47, -row, str(v), ha="center", va="center", fontsize=7,
                        color="white" if lab == "H" else TINTA)
        bx.text(col_j * 4.6 + 1.6, 0.85, etq, ha="center", fontsize=9.5, color=TINTA)
    for row, i in enumerate(sorted(C)):
        bx.text(-0.45, -row, r"$i=%d$" % i, ha="right", va="center", fontsize=8, color=col[i])
    bx.set_xlim(-1.6, 9.0)
    bx.set_ylim(-ncl + 0.2, 1.35)
    bx.axis("off")
    leg = [plt.Rectangle((0, 0), 1, 1, facecolor=GRIS, edgecolor=GRIS,
                         label="top half $H$"),
           plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=GRIS, label="omitted $g_i$"),
           plt.Rectangle((0, 0), 1, 1, facecolor=BANDA, edgecolor=REGLA, label="bottom half $L$")]
    bx.legend(handles=leg, loc="lower center", ncol=3, frameon=False, fontsize=7.5,
              bbox_to_anchor=(0.5, -0.14), handlelength=1.4, columnspacing=1.4)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_increments.pdf"))
    plt.close(fig)
    return a


# ================================================================== FIGURA 2 =====================
def _panel_reflexion(ax, beta, t, r, titulo):
    a = anatomia(beta, t, r)
    assert a is not None
    tau, S = a["tau"], a["S"]
    est = {}
    for v in S:
        est[v] = "K" if v in a["K"] else ("tie" if v in a["dif"] else "gcom")
    cmap = {"K": AZUL, "tie": MAGENTA, "gcom": NARANJA}
    mk = {"K": "o", "tie": "s", "gcom": "D"}
    lo, hi = min(S + [tau - max(S)]), max(S + [tau - min(S)])
    ax.axhline(0, color=REGLA, lw=0.9, zorder=1)
    # los arcos de la reflexion v <-> tau - v, cuando el companero esta en S
    for v in S:
        w = tau - v
        if w in S and w > v:
            m = (v + w) / 2.0
            ax.plot([v, m, w], [0.0, 0.42, 0.0], color=REGLA, lw=0.9, zorder=1,
                    solid_capstyle="round")
    ax.axvline(tau / 2.0, color=GRIS, lw=1.1, ls=(0, (4, 3)), zorder=1)
    ax.text(tau / 2.0 + 0.22, 0.86, r"$\tau/2=%.1f$" % (tau / 2.0), ha="left", fontsize=7.5,
            color=GRIS)
    for v in S:
        ax.scatter([v], [0], s=58, marker=mk[est[v]], color=cmap[est[v]],
                   edgecolor="white", linewidth=1.1, zorder=3)
        ax.text(v, -0.30, str(v), ha="center", fontsize=7, color=TINTA)
    for v, etq in ((min(S), r"$\min\mathcal{S}$"), (max(S), r"$\max\mathcal{S}$")):
        malo = v in a["gcom"]
        ax.annotate(etq, xy=(v, 0.06), xytext=(v, 0.62), ha="center", fontsize=8,
                    color=NARANJA if malo else TINTA,
                    arrowprops=dict(arrowstyle="-", lw=1.0,
                                    color=NARANJA if malo else TINTA), zorder=4)
    ax.set_ylim(-0.55, 1.02)
    ax.set_xlim(lo - 1.2, hi + 1.2)
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_title(titulo, loc="left", fontsize=9)
    return a


def fig_reflection():
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 3.5))
    a1 = _panel_reflexion(axes[0], (8, 7, 6, 5, 3, 2, 1, 0), 4, 2,
                          r"(a)  $[\det]_{D_1}=0$:  $T_{\mathrm{b}}=\tau-T_{\mathrm{a}}$, "
                          r"and both extremes of $\mathcal{S}$ avoid $g_{\mathrm{com}}$")
    # el titulo del panel (b) NO se escribe a mano: se LEE de que extremo cae en g_com.  Al
    # escribirlo a ojo puse "max" y el que cae es el MINIMO; que lo diga el dato.
    b = anatomia((12, 9, 7, 6, 3, 2, 1, 0), 4, 2)
    cual = [n for n, v in ((r"\min", min(b["S"])), (r"\max", max(b["S"]))) if v in b["gcom"]]
    a2 = _panel_reflexion(axes[1], (12, 9, 7, 6, 3, 2, 1, 0), 4, 2,
                          r"(b)  $|G|=2$ but $[\det]_{D_1}\neq0$:  no reflection, "
                          r"and $%s\mathcal{S}$ falls in $g_{\mathrm{com}}$" % "".join(cual))
    leg = [Line2D([], [], marker="o", ls="", color=AZUL, mec="white", ms=7,
                  label=r"$\mathcal{K}=T_{\mathrm{a}}\cap T_{\mathrm{b}}$"),
           Line2D([], [], marker="s", ls="", color=MAGENTA, mec="white", ms=7,
                  label=r"the tie $\{p_1,p_2,q_1,q_2\}$"),
           Line2D([], [], marker="D", ls="", color=NARANJA, mec="white", ms=7,
                  label=r"$g_{\mathrm{com}}$ (omitted by both)"),
           Line2D([], [], color=REGLA, lw=1.0, label=r"$v\leftrightarrow\tau-v$ inside $\mathcal{S}$")]
    fig.legend(handles=leg, loc="lower center", ncol=4, frameon=False, fontsize=7.5,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(OUT, "fig_reflection.pdf"))
    plt.close(fig)
    return a1, a2


# ================================================================== FIGURA 3 =====================
def _panel3d(ax, lam, r, titulo, zlim=None):
    """Barras 3D = magnitud; HUELLA en el suelo = signo.  La huella existe porque sin ella los
    coeficientes negativos, que son el asunto entero del panel (b), quedan invisibles al lado de
    los positivos grandes: la primera version de esta figura fallaba exactamente en eso."""
    c = expansion(list(lam), r)
    z0 = (zlim or (min([0] + list(c.values())), max([0] + list(c.values()))))[0]
    for (mu, v) in c.items():
        col = AZUL if v > 0 else NARANJA
        # la huella de signo, plana, en el plano de abajo
        ax.bar3d(mu[0] - 0.42, mu[1] - 0.42, z0, 0.84, 0.84, 1e-9, color=col,
                 edgecolor=col, linewidth=0.3, shade=False, alpha=0.5)
        ax.bar3d(mu[0] - 0.35, mu[1] - 0.35, 0, 0.7, 0.7, v, color=col,
                 edgecolor="white", linewidth=0.35, shade=True, alpha=0.95)
    for (mu, v) in c.items():
        if v < 0:
            ax.text(mu[0], mu[1], v - 0.35, r"$%d$" % v, color=NARANJA, fontsize=7,
                    ha="center", va="top", zorder=10)
    if zlim:
        ax.set_zlim(*zlim)
    ax.set_xlabel(r"$\mu_1$", labelpad=-6)
    ax.set_ylabel(r"$\mu_2$", labelpad=-6)
    ax.set_zlabel(r"$a_{\lambda\mu}$", labelpad=-8)
    ax.tick_params(labelsize=6.5, pad=-2)
    ax.view_init(elev=24, azim=-58)
    ax.set_title(titulo, loc="left", fontsize=8.5, pad=-2)
    ax.xaxis.pane.set_facecolor("white")
    ax.yaxis.pane.set_facecolor("white")
    ax.zaxis.pane.set_facecolor("white")
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_edgecolor(REGLA)
        axis._axinfo["grid"]["color"] = "#E5E3DE"
        axis._axinfo["grid"]["linewidth"] = 0.5
    return c


def fig_virtual(lam_gen, lam_vir, r=2):
    cg0, cv0 = expansion(list(lam_gen), r), expansion(list(lam_vir), r)
    todos = list(cg0.values()) + list(cv0.values())
    zl = (min(todos + [0]) - 0.6, max(todos + [0]) + 0.4)   # MISMA escala en los dos paneles
    fig = plt.figure(figsize=(7.0, 3.2))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    bx = fig.add_subplot(1, 2, 2, projection="3d")
    cg = _panel3d(ax, lam_gen, r,
                  r"(a)  $\lambda=(%s)$: a genuine character, all $a_{\lambda\mu}>0$"
                  % ",".join(map(str, lam_gen)), zlim=zl)
    cv = _panel3d(bx, lam_vir, r,
                  r"(b)  $\lambda=(%s)$: properly virtual" % ",".join(map(str, lam_vir)), zlim=zl)
    leg = [plt.Rectangle((0, 0), 1, 1, facecolor=AZUL, edgecolor="white", label=r"$a_{\lambda\mu}>0$"),
           plt.Rectangle((0, 0), 1, 1, facecolor=NARANJA, edgecolor="white", label=r"$a_{\lambda\mu}<0$"),
           plt.Rectangle((0, 0), 1, 1, facecolor=REGLA, edgecolor="white", alpha=0.5,
                         label="sign footprint, on the floor")]
    fig.legend(handles=leg, loc="lower center", ncol=2, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(OUT, "fig_virtual.pdf"))
    plt.close(fig)
    return cg, cv


# ================================================================== main =========================
if __name__ == "__main__":
    print("=" * 100)
    print("LAS TRES FIGURAS DE LA v2 -- todo calculado de las definiciones")
    print("=" * 100)

    a = fig_increments()
    print("")
    print("  fig_increments.pdf   beta=(8,7,6,5,3,2,1,0), t=4, r=2")
    print("     incrementos (valor, clase, k) : %s" % a["inc"])
    print("     tau = %d ;  |G| = %d ;  reflexion T_b = tau - T_a : %s"
          % (a["tau"], len(a["js"]), a["refl"]))

    a1, a2 = fig_reflection()
    print("")
    print("  fig_reflection.pdf")
    for et, x in (("(a)", a1), ("(b)", a2)):
        print("     %s tau=%-4d K=%-22s g_com=%-14s empate=%-20s reflexion=%s  extremos en g_com: %s"
              % (et, x["tau"], str(x["K"]), str(x["gcom"]), str(x["dif"]), x["refl"],
                 [v for v in (min(x["S"]), max(x["S"])) if v in x["gcom"]]))

    # se ELIGEN los dos lambda, no se ponen a dedo: el mas rico de cada tipo en el rango
    # criterio: el VIRTUAL se elige por masa negativa (que se vea lo que la figura ensena) y el
    # GENUINO por rango de coeficientes parecido, para que los dos paneles compartan escala z sin
    # que uno quede aplastado.  Se imprime el criterio y los dos elegidos.
    mejor_gen = mejor_vir = None
    cands_gen = []
    from itertools import combinations_with_replacement
    for lam in combinations_with_replacement(range(9, -1, -1), 6):
        c = expansion(list(lam), 2)
        if len(c) < 8:
            continue
        # "propiamente virtual" = LOS DOS signos.  Todo negativo es -(genuino), no virtual:
        # mi primer criterio decia "tiene algun negativo" y eligio (9,9,9,6,1,0), que es de esos.
        pos = [v for v in c.values() if v > 0]
        neg = [v for v in c.values() if v < 0]
        if pos and neg:
            clave = (min(len(pos), len(neg)), -min(c.values()), len(c))
            if mejor_vir is None or clave > mejor_vir[0]:
                mejor_vir = (clave, lam, c)
        elif pos and not neg:
            cands_gen.append((lam, c))
    # el GENUINO se elige DESPUES, con rango de coeficientes parecido al del virtual: los dos
    # paneles comparten escala z, y sin eso uno aplasta al otro y la figura deja de decir nada.
    rng = max(abs(v) for v in mejor_vir[2].values())
    mejor_gen = min(cands_gen, key=lambda x: (abs(max(x[1].values()) - rng), -len(x[1])))
    mejor_gen = (None, mejor_gen[0], mejor_gen[1])
    print("")
    print("     criterio -- virtual: LOS DOS signos, luego (min(#+,#-), |mas negativo|, #terminos);")
    print("                 genuino: rango de coeficientes mas parecido al del virtual (escala z comun)")
    cg, cv = fig_virtual(mejor_gen[1], mejor_vir[1], 2)
    print("")
    print("  fig_virtual.pdf   (los dos lambda se ELIGEN por numero de terminos, no a dedo)")
    print("     genuino : lambda=%s -> %s" % (str(mejor_gen[1]), sorted(cg.items(), reverse=True)))
    print("     virtual : lambda=%s -> %s" % (str(mejor_vir[1]), sorted(cv.items(), reverse=True)))
    print("")
    for f in ("fig_increments.pdf", "fig_reflection.pdf", "fig_virtual.pdf"):
        p = os.path.join(OUT, f)
        print("     %-24s %d bytes" % (f, os.path.getsize(p)))
    print("DONE")
