# -*- coding: utf-8 -*-
# EL NUMERADOR GKRS, DIBUJADO: un transversal de las clases plegadas.   16 de agosto de 2026.
#
# POR QUE ESTA FIGURA.  fig_division dibuja el COCIENTE y da por supuesto que el lector ya se
# imagina el numerador.  Pero el numerador es, por prop:transversal, un objeto que se ve: se
# reparten los R' residuos de v = 2 Lambda + 2 rho en las m'+1 clases plegadas modulo t, y un
# sumando de nu es elegir UN indice de cada clase no nula y NINGUNO de la clase 0.  Con eso, las
# tres partes con contenido de la proposicion se leen de un vistazo:
#
#     (ii)  |supp nu| = 2 prod_j n_j        -- una eleccion por columna, y dos quiralidades
#     (iii) nu = 0  <=>  una columna vacia  -- el criterio de anulacion, sin calcular nada
#     (v)   el top lo da el punto MAS BAJO de cada columna (el v mas pequeno)
#
# CONTROLES -- la figura SE NIEGA A DIBUJARSE si alguno falla
#   C1  para los dos Lambda dibujados, |supp nu| por enumeracion honesta de W^1 == 2 prod n_j.
#   C2  el transversal resaltado es admisible (tau^B != 0) y el conjunto NO transversal que se
#       dibuja tachado no lo es.
#   C3  el Lambda del panel derecho tiene nu vacio, y por la razon que se anuncia: una clase sin tocar.
#   C4  S_min (el punto mas bajo de cada columna) da el maximo de supp nu.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python fig_transversal.py

import itertools
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SECOND = "#52514e"; MUTED = "#898781"
PROVED = "#2a78d6"; VERIF = "#eb6834"; DEAD = "#b9b7ae"

plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"],
                     "mathtext.fontset": "dejavuserif", "font.size": 9,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
                     "savefig.facecolor": SURFACE})

T, R = 7, 2
MP, RP = (T - 1) // 2, (T - 1) // 2 + R


def plegar(v, t):
    v %= t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


def sgn_perm(perm):
    n, s, visto = len(perm), 1, [False] * len(perm)
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
    return s


def jacobi(a, n):
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


def tauB(A, t, mp):
    cl, ep = [], []
    for v in A:
        c_, e_ = plegar(v, t)
        cl.append(c_)
        ep.append(e_)
    if sorted(cl) != list(range(1, mp + 1)):
        return 0
    s = sgn_perm([c - 1 for c in cl])
    for e in ep:
        s *= e
    return int(jacobi((-2) % t, t) ** ((t + 3) // 2) * s)


def uve(Lam):
    return [2 * Lam[i] + 2 * (RP - i) - 1 for i in range(RP)]


def clases(v):
    d = {}
    for i, x in enumerate(v):
        d.setdefault(plegar(x, T)[0], []).append(i)
    return d


def nu_honesto(Lam):
    """nu por enumeracion de W^1, sin usar la lectura de transversales -- es el control."""
    v = uve(Lam)
    out = {}
    for perm in itertools.permutations(range(RP)):
        s = sgn_perm(list(perm))
        base = [v[perm[i]] for i in range(RP)]
        for eps in itertools.product((1, -1), repeat=RP):
            u = [base[i] * eps[i] for i in range(RP)]
            if not (all(u[i] > u[i + 1] for i in range(MP - 1)) and u[MP - 1] > 0):
                continue
            f = u[MP:]
            if not (all(f[i] > f[i + 1] for i in range(R - 1)) and f[R - 2] > abs(f[R - 1])):
                continue
            sg = s
            for e in eps:
                sg *= e
            tv = tauB(u[:MP], T, MP)
            if tv:
                out[tuple(u[MP:])] = out.get(tuple(u[MP:]), 0) + sg * tv
    return {k: val for k, val in out.items() if val}


# ------------------------------------------------------------------ elegir los dos ejemplos
VIVO = MUERTO = None
for Lam in sorted(itertools.product(range(4), repeat=RP), key=lambda L: sum(L)):
    if any(Lam[i] < Lam[i + 1] for i in range(RP - 1)):
        continue
    d = clases(uve(Lam))
    lleno = all(len(d.get(j, [])) for j in range(1, MP + 1))
    if lleno and VIVO is None and sum(len(d.get(j, [])) for j in range(1, MP + 1)) >= 4:
        VIVO = Lam
    if not lleno and MUERTO is None:
        MUERTO = Lam
    if VIVO and MUERTO:
        break
if VIVO is None or MUERTO is None:
    print("no se encontraron los dos ejemplos"); sys.exit(1)

fallos = []
for Lam in (VIVO, MUERTO):
    d = clases(uve(Lam))
    pred = 2
    for j in range(1, MP + 1):
        pred *= len(d.get(j, []))
    real = len(nu_honesto(Lam))
    if real != pred:
        fallos.append(("C1", Lam, pred, real))
if nu_honesto(MUERTO):
    fallos.append(("C3", MUERTO, "nu deberia ser vacio"))

v_v, d_v = uve(VIVO), clases(uve(VIVO))
S_min = [max(d_v[j]) for j in range(1, MP + 1)]
S_alt = [min(d_v[j]) for j in range(1, MP + 1)]
if tauB(sorted([v_v[i] for i in S_min], reverse=True), T, MP) == 0:
    fallos.append(("C2", "S_min no admisible"))
mala = None
for j in range(1, MP + 1):
    if len(d_v[j]) >= 2:
        mala = [d_v[j][0], d_v[j][1]] + [max(d_v[k]) for k in range(1, MP + 1) if k != j][:MP - 2]
        break
if mala and tauB(sorted([v_v[i] for i in mala], reverse=True), T, MP) != 0:
    fallos.append(("C2", "el conjunto NO transversal sale admisible"))

nu_v = nu_honesto(VIVO)
top = max(nu_v, key=lambda k: (sum(k), k))
Sc = [i for i in range(RP) if i not in S_min]
libre = sorted([v_v[i] for i in Sc], reverse=True)
if tuple(libre) != top:
    fallos.append(("C4", "S_min no da el maximo", libre, list(top)))

if fallos:
    print("LA FIGURA NO SE DIBUJA. Controles fallidos:")
    for f in fallos:
        print("   ", f)
    sys.exit(1)
print("  controles C1 C2 C3 C4 ok")
print("  vivo   Lambda=%s  v=%s  clases=%s  |supp nu|=%d" % (VIVO, v_v, {k: v for k, v in sorted(d_v.items())}, len(nu_v)))
print("  muerto Lambda=%s  v=%s  clases=%s" % (MUERTO, uve(MUERTO), {k: v for k, v in sorted(clases(uve(MUERTO)).items())}))

# ------------------------------------------------------------------ dibujo
fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.0))


def panel(ax, Lam, resaltar, titulo, pie):
    v, d = uve(Lam), clases(uve(Lam))
    for j in range(0, MP + 1):
        col = DEAD if j == 0 else SECOND
        ax.text(j, -0.75, ("class $0$" if j == 0 else "$%d$" % j), ha="center", va="top",
                fontsize=9, color=col)
        ax.plot([j, j], [-0.45, RP - 0.4], color="#e8e6e1", linewidth=10, zorder=0,
                solid_capstyle="round")
        if j and not d.get(j):
            ax.text(j, (RP - 1) / 2.0, "empty", ha="center", va="center", fontsize=9.5,
                    color=VERIF, rotation=90, style="italic")
    for i, x in enumerate(v):
        j = plegar(x, T)[0]
        dentro = resaltar is not None and i in resaltar
        color = PROVED if dentro else (DEAD if j == 0 else MUTED)
        ax.scatter([j], [RP - 1 - i], s=150 if dentro else 90, zorder=3,
                   facecolors=color, edgecolors="none")
        ax.text(j + 0.17, RP - 1 - i, "$V_%d\\!=\\!%d$" % (i + 1, x), ha="left", va="center",
                fontsize=7.6, color=INK if dentro else SECOND)
    ax.set_xlim(-0.62, MP + 0.72)
    ax.set_ylim(-1.15, RP - 0.25)
    ax.set_title(titulo, fontsize=9.5, color=INK)
    ax.text(0.5, -0.235, pie, transform=ax.transAxes, ha="center", va="top",
            fontsize=8.2, color=SECOND)
    ax.axis("off")


nj = [len(d_v[j]) for j in range(1, MP + 1)]
panel(axes[0], VIVO, set(S_min),
      "a transversal: one index per nonzero class, none from class $0$",
      r"$n_j=%s$, so $|\mathrm{supp}\,\nu|=2\prod_j n_j=%d$."
      r"  Filled: the transversal $S_{\min}$ of lowest $V$, which carries the top weight."
      % (",".join(str(x) for x in nj), len(nu_v)))
# NOTA: nada de numeros de proposicion aqui dentro.  Un numero escrito a mano en una figura envejece
# en silencio en cuanto se inserta un enunciado antes -- ya nos paso con los "§" de fig_thread.  La
# referencia va en el pie de LaTeX, que se renumera solo.
panel(axes[1], MUERTO, None,
      "one class unhit, and the whole numerator is zero",
      r"$\nu\equiv0$: no transversal exists, so nothing is left to divide."
      "\n" r"Read off the residues of $V=2(\Lambda+\rho)$; no character is evaluated.")

for a in axes:
    for s in a.spines.values():
        s.set_visible(False)
fig.tight_layout()
fig.subplots_adjust(bottom=0.22)
fig.savefig("fig_transversal.pdf")
print("  fig_transversal.pdf escrito")
