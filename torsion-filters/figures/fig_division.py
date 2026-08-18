# -*- coding: utf-8 -*-
# LA UNICA DIFICULTAD QUE QUEDA EN (L1): LA DIVISION POR Delta_t.   16 de agosto de 2026.
#
# POR QUE ESTA FIGURA Y NO OTRA.  Tras el Lema de inyectividad, el numerador nu vive en {0,+-1} POR
# TEOREMA.  Lo unico que sigue sin explicacion es por que dividir por  prod_j (z_j^{t/2}-z_j^{-t/2})
# no fabrica un coeficiente 2.  Y eso es exactamente lo que la formula NO ensena: nu tiene unos
# pocos terminos y c tiene muchos, y aun asi los dos viven en {+-1}.  c es una escalera con signo
# cuya diferencia finita es nu.  Eso se ve dibujado.
#
# TODO EL CALCULO ES ENTERO Y EN PYTHON PURO.  tau^B_t tiene forma cerrada (Cor. del filtro impar
# con su signo), luego nu no necesita ni caracteres ni Sage; y c se obtiene de nu por la sustitucion
# hacia atras, que es el objeto que la figura quiere ensenar.
#
# CONTROLES -- la figura SE NIEGA A DIBUJARSE si alguno falla
#   C1  para cada Lambda dibujado:  c . Delta_t == nu  exactamente.
#   C2  todos los coeficientes de c y de nu son +-1.
#   C3  el c de Lambda=(6,6,6) coincide con el que calculo SAGE por branching, guardado en
#       ../gates/gkrs_L1_DUMP.json.  Sin este cruce la figura solo se comprobaria a si misma.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python fig_division.py

import itertools
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, GRID = "#1b1b1b", "#8a8a8a", "#dedede"
PLUS, MINUS, SECOND = "#1f6f8b", "#c1462f", "#5a5a5a"


# ------------------------------------------------------------------ aritmetica de W(D_r)
def enderezar(x):
    r = len(x)
    a = [abs(int(v)) for v in x]
    if len(set(a)) != r:
        return None
    idx = sorted(range(r), key=lambda i: -a[i])
    s, visto = 1, [False] * r
    for i in range(r):
        if visto[i]:
            continue
        j, L = i, 0
        while not visto[j]:
            visto[j] = True
            j = idx[j]
            L += 1
        if L % 2 == 0:
            s = -s
    y = [int(x[i]) for i in idx]
    neg = sum(1 for v in y if v < 0)
    cero = any(v == 0 for v in y)
    y = [abs(v) for v in y]
    if neg % 2 == 1 and not cero:
        y[-1] = -y[-1]
    return (tuple(y), s)


def desplazar(x, paso):
    r = len(x)
    out = {}
    for eps in itertools.product((1, -1), repeat=r):
        sg = 1
        for e in eps:
            sg *= e
        e2 = enderezar(tuple(int(x[j]) + paso * eps[j] for j in range(r)))
        if e2 is None:
            continue
        out[e2[0]] = out.get(e2[0], 0) + sg * e2[1]
    return {k: v for k, v in out.items() if v != 0}


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


def plegar(v, t):
    v %= t
    if v == 0:
        return (0, 1)
    if 2 * v == t:
        return (t // 2, 1)
    return (v, 1) if v < t - v else (t - v, -1)


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
    """forma cerrada: tau^B_t = eps_t . delta(A),  eps_t = (-2/t)^{(t+3)/2}."""
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
    eps_t = jacobi((-2) % t, t) ** ((t + 3) // 2)
    return int(eps_t * s)


def nu_de(Lam, t, r):
    """nu(Lambda, .) indexado por x = 2(mu + rho_{D_r}), via GKRS."""
    mp = (t - 1) // 2
    Rp = mp + r
    v = [2 * Lam[i] + 2 * (Rp - i) - 1 for i in range(Rp)]
    out = {}
    for perm in itertools.permutations(range(Rp)):
        s = sgn_perm(list(perm))
        base = [v[perm[i]] for i in range(Rp)]
        for eps in itertools.product((1, -1), repeat=Rp):
            u = [base[i] * eps[i] for i in range(Rp)]
            if not (all(u[i] > u[i + 1] for i in range(mp - 1)) and u[mp - 1] > 0):
                continue
            f = u[mp:]
            if not (all(f[i] > f[i + 1] for i in range(r - 1)) and f[r - 2] > abs(f[r - 1])):
                continue
            sg = s
            for e in eps:
                sg *= e
            tv = tauB([u[i] - (2 * (mp - i) - 1) for i in range(mp)], t, mp)
            # ojo: tauB espera A = 2(eta+rho); aqui u[:mp] YA es ese vector
            tv = tauB(u[:mp], t, mp)
            if tv == 0:
                continue
            x = tuple(u[mp:])
            out[x] = out.get(x, 0) + sg * tv
    return {k: v2 for k, v2 in out.items() if v2 != 0}


def cabeza(d):
    return max(d, key=lambda k: (sum(k), k))


def dividir(nu, t, r, tope=20000):
    P = dict(nu)
    c = {}
    for _ in range(tope):
        P = {k: v for k, v in P.items() if v != 0}
        if not P:
            return c, {}
        y = cabeza(P)
        cand = None
        for eps in itertools.product((1, -1), repeat=r):
            e = enderezar(tuple(int(y[j]) - t * eps[j] for j in range(r)))
            if e is None:
                continue
            D = desplazar(e[0], t)
            if D and cabeza(D) == y:
                cand = (e[0], D)
                break
        if cand is None:
            return c, P
        x, D = cand
        if P[y] % D[y] != 0:
            return c, P
        cv = P[y] // D[y]
        c[x] = c.get(x, 0) + cv
        for k, v in D.items():
            nv = P.get(k, 0) - cv * v
            if nv == 0:
                P.pop(k, None)
            else:
                P[k] = nv
    return c, P


def dominantes(rk, cota):
    def rec(k, tope):
        if k == 0:
            yield ()
            return
        for a in range(tope, -1, -1):
            for resto in rec(k - 1, a):
                yield (a,) + resto
    return list(rec(rk, cota))


# ------------------------------------------------------------------ datos y controles
CASOS = [(3, 2, 6), (5, 2, 4), (7, 2, 3)]
TESTIGO = (3, 2, (6, 6, 6))

fallos = []
puntos = []
for (t, r, cota) in CASOS:
    mp, Rp = (t - 1) // 2, (t - 1) // 2 + r
    for Lam in dominantes(Rp, cota):
        nu = nu_de(Lam, t, r)
        if not nu:
            continue
        c, resto = dividir(nu, t, r)
        if resto:
            fallos.append(("resto no nulo", t, r, Lam))
            continue
        # C1
        lhs = {}
        for x, cv in c.items():
            for y, sg in desplazar(x, t).items():
                lhs[y] = lhs.get(y, 0) + cv * sg
        lhs = {k: v for k, v in lhs.items() if v != 0}
        if lhs != nu:
            fallos.append(("C1 c.Delta != nu", t, r, Lam))
        # C2
        if any(abs(v) != 1 for v in c.values()) or any(abs(v) != 1 for v in nu.values()):
            fallos.append(("C2 coeficiente no unitario", t, r, Lam))
        puntos.append((t, r, len(nu), len(c)))

# C3: cruce contra los numeros de SAGE
t0, r0, Lam0 = TESTIGO
nu0 = nu_de(Lam0, t0, r0)
c0, resto0 = dividir(nu0, t0, r0)
rhoD2 = [2 * (r0 - j) - 2 for j in range(r0)]
c0_mu = {tuple((x[j] - rhoD2[j]) // 2 for j in range(r0)): v for x, v in c0.items()}
ruta = os.path.join("..", "gates", "gkrs_L1_DUMP.json")
sage_c = None
if os.path.exists(ruta):
    for bloque in json.load(open(ruta)):
        for col in bloque.get("colisiones", []):
            if bloque["t"] == t0 and bloque["r"] == r0 and tuple(col["Lambda"]) == Lam0:
                sage_c = {eval(k): v for k, v in col["c"].items()}
if sage_c is None:
    fallos.append(("C3 no se encontro el testigo en el volcado de Sage", t0, r0, Lam0))
else:
    a = {k: abs(v) * (1 if v > 0 else -1) for k, v in c0_mu.items()}
    if a != sage_c and {k: -v for k, v in a.items()} != sage_c:
        fallos.append(("C3 discrepa con Sage", t0, r0, Lam0))

print("  control C3: testigo %s  |supp c| python = %d, sage = %s"
      % (str(Lam0), len(c0_mu), len(sage_c) if sage_c else "n/a"))
if fallos:
    print("  LA FIGURA NO SE DIBUJA. Controles fallidos:")
    for f in fallos[:8]:
        print("   ", f)
    sys.exit(1)
print("  controles: %d formas, C1 C2 C3 ok" % len(puntos))

# ------------------------------------------------------------------ dibujo
fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3))

# IZQUIERDA: una forma entera
ax = axes[0]
for x, v in c0.items():
    mu = tuple((x[j] - rhoD2[j]) // 2 for j in range(r0))
    ax.scatter([mu[0]], [mu[1]], s=46, marker="s", zorder=3,
               color=(PLUS if v > 0 else MINUS), edgecolors="none")
for x, v in nu0.items():
    ax.scatter([x[0] / 2.0], [x[1] / 2.0], s=150, marker="o", zorder=4,
               facecolors="none", edgecolors=(PLUS if v > 0 else MINUS), linewidths=1.6)
# las 2^r flechas desde un punto de c
x_ref = sorted(c0, key=lambda k: (sum(k), k))[-1]
mu_ref = tuple((x_ref[j] - rhoD2[j]) // 2 for j in range(r0))
for eps in itertools.product((1, -1), repeat=r0):
    y = tuple(x_ref[j] + t0 * eps[j] for j in range(r0))
    e = enderezar(y)
    if e is None:
        continue
    ax.annotate("", xy=(e[0][0] / 2.0, e[0][1] / 2.0), xytext=(mu_ref[0], mu_ref[1]),
                arrowprops=dict(arrowstyle="->", color=MUTED, linewidth=1.0, alpha=0.75),
                zorder=2)
ax.set_xlabel(r"$\mu_1$")
ax.set_ylabel(r"$\mu_2$")
ax.set_title(r"$t=3$, $r=2$, $\Lambda=(6,6,6)$: $\nu$ is the finite difference of $c$",
             fontsize=9.5, color=INK)
ax.text(0.5, -0.26,
        r"filled: $\mathrm{supp}\,c$ ($%d$ points).   ringed: $\mathrm{supp}\,\nu$ ($%d$)."
        r"   blue $+1$, red $-1$.   arrows: the $2^r$ shifts of one point"
        % (len(c0), len(nu0)),
        transform=ax.transAxes, ha="center", fontsize=8, color=SECOND)
ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

# DERECHA: el tamano de nu contra el de c, sobre toda la poblacion
ax = axes[1]
marcas = {3: "o", 5: "s", 7: "^"}
for t in sorted(set(p[0] for p in puntos)):
    xs = [p[2] for p in puntos if p[0] == t]
    ys = [p[3] for p in puntos if p[0] == t]
    ax.scatter(xs, ys, s=30, marker=marcas.get(t, "o"), alpha=0.75,
               facecolors="none", edgecolors=PLUS if t == 3 else (MINUS if t == 5 else SECOND),
               linewidths=1.1, label=r"$t=%d$" % t)
lim = max(max(p[3] for p in puntos), max(p[2] for p in puntos)) + 2
ax.plot([0, lim], [0, lim], color=MUTED, linewidth=0.9, linestyle="--", zorder=1)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
ax.set_xlabel(r"$|\mathrm{supp}\,\nu|$   (numerator: proved $\{0,\pm1\}$)")
ax.set_ylabel(r"$|\mathrm{supp}\,c|$   (quotient)")
ax.set_title("the quotient is far larger than the numerator, and still a unit",
             fontsize=9.5, color=INK)
ax.text(0.5, -0.26,
        r"every coefficient on both axes is $\pm1$, on all $%d$ shapes; the diagonal is $y=x$"
        % len(puntos),
        transform=ax.transAxes, ha="center", fontsize=8, color=SECOND)
ax.legend(loc="upper left", frameon=False, fontsize=8)
ax.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

for a in axes:
    for s in a.spines.values():
        s.set_color(MUTED)
        s.set_linewidth(0.6)

fig.tight_layout()
fig.subplots_adjust(bottom=0.20)
fig.savefig("fig_division.pdf")
print("  fig_division.pdf escrito")
print("  izquierda: |supp c| = %d, |supp nu| = %d" % (len(c0), len(nu0)))
print("  derecha  : %d formas; |supp nu| en [%d, %d], |supp c| en [%d, %d]"
      % (len(puntos), min(p[2] for p in puntos), max(p[2] for p in puntos),
         min(p[3] for p in puntos), max(p[3] for p in puntos)))
