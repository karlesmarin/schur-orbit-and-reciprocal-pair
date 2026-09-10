# -*- coding: utf-8 -*-
# FIGURA DEL VALOR CON SIGNO.   8 de septiembre de 2026.
# El resultado central del articulo --- el valor exacto con su signo, y su factorizacion --- no
# tenia ninguna figura.  Esta la da: los pesos supervivientes de Sp(6) en k=2 (d=2, q=4), en el
# espacio de (lambda_1, lambda_2, lambda_3).
#   izquierda:  color = SIGNO del caracter.  Se ve que el signo no es una funcion de la forma:
#               alterna en laminas, que es lo que dice F_q sobre los residuos ordenados.
#   derecha:    tamano y color = |chi|, que es el producto de dimensiones del Teorema de las dims.
# Los numeros del pie se calculan aqui y se imprimen al correr; no se escriben a mano.
# Run:  python fig_value.py        (desde gates/; escribe en ../note_prasad2/)
import sys
from math import gcd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

AZUL = "#1F4E79"
NARANJA = "#C55A11"


def dominants(m, top):
    def rec(k, hi):
        if k == 0:
            yield ()
            return
        for v in range(hi, -1, -1):
            for rest in rec(k - 1, v):
                yield (v,) + rest
    return rec(m, top)


def h_closed(q, d, N):
    base = [0] * (N + 1)
    k = 0
    while q * k <= N:
        c = 1
        for i in range(1, d):
            c = c * (k + i) // i
        base[q * k] = c
        k += 1
    return [base[n] - (2 * base[n - 1] if n >= 1 else 0) + (base[n - 2] if n >= 2 else 0)
            if d % 2 == 0 else base[n] - (base[n - 2] if n >= 2 else 0) for n in range(N + 1)]


def det_int(M):
    n = len(M)
    A = [r[:] for r in M]
    sg, prev = 1, 1
    for c in range(n):
        if A[c][c] == 0:
            p = next((r for r in range(c + 1, n) if A[r][c] != 0), None)
            if p is None:
                return 0
            A[c], A[p] = A[p], A[c]; sg = -sg
        for r in range(c + 1, n):
            for k in range(c + 1, n):
                A[r][k] = (A[r][k] * A[c][c] - A[r][c] * A[c][k]) // prev
            A[r][c] = 0
        prev = A[c][c]
    return sg * A[n - 1][n - 1]


def chi(lam, m, h):
    H = lambda n: (h[n] if 0 <= n < len(h) else 0)
    M = [[H(lam[i] - (i + 1) + (j + 1)) + H(lam[i] - (i + 1) - (j + 1) + 2)
          for j in range(m)] for i in range(m)]
    return det_int(M) // 2


m, k, TOP = 3, 2, 13
t = 2 * m + 2
d = gcd(k, t)
q = t // d
h = h_closed(q, d, 3 * (TOP + m) + 6)

pos, neg, vals = [], [], []
for lam in dominants(m, TOP):
    v = chi(lam, m, h)
    if v == 0:
        continue
    (pos if v > 0 else neg).append(lam)
    vals.append((lam, abs(v)))

vmax = max(v for _, v in vals)
fig = plt.figure(figsize=(11.2, 4.8))

ax = fig.add_subplot(1, 2, 1, projection="3d")
for datos, col, etq in ((pos, AZUL, r"$\chi>0$"), (neg, NARANJA, r"$\chi<0$")):
    ax.scatter([p[0] for p in datos], [p[1] for p in datos], [p[2] for p in datos],
               s=16, color=col, depthshade=False, label=etq)
ax.set_title(r"the sign of $\chi_\lambda(a_P^2)$ in $\mathrm{Sp}(6)$", fontsize=10)
ax.legend(fontsize=8, loc="upper left")

ax2 = fig.add_subplot(1, 2, 2, projection="3d")
sc = ax2.scatter([l[0] for l, _ in vals], [l[1] for l, _ in vals], [l[2] for l, _ in vals],
                 c=[v for _, v in vals], s=[8 + 60 * (v / vmax) ** 0.5 for _, v in vals],
                 cmap="viridis", depthshade=False)
cb = fig.colorbar(sc, ax=ax2, fraction=0.03, pad=0.02)
cb.ax.tick_params(labelsize=7)
cb.set_label(r"$|\chi_\lambda(a_P^2)|$", fontsize=8)
ax2.set_title(r"its modulus, $\kappa\,D_0\,D_{q/2}\prod_c D_c$", fontsize=10)

for a in (ax, ax2):
    a.set_xlabel(r"$\lambda_1$", fontsize=8)
    a.set_ylabel(r"$\lambda_2$", fontsize=8)
    a.set_zlabel(r"$\lambda_3$", fontsize=8)
    a.tick_params(labelsize=6)
    a.view_init(elev=20, azim=35)

fig.tight_layout()
fig.savefig("../note_prasad2/fig_value.pdf", bbox_inches="tight")

print("Sp(2m) con m=%d, k=%d: t=%d, d=%d, q=%d, lambda_1 <= %d" % (m, k, t, d, q, TOP))
print("supervivientes: %d   de los cuales positivos %d y negativos %d"
      % (len(vals), len(pos), len(neg)))
print("|chi| maximo: %d" % vmax)
print("ESCRITO ../note_prasad2/fig_value.pdf")
