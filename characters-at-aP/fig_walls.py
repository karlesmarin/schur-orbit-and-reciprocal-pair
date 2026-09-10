# -*- coding: utf-8 -*-
# FIGURA: el criterio como un punto que esquiva paredes modulo q.   8 de septiembre de 2026.
# Izquierda, C_2 en el plano: los q^2 pesos modulo q, las paredes de raiz, y los que sobreviven.
# Derecha, C_3 en el espacio: lo mismo una dimension mas arriba, que es donde deja de poder
# dibujarse la intuicion de "clases" y solo queda el arreglo.
# Los numeros del pie NO se escriben a mano: se calculan aqui y se imprimen al correr.
# Run:  python fig_walls.py       (desde gates/; escribe en ../note_prasad2/)
import itertools
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

AZUL = "#1F4E79"
NARANJA = "#C55A11"
GRIS = "#B8C4D0"

RAICES = json.load(open("root_data.json", encoding="utf-8"))
POR_TIPO = {d["tipo"]: d for d in RAICES}


def survivors(tipo, q):
    """devuelve r_q, los supervivientes, y los demas con su N_q (cuantas paredes tocan)."""
    d = POR_TIPO[tipo]
    n, R = d["rango"], d["raices"]
    r_q = sum(1 for c, v in R if (c * sum(v)) % q == 0)
    vivos, muertos = [], []
    for x in itertools.product(range(q), repeat=n):
        k = sum(1 for c, v in R
                if (c * sum(x[j] * v[j] for j in range(n))) % q == 0)
        (vivos if k == r_q else muertos).append((x, k))
    return r_q, [p for p, _ in vivos], muertos, d["exponentes"]


def prod(v):
    r = 1
    for x in v:
        r *= x
    return r


fig = plt.figure(figsize=(11.6, 4.6))

# ---------------------------------------------------------------- C2, dos valores de q
info = {}
for idx, q in enumerate((6, 7)):
    ax = fig.add_subplot(1, 3, idx + 1)
    r_q, vivos, muertos, exps = survivors("C2", q)
    info[q] = (r_q, len(vivos), prod([q - m for m in exps]))
    kmax = max(k for _, k in muertos)
    sc = ax.scatter([p[0] for p, _ in muertos], [p[1] for p, _ in muertos],
                    c=[k for _, k in muertos], s=[18 + 26 * k / kmax for _, k in muertos],
                    cmap="OrRd", vmin=0, vmax=kmax, zorder=2, edgecolors="none")
    ax.scatter([p[0] for p in vivos], [p[1] for p in vivos], s=62, color=AZUL, zorder=3)
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.03, ticks=range(1, kmax + 1))
    cb.ax.tick_params(labelsize=6)
    cb.set_label(r"$N_q(x)$", fontsize=7)
    ax.set_xticks(range(q)); ax.set_yticks(range(q))
    ax.set_xlim(-0.7, q - 0.3); ax.set_ylim(-0.7, q - 0.3)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=7)
    ax.set_title(r"$C_2$, $q=%d$   ($r_q=%d$)" % (q, r_q), fontsize=10)
    ax.set_xlabel(r"$c_1$", fontsize=9); ax.set_ylabel(r"$c_2$", fontsize=9)
    for s in ax.spines.values():
        s.set_color("#888888")

# ---------------------------------------------------------------- C3 en 3D
ax = fig.add_subplot(1, 3, 3, projection="3d")
Q3 = 7
r3, vivos3, muertos3, exps3 = survivors("C3", Q3)
ax.scatter([p[0] for p, _ in muertos3], [p[1] for p, _ in muertos3],
           [p[2] for p, _ in muertos3], s=7, color=GRIS, alpha=0.30, depthshade=False)
ax.scatter([p[0] for p in vivos3], [p[1] for p in vivos3], [p[2] for p in vivos3],
           s=26, color=AZUL, depthshade=False)
ax.set_title(r"$C_3$, $q=%d$   ($r_q=%d$)" % (Q3, r3), fontsize=10)
ax.set_xlabel(r"$c_1$", fontsize=8); ax.set_ylabel(r"$c_2$", fontsize=8)
ax.set_zlabel(r"$c_3$", fontsize=8)
ax.tick_params(labelsize=6)
ax.view_init(elev=20, azim=32)

fig.tight_layout()
fig.savefig("../note_prasad2/fig_walls.pdf", bbox_inches="tight")

print("C2 q=6 :  r_q=%d  supervivientes=%d   prod(q-m_i)=%d" % info[6])
print("C2 q=7 :  r_q=%d  supervivientes=%d   prod(q-m_i)=%d" % info[7])
print("C3 q=%d:  r_q=%d  supervivientes=%d   prod(q-m_i)=%d"
      % (Q3, r3, len(vivos3), prod([Q3 - m for m in exps3])))
print("ESCRITO ../note_prasad2/fig_walls.pdf")
