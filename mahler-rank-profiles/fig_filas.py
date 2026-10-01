# -*- coding: utf-8 -*-
# Copyright (c) 2026 Carles Marin. All rights reserved.
# Carles Marin <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
#
# QUE: figura 1 del paper III. Estado MEDIDO (pivotes de la matriz de masas, mahler_masas) de todas las filas
#      0 < J < p^a del caracter cuadratico (-n/.), (p, n) = (5, 75619), piso a = 4 (625 columnas), dibujado en la
#      rejilla J + 1 = p^s u (p no divide u): u en horizontal, s en vertical. Phi(J) = p(J+1) - 1 sube s en 1 y
#      conserva u; la componente de Teichmuller es b = u - 1 mod (p-1), independiente de s. Teorema G: J ausente sii
#      u < p r_b, asi que cada columna es entera ausente o entera presente, con umbral p r_b por clase de u.
#      r_b por la serie de Iwasawa (filas_negativas.serie_iwasawa). Control: se comprueba medido == Teorema G antes
#      de dibujar. Vueltas del piso 3 (J_k = 24 + 20(k-1)) marcadas: J_k + 1 = 5((p-1)k+1), y las copias k = pj+1
#      caen en la columna u = (p-1)j+1 un piso mas arriba.
#      Salida: fig_filas.pdf, fig_filas.png (en /out), fig_filas_OUT.txt (datos y control).
# FECHA: 2-oct-2026
# Uso (docker, /work = ext_clasicos, /out = note_iwasawa): sage -python fig_filas.py
import json
import sys

sys.path.insert(0, '/work')
import filas_negativas as FN  # noqa: E402
import mahler_masas as MM  # noqa: E402
from sage.all import kronecker_character  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

p, n, a = 5, 75619, 4
P = p ** a
chi = kronecker_character(-n)
ex = [-1 if v == 0 else (0 if v == 1 else 1) for v in chi.values()]
mu, g, f = FN.medida(ex, 2, n, p, 4)
r = [FN.lam(FN.serie_iwasawa(mu, p, 4, b, 4 * p)) for b in range(p - 1)]
piv = set(MM.cuadratico(n, p, a, P - 1)[1])
aus = {J for J in range(1, P) if J not in piv}


def u_s(J):
    x, s = J + 1, 0
    while x % p == 0:
        x //= p
        s += 1
    return x, s


pred = {J for J in range(1, P) if r[(J % (p - 1))] is None or u_s(J)[0] < p * r[J % (p - 1)]}
ok = aus == pred
print("(p,n,a) = (%d,%d,%d); r_b (b = 0..%d) = %s; medido == Teorema G en las %d filas positivas: %s" % (
    p, n, a, p - 2, r, P - 1, ok), flush=True)
assert ok

UMAX = 60
fig, ax = plt.subplots(figsize=(7.0, 2.9))
col_aus = {0: '#1f4e79', 2: '#7a3b8f'}   # componentes pares
col_imp = '#b8b8b8'                      # componentes impares (r = infinito)
col_pres = '#f2f2f2'
for J in range(1, P):
    u, s = u_s(J)
    if u > UMAX:
        continue
    b = J % (p - 1)
    if J in aus:
        c = col_imp if b % 2 else col_aus[b]
    else:
        c = col_pres
    ax.add_patch(Rectangle((u - 0.45, s - 0.42), 0.9, 0.84, facecolor=c,
                           edgecolor=('#9a9a9a' if c == col_pres else 'white'), linewidth=0.5))
# umbrales p r_b de las componentes pares
for b in (0, 2):
    if r[b] is not None:
        x = p * r[b]
        ax.plot([x, x], [-0.6, a + 0.5], color=col_aus[b], linewidth=1.2, linestyle='--')
        ax.text(x + 0.4, a + 0.45, r'$p\,r_{%d}=%d$' % (b, x), color=col_aus[b], fontsize=8, va='top')
# vueltas del piso 3: J_k + 1 = 5((p-1)k + 1)
for k in range(1, 9):
    J = p ** (3 - 2) * ((p - 1) * k + 1) - 1
    u, s = u_s(J)
    if u <= UMAX and J < P:
        ax.text(u, s, str(k), ha='center', va='center', fontsize=7,
                color='white' if J in aus else 'black', fontweight='bold')
ax.set_xlim(0, UMAX + 1)
ax.set_ylim(-0.6, a + 0.6)
ax.set_yticks(range(a + 1))
ax.set_ylabel(r'$s=v_p(J+1)$')
ax.set_xlabel(r'$u=(J+1)/p^{s}$')
ax.set_xticks([1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
leg = [Rectangle((0, 0), 1, 1, facecolor=col_aus[0]), Rectangle((0, 0), 1, 1, facecolor=col_aus[2]),
       Rectangle((0, 0), 1, 1, facecolor=col_imp), Rectangle((0, 0), 1, 1, facecolor=col_pres, edgecolor='#999999')]
ax.legend(leg, [r'dependent, $b=0$', r'dependent, $b=2$', r'dependent, $b$ odd', 'new'], fontsize=7, ncol=4,
          loc='upper center', bbox_to_anchor=(0.5, 1.22), frameon=False)
fig.tight_layout()
fig.savefig('/out/fig_filas.pdf')
fig.savefig('/out/fig_filas.png', dpi=200)
with open('/work/fig_filas_datos.json', 'w', encoding='utf-8') as fh:
    json.dump({'p': p, 'n': n, 'a': a, 'r_b': [None if x is None else int(x) for x in r],
               'ausentes': sorted(int(J) for J in aus)}, fh)
print("figura escrita: /out/fig_filas.pdf y .png (u <= %d)" % UMAX, flush=True)
