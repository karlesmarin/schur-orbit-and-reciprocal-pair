# -*- coding: utf-8 -*-
# ¿ES EL t-CORE EL DISCRIMINANTE?   Condicionar, que es lo que convierte una señal en criterio.
# 15 de agosto de 2026.
#
# DE DONDE SALE.  twelve_forms.sage encontro que las anulaciones residuales se concentran en TRES
# t-cores por configuracion (t=4: (5,2,1), (), (6,3,2,1);  t=6: (1,1,1), (2,1), (2,2,2,1)).
# Eso NO es todavia un criterio: falta el reciproco.  La pregunta correcta no es "que core tienen
# las nulas" sino "de las formas con ESE core, cuantas se anulan".
#
#   * si P(nula | core = c) es 1 para esos c y 0 para el resto  -> EL CORE ES EL CRITERIO.
#   * si es baja, el core solo es una condicion necesaria mas, como la paridad, y se dice.
#
# Es exactamente la leccion de [[an-aggregate-count-is-not-a-case]]: el conteo agregado no dice nada
# hasta que se condiciona.
#
# CONTROLES
#   C0  se imprime la tabla ENTERA de cores, no solo los de las nulas: si un core con 40 formas tiene
#       0 nulas, eso es tan informativo como el que las tiene todas.
#   C1  SEÑUELO: se repite el condicionamiento sobre un invariante que NO deberia decidir -- el
#       tamaño |lambda| -- para ver cuanta concentracion produce el azar en esta poblacion.  Si el
#       core no concentra mas que |lambda|, no hay señal.
#   C2  no vacuidad: n por core impreso siempre, y los cores con una sola forma se marcan, porque
#       P=1 sobre n=1 no es nada.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  MSYS_NO_PATHCONV=1 docker run --rm -v "E:/proyectos/Curiosity/research/orbit-pair/gates:/work" \
#         -w /work sage-normaliz:local sage core_conditional.sage

import itertools, json, sys
from collections import defaultdict

load("pob_helper.py")


def core_de(lam, tt):
    n = len(lam)
    beta = [lam[i] + (n - 1 - i) for i in range(n)]
    clases = defaultdict(list)
    for b in beta:
        clases[b % tt].append(b)
    nuevo = []
    for q in range(tt):
        nuevo += [q + tt * i for i in range(len(clases[q]))]
    nuevo.sort(reverse=True)
    core = tuple(nuevo[i] - (len(nuevo) - 1 - i) for i in range(len(nuevo)))
    return tuple(x for x in core if x)


for (t, r, W) in [(4, 2, 13), (6, 2, 13)]:
    N = t + 2 * r
    d = list(range(N - 1, -1, -1))
    print("=" * 118)
    print("EL t-CORE COMO DISCRIMINANTE   --   t=%d  r=%d  W=%d  N=%d" % (t, r, W, N))
    print("=" * 118)
    print("")
    porcore = defaultdict(lambda: [0, 0])       # core -> [total ocupadas, nulas]
    portam = defaultdict(lambda: [0, 0])        # |lambda| -> idem  (señuelo C1)
    for b in betas_py(t, r, W):
        if not occupied_py(b, t):
            continue
        z = phi_zero_py(b, t, r)
        if z is None:
            continue
        lam = tuple(x for x in (b[i] - d[i] for i in range(N)) if x != 0)
        c = core_de(lam, t)
        porcore[c][0] += 1
        porcore[c][1] += bool(z)
        portam[sum(lam)][0] += 1
        portam[sum(lam)][1] += bool(z)

    tot = sum(v[0] for v in porcore.values())
    totn = sum(v[1] for v in porcore.values())
    print("   %d formas ocupadas, %d nulas, repartidas en %d cores distintos" % (tot, totn, len(porcore)))
    print("")
    print("   core                      |   n | nulas | P(nula|core) | ")
    print("   " + "-" * 84)
    filas = sorted(porcore.items(), key=lambda kv: (-kv[1][1], -kv[1][0]))
    for c, (n, nl) in filas:
        if nl == 0 and n < 12:
            continue                              # se resumen abajo
        marca = "  <- n=1, no dice nada" if n == 1 else ""
        print("   %-25s | %3d | %5d | %11.1f %% |%s" % (str(c), n, nl, 100.0 * nl / n, marca))
    resto_n = sum(n for c, (n, nl) in filas if nl == 0 and n < 12)
    resto_c = sum(1 for c, (n, nl) in filas if nl == 0 and n < 12)
    print("   %-25s | %3d | %5d | %11.1f %% |  (%d cores agrupados)"
          % ("... todos los demas", resto_n, 0, 0.0, resto_c))
    print("")
    con = [c for c, (n, nl) in porcore.items() if nl > 0]
    n_con = sum(porcore[c][0] for c in con)
    print("   RECIPROCO:  los %d cores que tienen alguna nula suman %d formas ocupadas, de las que"
          % (len(con), n_con))
    print("               %d se anulan  ->  P(nula | core con alguna nula) = %.1f %%"
          % (totn, 100.0 * totn / max(1, n_con)))
    if n_con == totn:
        print("               *** TODAS las formas de esos cores se anulan: EL CORE ES EL CRITERIO ***")
    else:
        print("               -> el core NO es suficiente: %d formas comparten core con una nula y NO"
              % (n_con - totn))
        print("                  se anulan.  Es una condicion necesaria mas, como la paridad.")
    print("")
    # ------------------------------------------------------------------ C1 señuelo --------------
    con_t = [k for k, (n, nl) in portam.items() if nl > 0]
    n_con_t = sum(portam[k][0] for k in con_t)
    print("   C1  SEÑUELO, el mismo condicionamiento sobre |lambda| (que no deberia decidir):")
    print("       %d valores de |lambda| tienen alguna nula, suman %d formas, P = %.1f %%"
          % (len(con_t), n_con_t, 100.0 * totn / max(1, n_con_t)))
    print("       Si el core no concentra MUCHO mas que esto, no hay señal en el core.")
    print("")
    sys.stdout.flush()

print("=" * 118)
print("DONE")
