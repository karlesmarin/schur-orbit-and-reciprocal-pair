# -*- coding: utf-8 -*-
# LA ESCALERA CON SOPORTE, EN TODAS LAS (t, r).  13 de agosto de 2026.
#
# POR QUE.  k_vs_m.py midio K(W) = W/2 - 3 en t=4 r=2 y de ahi salio "no hay ninguna cota de
# profundidad, se cae el instrumento entero".  witness_family.py enseno que esa profundidad cuenta
# peldanos de la escalera COMPLETA (D, D-2, D-4, ... se cuenten o no), y que la mayor parte de ellos
# NO TIENEN NI UN MONOMIO.  Separadas las dos cosas, en t=4 r=2:
#
#     W                 18  22  26  30  34  38  42
#     K  completa        6   8  10  12  14  16  18     crece, la recta
#     K_sop  cancelan    2   2   3   3   4   4   4     SE PARA EN 4
#
# y las dos familias explicitas de contact_order.py llegan a anchura 114 con prof = 54 y siguen
# cancelando 4.  O sea el enunciado "alguno de los K primeros estratos es no nulo" es falso para todo
# K en la escalera completa y parece CIERTO CON K = 5 en la escalera con soporte.
#
# ESTE GUION ES EL QUE PUEDE MATARLO.  Todo lo anterior es t=4 r=2.  Aqui se barre el SWEEP entero
# -- las siete configuraciones de k_vs_m.py -- imprimiendo K y K_sop lado a lado, y se extiende
# t=4 r=2 dos anchuras mas (46, 50).  Si K_sop pasa de 4 en cualquier sitio, el 4 es otro artefacto
# de rango y esta columna lo dice; si se para en 4 en las siete, hay un enunciado que intentar probar.
#
# DEFINICION, sin ambiguedad.  Para un superviviente, con D el grado maximo y 'first' el primer grado
# no nulo:
#     prof      = D - first                              peldanos de la escalera completa
#     cancelan  = #{ d en D, D-2, ..., first+2 : d TIENE bloque de Laplace }   los que se anulan de
#                                                                             verdad, con monomios
# y K, K_sop son el maximo de cada uno sobre los supervivientes de esa anchura.  cancelan >= 1
# siempre, porque el estrato de arriba se anula por definicion de superviviente.
#
# EL CRITERIO NO SE REESCRIBE: probe() sale de ejecutar el preambulo de k_vs_m.py, mismos bytes.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python support_ladder.py

import itertools
import json
import os
import sys
import time
from collections import Counter, defaultdict

from survivors_wide import scan

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
assert "def probe(" in _head and "def shapes_of_width(" in _head, "k_vs_m.py cambio de forma"
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]
shapes_of_width = _ns["shapes_of_width"]

# el SWEEP de k_vs_m.py, con t=4 r=2 extendido a 50 (era 36) porque es donde vive la saturacion
SWEEP = [(4, 2, 50), (6, 2, 26), (8, 2, 24), (6, 3, 24), (8, 3, 24), (6, 4, 23), (10, 2, 23)]
OUT_JSON = "support_ladder_SWEEP.json"

# testigos ya publicados (witness_family_OUT.txt / contact_order_OUT.txt): prof y cancelan
TESTIGOS = [
    ((18, 17, 11, 8, 7, 6, 1, 0), 6, 2),
    ((26, 17, 15, 12, 11, 10, 9, 0), 10, 3),
    ((38, 23, 21, 18, 17, 16, 15, 0), 16, 4),
    ((114, 61, 59, 56, 55, 54, 53, 0), 54, 4),
]

# ===================================================================== C0 ========================
print("=" * 116)
print("C0  ACEPTACION -- fatal.  El criterio contra scan(), y los testigos ya publicados re-medidos")
print("=" * 116)
print("")
bad = 0
for (t, r, M) in [(4, 2, 15), (6, 2, 17), (6, 3, 18)]:
    n_ref, cont_ref, sv_ref = scan(t, r, M)
    mine, mb = Counter(), []
    for comb in itertools.combinations(range(M + 1), t + 2 * r):
        beta = tuple(sorted(comb, reverse=True))
        rec = probe(beta, t, r, deep=False)
        if rec is None:
            continue
        mine[(rec['e'] == t, rec['surv'])] += 1
        if rec['surv']:
            mb.append(beta)
    ok = (sum(mine.values()) == n_ref and mine == cont_ref
          and sorted(mb) == sorted(x['beta'] for x in sv_ref))
    bad += not ok
    print("  C0a  t=%d r=%d M=%d : objetivo %d/%d, betas %d/%d   %s"
          % (t, r, M, sum(mine.values()), n_ref, len(mb), len(sv_ref), "ok" if ok else "*** FALLA ***"))
    sys.stdout.flush()
print("")
for (beta, p, c) in TESTIGOS:
    rec = probe(beta, 4, 2)
    ok = rec is not None and rec['surv'] and rec['prof'] == p and rec['vac_cancelan'] == c
    bad += not ok
    print("  C0b  %-40s prof %s/%d, cancelan %s/%d   %s"
          % (str(beta), rec['prof'] if rec else "-", p,
             rec['vac_cancelan'] if rec else "-", c, "ok" if ok else "*** FALLA ***"))
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")

# ===================================================================== N1 ========================
print("")
print("=" * 116)
print("N1  LAS DOS ESCALERAS, CONFIGURACION A CONFIGURACION.  Cada linea se imprime al terminar su W")
print("=" * 116)
DATA = defaultdict(dict)
for (t, r, Wmax) in SWEEP:
    N = t + 2 * r
    print("")
    print("  t=%d r=%d  (N=%d, tope W=%d por RELOJ)" % (t, r, N, Wmax))
    print("     W    formas  objetivo  surv |    K   K_sop |  histograma de los que CANCELAN")
    print("  " + "-" * 104)
    Kacc = Kacc_sop = 0
    t0 = time.time()
    for W in range(N - 1, Wmax + 1):
        n_sh = n_tg = 0
        can = Counter()
        K = K_sop = 0
        ns = 0
        for beta in shapes_of_width(W, N):
            n_sh += 1
            rec = probe(beta, t, r)
            if rec is None:
                continue
            n_tg += 1
            if not rec['surv'] or rec['prof'] is None:
                continue
            ns += 1
            K = max(K, rec['prof'])
            K_sop = max(K_sop, rec['vac_cancelan'])
            can[rec['vac_cancelan']] += 1
        nuevo = K_sop > Kacc_sop
        Kacc = max(Kacc, K)
        Kacc_sop = max(Kacc_sop, K_sop)
        DATA["%d_%d" % (t, r)][W] = dict(formas=n_sh, objetivo=n_tg, surv=ns, K=K, K_sop=K_sop,
                                         can={str(k): v for k, v in can.items()})
        print("  %4d %9d %9d %5d | %4d %6d  |  %-40s%s"
              % (W, n_sh, n_tg, ns, K, K_sop,
                 " ".join("c%d:%d" % (k, v) for k, v in sorted(can.items())) or "-",
                 "   <-- K_sop SUBE a %d" % K_sop if nuevo else ""))
        sys.stdout.flush()
        json.dump(dict(DATA), open(OUT_JSON, "w"), indent=1)
    print("     %d anchuras en %.0f s   |   K(<=Wmax) = %d,  K_sop(<=Wmax) = %d"
          % (Wmax - N + 2, time.time() - t0, Kacc, Kacc_sop))
    sys.stdout.flush()

# ===================================================================== N2 ========================
print("")
print("=" * 116)
print("N2  VEREDICTO -- K_sop acumulado por configuracion")
print("=" * 116)
print("")
print("     t   r  |  K max  |  K_sop max  |  anchura donde K_sop llega a su maximo")
print("  " + "-" * 84)
peor = 0
for (t, r, Wmax) in SWEEP:
    d = DATA["%d_%d" % (t, r)]
    if not d:
        continue
    Km = max(x['K'] for x in d.values())
    Ks = max(x['K_sop'] for x in d.values())
    donde = min(W for W in d if d[W]['K_sop'] == Ks)
    peor = max(peor, Ks)
    print("    %2d  %2d  | %6d  | %10d  |  W = %d  (tope del barrido %d)" % (t, r, Km, Ks, donde, Wmax))
print("")
print("  K_sop MAXIMO EN TODO EL BARRIDO : %d" % peor)
print("")
if peor <= 4:
    print("  Ninguna configuracion pasa de 4.  El enunciado que sobrevive es")
    print("      Phi_t != 0  =>  alguno de los %d PRIMEROS ESTRATOS CON SOPORTE es no nulo" % (peor + 1))
    print("  y la recta K = W/2 - 3 es, en sus tres cuartas partes, hueco de espectro.")
else:
    print("  ALGUNA configuracion pasa de 4: el 4 era de t=4 r=2 y no del objeto.  Se dice y se mira")
    print("  donde, que es lo que este guion existe para poder decir.")
print("")
print("  ALCANCE, sin adornarlo: los topes de W son de RELOJ.  Un 'no pasa de 4 hasta el tope' no es")
print("  'no pasa de 4': es un null con alcance, y el alcance es ese.  Lo que SI esta exhibido, y no")
print("  depende de ningun tope, son las dos familias de contact_order.py hasta anchura 114.")
print("")
print("DONE")
