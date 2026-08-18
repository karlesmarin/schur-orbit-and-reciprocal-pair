# -*- coding: utf-8 -*-
# LA RECTA K = W/2 - 3: ¿en que escalera?  Y QUIENES son los testigos.  13 de agosto de 2026.
#
# POR QUE.  k_vs_m.py midio, en t=4 r=2, cinco puntos sobre K = W/2 - 3 (W = 18, 22, 26, 30, 34), y de
# ahi salio el veredicto "no hay ninguna cota de profundidad, el instrumento entero se cae".  Antes de
# creerse eso hay que mirar QUE cuenta esa profundidad, porque hay DOS escaleras y no son la misma:
#
#   LA ESCALERA COMPLETA.  Los estratos son D, D-2, D-4, ..., se cuenten o no.  prof = D - first los
#   cuenta TODOS.  Es el convenio de depth.py:205 y el de los enunciados.
#   LA ESCALERA CON SOPORTE.  Los grados que de verdad tienen algun bloque de Laplace.  Un grado sin
#   NI UN monomio esta vacio por contabilidad, no porque algo cancele.
#
# La diferencia no es cosmetica y la propia cabecera de measure() en depth_histogram.py ya la habia
# senalado -- "no son lo mismo para quien quiera probar una cota" -- y luego la tabla de k_vs_m.py
# imprime prof y NO imprime el desglose, que si calcula (vac_cancelan / vac_sin_soporte) y archiva.
# Este guion imprime la columna que falta.  Si la profundidad CON SOPORTE se queda acotada mientras la
# completa crece como W/2 - 3, entonces la recta es una medida de lo HUECA que se vuelve la escalera
# al ensancharse beta, y la cota de profundidad no esta muerta: estaba enunciada en la escalera
# equivocada.  Si tambien crece, la cota esta muerta de verdad y esta columna lo dice igual de claro.
#
# EL TECHO.  Los grados de un T son 2*sum_S(T) - sum(T) sobre los C(2r,r) bloques, o sea van de +deg(T)
# a -deg(T).  El techo de prof es 2*D, y D crece con W.  Se mide tambien, para saber si K esta pegado
# al techo o muy por debajo.
#
# N2 son los TESTIGOS, uno a uno, con sus betas: si los de profundidad maxima de cada anchura son UNA
# FAMILIA (el mismo patron desplazado), la recta deja de ser un ajuste y pasa a ser un enunciado con
# testigo explicito, que es lo que se puede intentar demostrar.
#
# N3 extiende DOS anchuras por encima de donde llego k_vs_m.py (W = 38, 42): la ley se escribio viendo
# 18/22/26 y acerto 30/34; que acierte 38 y 42 es prediccion, no ajuste.
#
# EL CRITERIO NO SE REESCRIBE.  probe() y shapes_of_width() se toman de k_vs_m.py ejecutando su
# PREAMBULO (todo lo anterior al bloque C0): k_vs_m.py es un script recto, sin guard __main__, asi que
# no se puede importar -- pero si ejecutar su cabecera.  Son los MISMOS BYTES, no una copia.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python witness_family.py

import itertools
import json
import os
import sys
import time
from collections import Counter

from second_stratum import setup, all_transversals
from depth_histogram import stratify
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

T_, R_ = 4, 2
N_ = T_ + 2 * R_
WIDTHS = [18, 22, 26, 30, 34, 38, 42]
OUT_JSON = "witness_family_WITNESS.json"


def window_of(beta, t, r):
    """(D, min grado con soporte, num de grados CON SOPORTE, num de peldanos de la escalera completa)."""
    cl, E, Cd = setup(beta, t)
    tr = all_transversals(beta, cl, r, t)
    B = stratify([(x[2], x[1]) for x in tr], r)
    D = max(x[3] for x in tr)
    return D, min(B), len(B), D + 1


# ===================================================================== C0 ========================
print("=" * 108)
print("C0  ACEPTACION -- el criterio es el de k_vs_m.py, byte a byte, y se re-firma contra scan()")
print("=" * 108)
print("")
bad = 0
for (t, r, M) in [(4, 2, 15), (6, 3, 18)]:
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
print("       C0a %s" % ("PASA" if not bad else "FALLA"))
print("")

# C0b  LA SIMETRIA de la ventana es fatal (min = -D): es lo que hace que el techo sea 2*D.
#      LA DENSIDAD no se exige -- SE MIDE.  Yo la habia dado por supuesta y la primera corrida de este
#      guion me tumbo la suposicion en el primer beta: (20,19,18,5,3,2,1,0) tiene D=36, o sea 37
#      peldanos, y solo 29 con soporte.  La escalera es HUECA, y eso es justo el asunto del guion.
print("  C0b  la ventana es simetrica [-D, D] (fatal), y la DENSIDAD de la escalera se mide")
seen = 0
huecos = []
for beta in shapes_of_width(20, N_):
    rec = probe(beta, T_, R_, deep=False)
    if rec is None:
        continue
    D, mn, ns, nl = window_of(beta, T_, R_)
    if mn != -D:
        bad += 1
        print("       *** FALLA *** beta=%s  D=%d  min=%d  (la ventana NO es simetrica)"
              % (str(beta), D, mn))
        break
    huecos.append((ns, nl))
    seen += 1
    if seen >= 40:
        break
if huecos:
    dens = [a / float(b) for a, b in huecos]
    print("       simetria: ok en %d betas objetivo (W=20)" % seen)
    print("       densidad: %.1f%% de los peldanos tienen soporte de media (min %.1f%%, max %.1f%%)"
          % (100 * sum(dens) / len(dens), 100 * min(dens), 100 * max(dens)))
    print("       o sea la escalera completa NO es la escalera con soporte, y medido esta.")
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")

# ===================================================================== N1, N2 ====================
print("")
print("=" * 108)
print("N1  LAS DOS ESCALERAS   |   N2  los TESTIGOS de profundidad maxima, uno a uno")
print("=" * 108)
print("")
print("  prof      = D - first, peldanos de la escalera COMPLETA (el convenio de k_vs_m.py)")
print("  prof_sop  = de esos peldanos, los que TIENEN soporte y aun asi se anulan (los que cancelan)")
print("  K, K_sop  = el maximo de cada uno sobre los supervivientes de esa anchura")
print("")
DATA = {}
for W in WIDTHS:
    t0 = time.time()
    best = best_sop = -1
    wit = []
    hist = Counter()
    hist_sop = Counter()
    Dmax_pop = 0
    n_sh = n_tg = n_sv = 0
    for beta in shapes_of_width(W, N_):
        n_sh += 1
        rec = probe(beta, T_, R_)
        if rec is None:
            continue
        n_tg += 1
        if not rec['surv'] or rec['prof'] is None:
            continue
        n_sv += 1
        hist[rec['prof']] += 1
        hist_sop[rec['vac_cancelan']] += 1
        best_sop = max(best_sop, rec['vac_cancelan'])
        D = window_of(beta, T_, R_)[0]
        Dmax_pop = max(Dmax_pop, D)
        if rec['prof'] > best:
            best, wit = rec['prof'], []
        if rec['prof'] == best:
            wit.append((beta, D, rec['prof'], rec['vac_cancelan'], rec['vac_sin_soporte']))
    ley = W // 2 - 3
    print("  W=%d   formas %d, objetivo %d, supervivientes %d, %.0f s"
          % (W, n_sh, n_tg, n_sv, time.time() - t0))
    print("     K (escalera completa)  : %-4d   ley K = W/2 - 3 predice %-4d  %s"
          % (best, ley, "ACIERTA" if best == ley else "*** FALLA LA LEY ***"))
    print("     K_sop (con soporte)    : %-4d   <-- LA COLUMNA QUE FALTABA" % best_sop)
    print("     techo trivial 2*D      : %-4d   (D maximo de la poblacion objetivo = %d)"
          % (2 * Dmax_pop, Dmax_pop))
    print("     histograma prof        : %s"
          % " ".join("p%d:%d" % (k, v) for k, v in sorted(hist.items())))
    print("     histograma prof_sop    : %s"
          % " ".join("s%d:%d" % (k, v) for k, v in sorted(hist_sop.items())))
    print("     testigos de prof %d  (%d de ellos):" % (best, len(wit)))
    for (beta, D, p, vc, vs) in wit[:24]:
        print("        beta = %-40s D=%-4d prof=%-3d  de los %d peldanos vacios: %d CANCELAN, %d sin soporte"
              % (str(beta), D, p, p // 2, vc, vs))
    if len(wit) > 24:
        print("        ... y %d mas" % (len(wit) - 24))
    print("")
    sys.stdout.flush()
    DATA[W] = dict(K=best, K_sop=best_sop, ley=ley, Dmax=Dmax_pop,
                   n_sh=n_sh, n_tg=n_tg, n_sv=n_sv,
                   hist={str(k): v for k, v in hist.items()},
                   hist_sop={str(k): v for k, v in hist_sop.items()},
                   testigos=[[list(b), D, p, vc, vs] for (b, D, p, vc, vs) in wit])
    json.dump(DATA, open(OUT_JSON, "w"), indent=1)

# ===================================================================== N3 ========================
print("")
print("=" * 108)
print("N3  VEREDICTO -- las dos escaleras, una al lado de la otra")
print("=" * 108)
print("")
print("     W    | K completa | ley W/2-3 | K_sop | techo 2D | K/techo")
print("  " + "-" * 72)
for W in WIDTHS:
    d = DATA[W]
    print("  %4d   | %10d | %9d | %5d | %8d | %.3f"
          % (W, d['K'], d['ley'], d['K_sop'], 2 * d['Dmax'],
             d['K'] / (2.0 * d['Dmax']) if d['Dmax'] else 0))
print("")
fallos = [W for W in WIDTHS if DATA[W]['K'] != DATA[W]['ley']]
print("  la ley K = W/2 - 3 %s"
      % ("ACIERTA LAS %d ANCHURAS (38 y 42 son PREDICCION, no ajuste)" % len(WIDTHS) if not fallos
         else "FALLA en W = %s" % fallos))
sop = [DATA[W]['K_sop'] for W in WIDTHS]
print("  K_sop = %s" % sop)
if max(sop) == min(sop):
    print("  K_sop es CONSTANTE en todo el rango: la recta mide lo HUECA que se vuelve la escalera,")
    print("  no una profundidad de cancelacion.  La cota estaba enunciada en la escalera equivocada.")
elif max(sop) - min(sop) < (max(DATA[W]['K'] for W in WIDTHS) - min(DATA[W]['K'] for W in WIDTHS)):
    print("  K_sop crece MUCHO MAS DESPACIO que K: la recta es en su mayor parte hueco de espectro.")
else:
    print("  K_sop crece igual que K: la cancelacion es real y la cota esta muerta de verdad.")
print("")
print("DONE")
