# -*- coding: utf-8 -*-
# LA PREGUNTA QUE ABRE S4.  12 de agosto de 2026 (noche).
#
# two_strata_audit.py encontro que los contraejemplos NO son 8 ni son solo de t=6: la poblacion CRECE
# con M (4, 8, 16, 24 en t=6,r=3 para M=16..19) y aparece tambien en t=4 (20 en t=4,r=2,M=19).  Eso
# obliga a re-preguntar lo unico que quedaba en pie:
#
#     "Phi_t != 0  =>  alguno de los TRES PRIMEROS estratos es no nulo"     (la noche del 12,
#      medido en 100 288 formas, profundidad <= 4)
#
# Si algun contraejemplo del enunciado de dos estratos tiene su primer estrato no nulo POR DEBAJO de
# Dmax - 4, ese enunciado tambien muere y no queda ninguna ruta en pie.  Estas formas son la peor
# poblacion posible para el, o sea el test mas duro que se le puede hacer.
#
# COLUMNAS
#   P1  profundidad por arriba (Dmax - primer grado total no nulo) de CADA contraejemplo.
#       Prediccion de la noche del 12: <= 4.
#   P2  profundidad por abajo (primer grado |.| no nulo - Dmin).
#   P3  Phi_t != 0 en todos (si alguno fuera 0 no seria contraejemplo, seria algo peor).
#   P4  e = t en todos, o sea S = beta entero.  Es la firma que V2_RESUME apunta como semilla.
#   P5  cuantos son ESENCIALMENTE distintos: forma canonica bajo traslacion + complemento.
#   C0  CONTROL: la misma cuenta sobre formas que NO son contraejemplos debe dar profundidad 0
#       (el estrato de arriba sobrevive), o el medidor de profundidad esta mal.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python two_strata_depth.py

import itertools
import sys
from collections import defaultdict

sys.path.insert(0, ".")
from two_strata_audit import analyse, is_fail, full_expansion, lam_of   # noqa: E402

CONFIGS = [(4, 2, 19), (4, 3, 17), (6, 3, 19)]


def canon(beta):
    """forma canonica bajo traslacion (min -> 0) y complemento (M - beta)."""
    b = tuple(sorted((v - min(beta) for v in beta), reverse=True))
    c = tuple(sorted((max(b) - v for v in b), reverse=True))
    return min(b, c)


def depth_of(a, r):
    FE = full_expansion([(x[1], x[0]) for x in a['tr']], r)
    if not FE:
        return None
    mx = max(sum(k) for k in FE)
    mn = min(sum(abs(x) for x in k) for k in FE)
    return (a['Dmax'] - mx, mn - a['Dmin'], len(FE))


print("=" * 108)
print("P1-P5  la profundidad de TODOS los contraejemplos, con la expansion entera")
print("=" * 108)
print("")
print("     t   r    M  beta                                        |lam|  e  t  Dmax  prof.arriba"
      "  prof.abajo")
print("  " + "-" * 104)

allc = []
for (t, r, M) in CONFIGS:
    N = t + 2 * r
    for comb in itertools.combinations(range(M + 1), N):
        beta = sorted(comb, reverse=True)
        a = analyse(beta, t, r)
        if not is_fail(a):
            continue
        d = depth_of(a, r)
        allc.append((t, r, M, tuple(beta), a, d))
    sys.stdout.flush()

seen = set()
for (t, r, M, beta, a, d) in allc:
    c = canon(beta)
    mark = "" if c in seen else "  <- nueva"
    seen.add(c)
    lam = lam_of(beta)
    print("  %4d %3d %4d  %-42s %5d %2d %2d %5d %12s %11s%s"
          % (t, r, M, str(list(beta)), sum(lam), a['e'], t, a['Dmax'],
             "Phi=0!" if d is None else str(d[0]), "-" if d is None else str(d[1]), mark))

print("")
print("  P3  Phi_t != 0 en todos                 : %d de %d"
      % (sum(1 for x in allc if x[5] is not None), len(allc)))
print("  P4  e = t en todos                      : %d de %d"
      % (sum(1 for x in allc if x[4]['e'] == x[0]), len(allc)))
print("  P5  esencialmente distintos (canonicos) : %d de %d" % (len(seen), len(allc)))
prof = defaultdict(int)
for x in allc:
    if x[5] is not None:
        prof[x[5][0]] += 1
print("")
print("  P1  PROFUNDIDAD POR ARRIBA (Dmax - primer grado no nulo).  Prediccion de la noche: <= 4")
for k in sorted(prof):
    print("        %2d : %d formas%s" % (k, prof[k], "   *** POR DEBAJO DE 4 ***" if k > 4 else ""))
prof2 = defaultdict(int)
for x in allc:
    if x[5] is not None:
        prof2[x[5][1]] += 1
print("")
print("  P2  PROFUNDIDAD POR ABAJO (primer grado no nulo - Dmin)")
for k in sorted(prof2):
    print("        %2d : %d formas" % (k, prof2[k]))

# ---------------------------------------------------------------- C0 ----------------------------
print("")
print("=" * 108)
print("C0  CONTROL: en formas que NO son contraejemplos la profundidad por arriba debe ser 0")
print("=" * 108)
print("")
n = bad = 0
for comb in itertools.combinations(range(16), 8):
    beta = sorted(comb, reverse=True)
    a = analyse(beta, 4, 2)
    if a is None or a['crit'] or a['top_zero']:
        continue
    n += 1
    if n > 300:
        break
    d = depth_of(a, 2)
    if d is None or d[0] != 0:
        bad += 1
print("  %d formas con [Phi]_top != 0 : %d con profundidad != 0  (debe ser 0)" % (min(n, 300), bad))

print("")
if any(x[5] is not None and x[5][0] > 4 for x in allc):
    print("VEREDICTO: hay contraejemplos POR DEBAJO del tercer estrato -- la ruta de la noche del 12")
    print("           tambien cae.")
elif bad:
    print("VEREDICTO SUSPENDIDO: el control C0 fallo, el medidor de profundidad no es de fiar.")
else:
    print("VEREDICTO: todos los contraejemplos caen dentro de los TRES PRIMEROS estratos.")
    print("           La ruta de la noche del 12 sobrevive a su peor poblacion.")
print("DONE")
