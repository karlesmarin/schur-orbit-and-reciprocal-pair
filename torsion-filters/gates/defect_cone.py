# -*- coding: utf-8 -*-
# LA LEY DEL DEFECTO INVARIANTE, y el cono de generadores.  13 de agosto de 2026.
#
# DE DONDE SALE.  contact_order.py encontro que los testigos profundos de t=4 r=2 se generan sumando
# v = (4,2,2,2,2,2,2,0), que la familia aguanta 25 miembros hasta anchura 114, y que el defecto de
# concentricidad  d_i(beta) = beta_i + beta_{N-1-i} - C  vale [0,0,-1,-3] en TODOS los miembros.  Que
# sea constante no es suerte: sumale sus pares opuestos al generador y salen  4+0 = 2+2 = 2+2 = 2+2 = 4.
# EL GENERADOR ES EL MISMO CONCENTRICO.  Y entonces, para cualquier v con v_i + v_{N-1-i} = c:
#
#     d_i(beta + j v) = (beta_i + j v_i) + (beta_{N-1-i} + j v_{N-1-i}) - (C + j c) = d_i(beta)
#
#     LEY DEL DEFECTO INVARIANTE.  El defecto de concentricidad es invariante al sumar cualquier
#     vector concentrico.  Los v no-crecientes con v_i + v_{N-1-i} constante forman un CONO; cada
#     rayo genera una familia de defecto fijo y escala creciente, o sea de orden de contacto -> infinito
#     con el lugar donde Phi_t se anula.
#
# ESO ES UNA PREDICCION, no un resumen: si la ley es del DEFECTO y no del vector concreto, cualquier
# otro rayo del cono tiene que generar familia tambien.  Y hay uno especialmente limpio,
#
#     v* = (4,4,4,4,0,0,0,0)
#
# concentrico (4+0 en los cuatro pares) y ademas  ==  0 (mod 4), o sea que NO MUEVE LAS CLASES DE
# RESIDUOS -- al reves que v, que las desplaza 2 y tiene periodo 2.  Si v* genera familia, la ley es
# del defecto; si NO la genera, la ley es falsa como esta escrita y el generador original tenia algo
# mas dentro que la concentricidad.  Las dos salidas estan escritas ANTES de correr.
#
# QUE SE MIDE
#   N1  el CONO entero para N=8 con entradas <= 6: todos los v no-crecientes, concentricos, y se
#       itera cada uno desde la misma semilla.  Se cuenta cuantos generan familia.
#   N2  el SENUELO, y es el control que decide: los v NO concentricos, mismo tamano, misma semilla.
#       Si tambien generasen familias, la concentricidad de v no seria la causa y la ley sobraria.
#   N3  el defecto medido miembro a miembro, que es lo que la ley predice constante.
#
# EL CRITERIO NO SE REESCRIBE: probe() sale de ejecutar el preambulo de k_vs_m.py, mismos bytes.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python defect_cone.py

import itertools
import json
import os
import sys
from collections import Counter

from survivors_wide import scan

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
assert "def probe(" in _head, "k_vs_m.py cambio de forma"
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]

T_, R_ = 4, 2
N_ = 8
SEMILLA = (18, 17, 11, 8, 7, 6, 1, 0)
SEMILLA_B = (18, 13, 11, 8, 7, 6, 5, 0)
JMAX = 20
VMAX = 6
PREDICHO = (4, 4, 4, 4, 0, 0, 0, 0)


def defecto(beta):
    C = beta[0] + beta[-1]
    return tuple(beta[i] + beta[len(beta) - 1 - i] - C for i in range(len(beta) // 2))


def concentrico(v):
    s = {v[i] + v[len(v) - 1 - i] for i in range(len(v) // 2)}
    return len(s) == 1


def itera(beta0, v, jmax=JMAX):
    """[(j, W, prof, cancelan, defecto)] mientras siga siendo superviviente."""
    filas = []
    for j in range(jmax + 1):
        beta = tuple(b + j * x for b, x in zip(beta0, v))
        if any(beta[i] <= beta[i + 1] for i in range(len(beta) - 1)):
            break
        rec = probe(beta, T_, R_)
        if rec is None or not rec['surv'] or rec['prof'] is None:
            break
        filas.append((j, beta[0] - beta[-1], rec['prof'], rec['vac_cancelan'], defecto(beta)))
    return filas


def aritmetica(filas):
    """(True, paso) si prof crece en progresion aritmetica de paso > 0."""
    p = [f[2] for f in filas]
    if len(p) < 4:
        return False, 0
    d = p[1] - p[0]
    return (d > 0 and all(p[i] - p[i - 1] == d for i in range(1, len(p)))), d


# ===================================================================== C0 ========================
print("=" * 112)
print("C0  ACEPTACION -- fatal")
print("=" * 112)
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
    print("  C0a  probe() == scan()  t=%d r=%d M=%d : %s" % (t, r, M, "ok" if ok else "*** FALLA ***"))

f = itera(SEMILLA, (4, 2, 2, 2, 2, 2, 2, 0))
ok = len(f) == JMAX + 1 and [x[2] for x in f] == [6 + 2 * j for j in range(JMAX + 1)]
bad += not ok
print("  C0b  la familia YA publicada (v = (4,2,2,2,2,2,2,0)) se re-mide: %d miembros, prof %d..%d   %s"
      % (len(f), f[0][2] if f else -1, f[-1][2] if f else -1, "ok" if ok else "*** FALLA ***"))
print("")
if bad:
    print("  C0 FALLA -- el resto NO vale.")
    print("DONE (veredicto suspendido)")
    raise SystemExit(1)
print("  C0 PASA")

# ===================================================================== N1 ========================
print("")
print("=" * 112)
print("N1  LA PREDICCION, escrita en la cabecera ANTES de correr:  v* = (4,4,4,4,0,0,0,0)")
print("=" * 112)
print("")
print("     v* es concentrico (4+0 en los cuatro pares) y == 0 (mod 4): NO mueve las clases de residuos.")
print("     Si la ley es del DEFECTO, tiene que generar familia.")
print("")
for nombre, sem in [("semilla A", SEMILLA), ("semilla B", SEMILLA_B)]:
    filas = itera(sem, PREDICHO)
    ar, d = aritmetica(filas)
    print("  %s = %s" % (nombre, str(sem)))
    if not filas:
        print("     v* NO genera familia: se rompe en el primer paso.  LA PREDICCION FALLA.")
        continue
    print("        j    W   prof  cancelan  defecto")
    for (j, W, p, c, dd) in filas:
        print("     %4d %4d %6d %9d   %s" % (j, W, p, c, str(dd)))
    defs = {x[4] for x in filas}
    print("     -> %d miembros; prof %s; defecto %s"
          % (len(filas), "en progresion aritmetica de paso %d  LA PREDICCION ACIERTA" % d if ar
             else "NO en progresion aritmetica  LA PREDICCION FALLA A MEDIAS",
             "CONSTANTE %s" % str(defs.pop()) if len(defs) == 1 else "VARIA: %s" % sorted(defs)))
    print("")
sys.stdout.flush()

# ===================================================================== N2 ========================
print("=" * 112)
print("N2  EL CONO ENTERO contra EL SENUELO -- el control que decide si la causa es la concentricidad")
print("=" * 112)
print("")
todos = [v for v in itertools.product(range(VMAX + 1), repeat=N_)
         if all(v[i] >= v[i + 1] for i in range(N_ - 1)) and sum(v) > 0]
cono = [v for v in todos if concentrico(v)]
senuelo = [v for v in todos if not concentrico(v)]
print("     vectores no-crecientes con entradas 0..%d : %d, de los cuales %d concentricos y %d no"
      % (VMAX, len(todos), len(cono), len(senuelo)))
print("")
res = {}
for etiqueta, lista in [("CONO (v concentrico)", cono), ("SENUELO (v NO concentrico)", senuelo)]:
    viven = []
    defvar = 0
    for v in lista:
        filas = itera(SEMILLA, v, jmax=8)
        ar, d = aritmetica(filas)
        if ar:
            viven.append((v, len(filas), d, filas[-1][3]))
            if len({x[4] for x in filas}) > 1:
                defvar += 1
    res[etiqueta] = (len(lista), len(viven), viven, defvar)
    print("  %-28s : %d de %d generan familia (prof en progresion aritmetica, >=4 miembros)"
          % (etiqueta, len(viven), len(lista)))
    for (v, n, d, c) in viven[:14]:
        print("       %-26s  %2d miembros, prof +%d por paso, cancelan %d" % (str(v), n, d, c))
    if len(viven) > 14:
        print("       ... y %d mas" % (len(viven) - 14))
    if viven:
        print("       de esas, %d tienen el defecto VARIANDO (la ley predice 0 en el cono)" % defvar)
    print("")
sys.stdout.flush()

# ===================================================================== N3 ========================
print("=" * 112)
print("N3  VEREDICTO")
print("=" * 112)
print("")
nc, vc, _, dvc = res["CONO (v concentrico)"]
ns, vs, _, dvs = res["SENUELO (v NO concentrico)"]
print("     cono     : %d de %d concentricos generan familia  (%.0f%%)" % (vc, nc, 100.0 * vc / nc))
print("     senuelo  : %d de %d NO concentricos generan familia  (%.0f%%)"
      % (vs, ns, 100.0 * vs / ns if ns else 0))
print("     defecto que varia dentro del cono: %d  (la ley predice 0)" % dvc)
print("")
if dvc == 0 and vc and (not ns or vc / float(nc) > 3 * vs / float(ns)):
    print("     LA LEY AGUANTA: el defecto es invariante en TODO el cono, y ser concentrico separa")
    print("     de verdad -- el senuelo no genera familias al mismo ritmo.  Luego 'para todo K existe")
    print("     un beta de profundidad > K' no es una patologia: es que el lugar de ceros tiene un")
    print("     CONO de direcciones que lo aproximan sin tocarlo, y la profundidad es el orden de")
    print("     contacto a lo largo de cada rayo.")
elif dvc:
    print("     LA LEY ES FALSA COMO ESTA ESCRITA: hay %d vectores concentricos con defecto variable." % dvc)
else:
    print("     EL SENUELO GENERA FAMILIAS AL MISMO RITMO: la concentricidad de v no es la causa, y")
    print("     la ley sobra.  Se dice, que para eso esta la columna.")
print("")
print("     Lo que esto NO es: una prueba.  El defecto invariante SI esta probado (dos lineas de")
print("     algebra, arriba), pero 'defecto constante => la familia sigue siendo superviviente' es")
print("     MEDIDO, no probado: el criterio de superviviente no es solo el defecto.")
print("")
json.dump({k: [v[0], v[1], [[list(x[0]), x[1], x[2], x[3]] for x in v[2]], v[3]]
           for k, v in res.items()}, open("defect_cone_CONE.json", "w"), indent=1)
print("DONE")
