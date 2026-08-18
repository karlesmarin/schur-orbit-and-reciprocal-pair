# -*- coding: utf-8 -*-
# CADA NUMERO DEL PAPER, CONTRA SU VOLCADO.
#
# El inventario de _claimaudit.py lista afirmaciones para leerlas a mano.  Esto hace lo contrario y
# es lo que de verdad ataja el fallo: coge los numeros que los gates de hoy dejaron en sus JSON y
# comprueba que el paper los cita tal cual.  Un numero que cambio al ampliar una medida y no se
# actualizo en el texto sale aqui, y no sale de ninguna otra manera.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _dumpcheck.py

import io
import json
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

G = os.path.join("..", "gates")
tex = io.open("orbit_pair_ii.tex", encoding="utf-8").read()
plano = tex.replace("\\,", "").replace(",", "")


def hay(n):
    """.aparece el entero n en el paper, como numero suelto?"""
    return re.search(r"(?<!\d)%d(?!\d)" % n, plano) is not None


def carga(nombre):
    p = os.path.join(G, nombre)
    if not os.path.exists(p):
        return None
    return json.load(open(p))


CHEQUEOS = []


def chequea(etiqueta, valores):
    for v in valores:
        CHEQUEOS.append((etiqueta, int(v), hay(int(v))))


d = carga("divided_differences_DUMP.json")
if d:
    chequea("divided_differences C0", [d["C0"][0], d["C0"][1]])
    chequea("divided_differences senuelos",
            [d["senuelos"]["par"][0], d["senuelos"]["sin_signo"][0],
             d["senuelos"]["solo_primero"][0]])
    chequea("divided_differences hist",
            [d["hist_terminos"].get("1", 0), d["hist_terminos"].get("2", 0),
             d["hist_terminos"].get("4", 0)])

d = carga("_probe_cancelacion_DUMP.json")
if d:
    chequea("cancelacion caja", [d["P2"][0], d["P2"][1]])

d = carga("_probe_involucion_DUMP.json")
if d:
    chequea("involucion I3/I4", [d["I3"][0], d["I3"][1], d["I4"][0], d["I4"][1]])

d = carga("_probe_weyl_doble_DUMP.json")
if d:
    chequea("weyl doble trivial", [d["trivial"], d["tot"]])

d = carga("_probe_afin_DUMP.json")
if d:
    a6 = d.get("A6_det_vs_delta", {})
    chequea("afin A6", [v for v in a6.values()])
    chequea("afin A5", [d["A5"][0]])

d = carga("_probe_toggle_DUMP.json")
if d:
    chequea("toggle", [d["T1"], d["n_pares"], d["T7_completo"][0]])
    chequea("toggle formas", list(d["forma"].values()))
    chequea("toggle signos", list(d["signos_en_el_toggle"].values()))

d = carga("_probe_multiplete_DUMP.json")
if d:
    chequea("multiplete", [d["M0"][0], d["M0"][1], d["N1_senuelo"][0], d["M2"][0], d["M2"][1]])

d = carga("unimodularidad_barrido_DUMP.json")
if d:
    chequea("barrido", [d["D0"][0], d["D0"][1], d["TU"][0], d["TU"][1], d["P0"][0]])
    chequea("barrido senuelos",
            [d["senuelos"]["par"][1], d["senuelos"]["sin_signo"][0], d["senuelos"]["sin_signo"][1]])

d = carga("L3_caja_grande_DUMP.json")
if d:
    for bloque in d:
        chequea("L3 t=%d" % bloque["t"], [bloque["L3"], bloque["L2"]])

print("NUMEROS DE LOS VOLCADOS DE HOY, BUSCADOS EN EL PAPER")
print("=" * 92)
faltan = [c for c in CHEQUEOS if not c[2]]
for (et, v, ok) in CHEQUEOS:
    if not ok:
        print("   FALTA  %-30s %d" % (et, v))
print("")
print("  comprobados: %d   presentes: %d   ausentes: %d"
      % (len(CHEQUEOS), len(CHEQUEOS) - len(faltan), len(faltan)))
print("")
print("  NOTA: 'ausente' no es siempre un fallo -- hay cifras de los volcados que el paper no cita")
print("  a proposito.  Es una lista para MIRAR, no un veredicto.")
