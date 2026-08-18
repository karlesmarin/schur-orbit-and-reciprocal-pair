# -*- coding: utf-8 -*-
# Reimprime la tabla del cuadrado conmutativo DESDE commuting_square_DUMP.json.
#
# POR QUE EXISTE ESTE FICHERO.  La corrida de commuting_square.sage termino (el volcado se escribe
# en la ultima linea del guion y tiene las cinco filas), pero el _OUT.txt capturado se quedo en la
# PRIMERA fila -- 616 bytes.  O sea: el dato esta completo, el registro impreso no.  Esto no vuelve
# a calcular nada: lee el volcado y lo formatea, para que el fichero que se cita tenga las cinco.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _reprint_square.py
import json

D = json.load(open("commuting_square_DUMP.json"))
L = []
L.append("=" * 122)
L.append("EL CUADRADO CONMUTATIVO   --   t=4, r=2, N=8")
L.append("=" * 122)
L.append("")
L.append("  (tabla reconstruida de commuting_square_DUMP.json: la corrida termino, el stdout")
L.append("   capturado se quedo en la primera fila.  Ningun numero se recalcula aqui.)")
L.append("")
L.append("  beta                         | lambda            | #nu | eps 0,+-1 | C0 izq | C0 der | (*) coef a coef | señuelo t'")
L.append("  " + "-" * 118)
for d in D:
    L.append("  %-28s | %-17s | %3d | %-9s | %-6s | %-6s | %-15s | %s" % (
        str(tuple(d["beta"])), str(tuple(d["lambda"])), d["n_eps"],
        "si" if d["eps_pm1"] else "NO",
        "ok" if d["C0_izq"] else "FALLA",
        "ok" if d["C0_der"] else "FALLA",
        "IGUAL" if d["igual"] else "DISCREPA",
        "discrepa" if d["senuelo_discrepa"] else "EMPATA"))
L.append("  " + "-" * 118)
n = len(D)
L.append("  TOTAL  %d formas | (*) IGUAL en %d | C0 izq %d | C0 der %d | eps en {0,+-1} %d | señuelo discrepa en %d" % (
    n, sum(1 for d in D if d["igual"]), sum(1 for d in D if d["C0_izq"]),
    sum(1 for d in D if d["C0_der"]), sum(1 for d in D if d["eps_pm1"]),
    sum(1 for d in D if d["senuelo_discrepa"])))
L.append("")
L.append("=" * 122)
L.append("DONE")
open("commuting_square_OUT.txt", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
