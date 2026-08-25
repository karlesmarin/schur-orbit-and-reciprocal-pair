# -*- coding: utf-8 -*-
"""Deja un .sage REPL-safe: ninguna linea en blanco dentro de NINGUN bloque compuesto de nivel
superior (def, for, while, if, with, class, try).  Mi primer verificador solo miraba `def`, y por eso
se colo una linea en blanco dentro de un `for` que corto el bucle y dejo correr solo la ultima
iteracion.  Uso:  python repl_safe.py <fichero> [--fix]"""
import io
import re
import sys

ABRE = re.compile(r"^(def|for|while|if|with|class|try)\b")


def analiza(lineas):
    """Devuelve las posiciones (1-based) de lineas en blanco dentro de un bloque de nivel superior."""
    malas = []
    dentro = False
    for k, l in enumerate(lineas):
        if ABRE.match(l):
            dentro = True
            continue
        if l and not l.startswith((" ", "\t")):
            dentro = False
            continue
        if dentro and l.strip() == "":
            nxt = next((x for x in lineas[k + 1:] if x.strip() != ""), "")
            if nxt.startswith((" ", "\t")):
                malas.append(k + 1)
            else:
                dentro = False
    return malas


p = sys.argv[1]
fix = "--fix" in sys.argv
lineas = io.open(p, encoding="utf-8").read().split("\n")
malas = analiza(lineas)
print("%s: %d lineas en blanco dentro de bloques de nivel superior%s"
      % (p, len(malas), ("  -> %s" % malas[:12]) if malas else ""))
tope = [k + 1 for k, l in enumerate(lineas) if re.match(r"^(else|elif|except|finally)\b", l)]
print("   else/elif/except a nivel superior: %s" % (tope if tope else "ninguno"))
if fix and malas:
    fuera = set(malas)
    nuevo = [l for k, l in enumerate(lineas) if (k + 1) not in fuera]
    io.open(p, "w", encoding="utf-8", newline="").write("\n".join(nuevo))
    print("   ARREGLADO: %d lineas quitadas; quedan %d" % (len(malas), len(analiza(nuevo))))
