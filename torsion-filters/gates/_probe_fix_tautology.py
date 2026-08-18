# -*- coding: utf-8 -*-
# .QUE ARREGLO SIRVE PARA SOCRATES?  UNA SOLA PRUEBA, TODOS LOS CANDIDATOS.
#
# `socrates.tools.tautology.is_tautology` falla en 9 de 16 identidades ciclotomicas del Paper II:
# su `sympy.simplify` no cierra raices de orden 3, 5, 6, 7, 9, 10 (si cierra 4 y 8, porque SymPy
# autoevalua esas a I y a sqrt(2)/2(1+I)).
#
# Antes de tocar la herramienta hay que saber QUE arreglo sirve, y no averiguarlo a base de sondas
# sueltas.  Esto compara de una vez cuatro candidatos sobre TODOS los casos, identidades y señuelos:
#
#   A  simplify            -- lo que hace hoy, como linea base
#   B  equals(0)           -- el test de igualdad de SymPy, probabilistico
#   C  coeficientes con equals(0)  -- expandir en la variable libre y decidir coeficiente a coeficiente
#   D  minimal_polynomial  -- exacto: un numero algebraico es 0 sii su polinomio minimo es x
#
# El criterio no es "cual dice TRIVIAL mas veces" sino cual acierta en las dos direcciones: TRIVIAL
# en las identidades y NO TRIVIAL en los señuelos.  Un arreglo que ademas apruebe los señuelos es
# peor que el fallo que arregla.
#
# Se mide tambien el TIEMPO, porque un decisor exacto que tarda minutos no sirve dentro de un bucle.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python _probe_fix_tautology.py

import io
import json
import sys
import time

import sympy
from sympy import I, exp, pi, symbols

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

z, x = symbols("z x")


def zeta(t):
    return exp(2 * pi * I / t)


def orbita(t):
    return sympy.prod([1 - zeta(t) ** k * z for k in range(t)]) - (1 - z ** t)


def cross(t):
    mp = (t - 1) // 2
    P = (1 - z) * sympy.prod([(1 - zeta(t) ** i * z) * (1 - zeta(t) ** (-i) * z)
                              for i in range(1, mp + 1)])
    return P - (1 - z ** t)


def even(t):
    m = (t - 2) // 2
    P = sympy.prod([(1 - zeta(t) ** i * z) * (1 - zeta(t) ** (-i) * z)
                    for i in range(1, m + 1)])
    return sympy.expand(P * (1 - z) * (1 + z)) - (1 - z ** t)


CASOS = []
for t in (3, 4, 5, 7, 8, 9):
    CASOS.append(("orbita t=%d" % t, orbita(t), True))
for t in (3, 5, 7, 9):
    CASOS.append(("cross t=%d" % t, cross(t), True))
for t in (4, 6, 8, 10):
    CASOS.append(("even t=%d" % t, even(t), True))
# señuelos: variantes deliberadamente falsas
CASOS.append(("SEÑUELO orbita t=5 con 1-z^6",
              sympy.prod([1 - zeta(5) ** k * z for k in range(5)]) - (1 - z ** 6), False))
CASOS.append(("SEÑUELO cross t=7 sin (1-z)",
              cross(7) / (1 - z), False))
CASOS.append(("SEÑUELO even t=6 sin (1+z)",
              sympy.expand(sympy.prod([(1 - zeta(6) ** i * z) * (1 - zeta(6) ** (-i) * z)
                                       for i in range(1, 3)]) * (1 - z)) - (1 - z ** 6), False))
CASOS.append(("SEÑUELO z - z**2", z - z ** 2, False))


def A_simplify(e):
    return sympy.simplify(e) == 0


def B_equals(e):
    try:
        return bool(e.equals(0))
    except Exception:
        return False


def C_coefs_equals(e):
    try:
        p = sympy.Poly(sympy.expand(e), z)
    except Exception:
        return B_equals(e)
    for c in p.all_coeffs():
        if c == 0:
            continue
        try:
            if not bool(c.equals(0)):
                return False
        except Exception:
            return False
    return True


def D_minpoly(e):
    try:
        p = sympy.Poly(sympy.expand(e), z)
    except Exception:
        return False
    for c in p.all_coeffs():
        if c == 0:
            continue
        if getattr(c, "is_number", False) is not True:
            return False
        try:
            mp = sympy.minimal_polynomial(c, x)
        except Exception:
            return False
        if mp != x:
            return False
    return True


CAND = [("A simplify (hoy)", A_simplify), ("B equals(0)", B_equals),
        ("C coefs+equals", C_coefs_equals), ("D minimal_polynomial", D_minpoly)]

print("=" * 104)
print("QUE ARREGLO SIRVE: cuatro candidatos contra %d casos (identidades y señuelos)" % len(CASOS))
print("=" * 104)
print("")
print("  %-30s %-10s %-10s %-10s %-10s" % ("caso", *[c[0][:10] for c in CAND]))
print("  " + "-" * 84)

marcador = {n: [0, 0.0] for (n, _) in CAND}
for (etq, e, espera) in CASOS:
    fila = []
    for (n, f) in CAND:
        t0 = time.time()
        try:
            v = bool(f(e))
        except Exception:
            v = None
        dt = time.time() - t0
        marcador[n][1] += dt
        acierta = (v == espera)
        marcador[n][0] += 1 if acierta else 0
        fila.append(("ok " if acierta else "!! ") + ("triv" if v else "no"))
    print("  %-30s %-10s %-10s %-10s %-10s" % (etq[:30], *fila))

print("")
print("  " + "-" * 84)
for (n, _) in CAND:
    a, dt = marcador[n]
    print("  %-24s aciertos %2d de %2d      tiempo total %6.1f s" % (n, a, len(CASOS), dt))
print("")
print("  CRITERIO: gana el que acierte en LAS DOS direcciones.  Un arreglo que tambien apruebe los")
print("  señuelos es peor que el fallo que arregla.")

json.dump({n: {"aciertos": marcador[n][0], "de": len(CASOS), "segundos": round(marcador[n][1], 2)}
           for (n, _) in CAND},
          io.open("_probe_fix_tautology_DUMP.json", "w", encoding="utf-8"), indent=1)
print("")
print("DONE")
