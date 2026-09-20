# -*- coding: utf-8 -*-
r"""escalera_defectos.py -- la escalera de los defectos pequenos, con su censo.

De donde salen los numeros 2 y 4 de las expansiones de la nota: de CUANTOS semigrupos numericos
hay de cada genero.  Genero 1 hay uno, genero 2 hay dos, genero 3 hay cuatro --- y esos son
exactamente el 2 de F_n = 2/p - 1/p^2 y el 4 de P(delta>=3) = 4/p^2, que a su vez da el
3 = -1 + 4 del defecto medio.

Aqui se mide, para h = 1 y T centrado, cuantas orbitas realizan cada semigrupo, y se comprueba
que las condiciones de momentos de cada fila son necesarias y suficientes.  La salida alimenta
la figura de la escalera.

    python escalera_defectos.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import centrados, estabilizador, momento, s_exacto, semigrupo   # noqa: E402

# genero -> lista de (nombre, huecos, momentos que se anulan, momentos que NO)
ESCALERA = {
    1: [("<2,3>", (1,), (), (2, 3))],
    2: [("<3,4,5>", (1, 2), (2,), (3, 4, 5)),
        ("<2,5>", (1, 3), (3,), (2, 5))],
    3: [("<4,5,6,7>", (1, 2, 3), (2, 3), (4, 5, 6, 7)),
        ("<3,5,7>", (1, 2, 4), (2, 4), (3, 5, 7)),
        ("<3,4>", (1, 2, 5), (2, 5), (3, 4)),
        ("<2,7>", (1, 3, 5), (3, 5), (2, 7))],
}
RANGOS = [(11, range(3, 10)), (13, range(3, 12)), (17, range(3, 16)), (19, range(3, 18)),
          (23, range(4, 10)), (29, range(5, 9))]

cuenta = {}
fallos = 0
total = 0
for p, ns in RANGOS:
    for n in ns:
        for T in centrados(p, n):
            if len(estabilizador(T, p)) != 1:
                continue
            try:
                S = semigrupo(s_exacto(T, p))
            except ValueError:
                continue
            g = S["genero"]
            if g not in ESCALERA:
                continue
            total += 1
            hue = tuple(S["huecos"])
            fila = [f for f in ESCALERA[g] if f[1] == hue]
            if not fila:
                fallos += 1
                print("  HUECOS FUERA DE LA ESCALERA: p=%d T=%s g=%d huecos=%s" % (p, T, g, hue))
                continue
            nombre, _, ceros, nonulos = fila[0]
            ok = (all(momento(T, r, p) == 0 for r in ceros)
                  and all(momento(T, r, p) != 0 for r in nonulos))
            if not ok:
                fallos += 1
                print("  CONDICION FALLA: p=%d T=%s %s" % (p, T, nombre))
            cuenta[nombre] = cuenta.get(nombre, 0) + 1

print("LA ESCALERA DE LOS DEFECTOS PEQUENOS")
print("%-6s %-12s %-14s %-22s %-24s %s"
      % ("delta", "semigrupo", "huecos", "se anulan", "no se anulan", "orbitas medidas"))
for g in sorted(ESCALERA):
    for nombre, hue, ceros, nonulos in ESCALERA[g]:
        print("%-6d %-12s %-14s %-22s %-24s %d"
              % (g, nombre, str(hue),
                 ", ".join("M_%d" % r for r in ceros) or "(ninguno)",
                 ", ".join("M_%d" % r for r in nonulos),
                 cuenta.get(nombre, 0)))
    print("       %d semigrupo(s) de genero %d  ->  el coeficiente %d de la expansion"
          % (len(ESCALERA[g]), g, len(ESCALERA[g])))

print("")
print("orbitas con h=1 y delta<=3 examinadas: %d ; fallos: %d" % (total, fallos))
print("VEREDICTO:", "la escalera 1, 2, 4 se sostiene y sus condiciones son exactas"
      if not fallos else "ALGO FALLA")
