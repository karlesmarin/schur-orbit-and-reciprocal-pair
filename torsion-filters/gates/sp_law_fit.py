# -*- coding: utf-8 -*-
# ¿DE QUE DEPENDE mu_max?  El ajuste, con su control de vacuidad.  15 de agosto de 2026.
#
# ENTRADA: sp_law_DUMP.json, las 435 formas criticas con g_com asimetrico y mu_max unico, calculadas
# en Sage (sp_law.sage) y validadas contra el bialternante (sp_expansion.sage, C0).
#
# LO QUE SE PRUEBA
#   N1  formulas cerradas ingenuas para mu_max.  Ninguna pasa, y se dice el porcentaje de cada una en
#       vez de callarlas: un 51 % es informacion, no un fracaso mudo.
#   N2  la pregunta principiada, no adivinada: Littlewood dice que el VALOR de s_lambda(mu_t) lo fijan
#       el t-core y el t-cociente.  ¿Es mu_max funcion del cociente?
#
# EL CONTROL QUE ESTE GUION EXISTE PARA NO SALTARSE
#   C1  VACUIDAD.  "ninguna clase tiene dos mu distintos" es una afirmacion VACIA si cada clase tiene
#       un solo elemento.  Hay que imprimir el tamaño de las clases ANTES de leer el veredicto.  La
#       primera version de este analisis dio "SI es funcion del cociente" con 435 clases para 435
#       formas -- o sea, sin un solo par que comparar.  Es el mismo error que ya nos costo una vuelta
#       con el empate de CANCEL: un empate no es un acuerdo, y una clase de un elemento no es un test.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python sp_law_fit.py

import json
from collections import defaultdict

D = json.load(open('sp_law_DUMP.json'))
t, N = 4, 8


def quotient(beta, t):
    """(t-cociente, vector de ocupacion).  La ocupacion determina el t-core."""
    run = defaultdict(list)
    for v in beta:
        run[v % t].append((v - v % t) // t)
    qs = []
    for k in range(t):
        b = sorted(run.get(k, []), reverse=True)
        lam = tuple(b[i] - (len(b) - 1 - i) for i in range(len(b)))
        qs.append(tuple(x for x in lam if x > 0))
    return tuple(qs), tuple(len(run.get(k, [])) for k in range(t))


print("=" * 100)
print("N1  FORMULAS CERRADAS PARA mu_max -- ninguna pasa, y se dice cuanto acierta cada una")
print("=" * 100)
print("")
inv = [("lambda_1", lambda d: d['lam'][0]),
       ("(maxS - minS)//2", lambda d: (d['S'][-1] - d['S'][0]) // 2),
       ("C//2", lambda d: d['C'] // 2),
       ("|S|", lambda d: len(d['S'])),
       ("C//2 - 1", lambda d: d['C'] // 2 - 1),
       ("C - max(g_com)", lambda d: d['C'] - max(d['gcom'])),
       ("min(g_com)", lambda d: min(d['gcom']))]
print("    mu_0 == ...                 aciertos")
for nom, f in inv:
    ok = sum(1 for d in D if f(d) == d['mu'][0])
    print("    %-22s %4d de %d  (%3.0f %%)" % (nom, ok, len(D), 100.0 * ok / len(D)))
print("")
par = sum(1 for d in D if (sum(d['mu']) - sum(d['lam'])) % 2 == 0)
print("    |mu| = |lambda| mod 2 : %d de %d" % (par, len(D)))
print("    -- y esto NO es un hallazgo: z <-> 1/z mueve el grado de dos en dos, luego la paridad")
print("       del grado total esta fijada de antemano.  Se deja escrito para no citarlo como ley.")
print("")
print("    A = +1 en %d,  A = -1 en %d" % (sum(1 for d in D if d['A'] == 1),
                                           sum(1 for d in D if d['A'] == -1)))

print("")
print("=" * 100)
print("N2 + C1  ¿ES mu_max FUNCION DEL t-COCIENTE?  -- con el control de vacuidad DELANTE")
print("=" * 100)
print("")
for etiqueta, clave in [("cociente + ocupacion", lambda d: quotient(d['beta'], t)),
                        ("solo el cociente", lambda d: quotient(d['beta'], t)[0]),
                        ("solo la ocupacion", lambda d: quotient(d['beta'], t)[1])]:
    mapa = defaultdict(set)
    tam = defaultdict(int)
    for d in D:
        k = clave(d)
        mapa[k].add(tuple(d['mu']))
        tam[k] += 1
    comparables = sum(1 for k in tam if tam[k] > 1)
    choques = sum(1 for k, v in mapa.items() if len(v) > 1)
    print("    agrupando por %-22s : %3d clases, %3d con MAS DE UNA forma" % (etiqueta, len(mapa), comparables))
    if not comparables:
        print("       *** VACIO: cada clase tiene una sola forma.  NO HAY NADA QUE COMPARAR y no se")
        print("           puede decir 'si es funcion'.  Esta fila no dice nada. ***")
    else:
        print("       clases con mas de un mu_max : %d  ->  %s"
              % (choques, "es funcion" if not choques else "NO es funcion"))
    print("")

print("=" * 100)
print("N3  DONDE VIVE LO QUE FALTA -- lambda es 4-core (cociente vacio), y mu_max sigue moviendose")
print("=" * 100)
print("")
print("    Con el cociente VACIO, lo unico que queda es que clases residuales son las de exceso.")
print("")
print("    ocupacion    | clases de exceso | mu_max")
vis = set()
filas = []
for d in D:
    q, occ = quotient(d['beta'], t)
    if any(q):
        continue
    E = tuple(k for k in range(t) if occ[k] >= 2)
    if (occ, E) in vis:
        continue
    vis.add((occ, E))
    filas.append((occ, E, tuple(d['mu'])))
for occ, E, mu in sorted(filas)[:12]:
    print("    %-12s | %-16s | %s" % (str(occ), str(E), str(mu)))
print("")
print("    mu_max se mueve con QUE clases son las de exceso, no solo con cuantas.  O sea que lo que")
print("    falta vive en el bloque de raices de la unidad, que es justo la mitad que ni nosotros ni")
print("    la consulta hemos tocado.")
print("")
print("=" * 100)
