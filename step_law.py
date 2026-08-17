# -*- coding: utf-8 -*-
# NOTA DEL PAQUETE.  La unica diferencia con el fichero de trabajo es la linea de import:
# las rutinas compartidas viven aqui en criterion_base.py.  Todo lo demas es identico.
"""LA LEY DE PASO -- el Lema 4.4 del paper, para clases de tamano cualquiera.

El paper enuncia el movimiento de columna solo para clases de DOS elementos, donde no hay nada
entre medias.  La seccion 8 trabaja con clases de tamano cualquiera y lo necesita en general:

    g y g' coinciden fuera de una clase y eligen alli u > v.  Entonces

        w(g') / w(g) = (-1)^{1 + B + M},

    con B = #{beta_j : v < beta_j < u}  y  M = #{x en P : v < x < u},

    donde P es el conjunto CONGELADO entero -- un valor por cada una de las t clases, o sea las
    elecciones g MAS los singletons.  Leer M solo sobre g es FALSO, y no por poco: el testigo
    minimo es beta=(4,3,2,1,0) con t=3, moviendo la clase del 1 de u=4 a v=1, donde el valor 2
    esta en medio y vive en una clase singleton.  El paper lo decia mal (2026-08-15) y este
    guion ya media bien -- el codigo usaba P y esta linea decia g.

La prueba es de la definicion: la columna congelada se desplaza B+1 posiciones, que es lo que cambia
la suma de columnas; y en la palabra de residuos la letra cruza exactamente las M letras congeladas
que hay en medio, cada cruce cambiando inv en uno.  Las B-M columnas intermedias que NO estan
congeladas no aportan letra.

SENUELO, y es el que hay que batir: (-1)^{p-q}, contando solo pasos dentro de la clase.  Ignora los
dos conteos y acierta aproximadamente la mitad de las veces -- lo cual basta para colarse en un
rango corto.  Fue mi primera version y de ahi salio una descomposicion entera equivocada.

Authors: Carles Marin, Claude (AI assistant).
"""
import sys
sys.path.insert(0, r"E:\proyectos\Curiosity\research\orbit-pair\gates")
from criterion_base import setup, all_transversals
from criterion_control import betas

CFG = [(4, 2, 13), (6, 2, 13), (3, 2, 12), (2, 2, 12), (4, 3, 12), (2, 3, 12), (5, 2, 12)]
print("=" * 96)
print("LA LEY DE PASO, corregida.  g y g' difieren solo en la clase j, valores u > v.")
print("   B = #{valores de beta estrictamente entre v y u}")
print("   M = #{valores de P  estrictamente entre v y u}   (P = los congelados, uno por clase)")
print("   CANDIDATA:  w(g')/w(g) = (-1)^{1 + B + M}")
print("   SENUELO   :  w(g')/w(g) = (-1)^{p-q}   (la que acabo de refutar; debe seguir fallando)")
print("=" * 96)
print("  t  r     casos    candidata   senuelo")
TA = TB = TOT = 0
for (t, r, W) in CFG:
    n = a = b_ = 0
    for beta in betas(t, r, W):
        st = setup(beta, t)
        if st is None:
            continue
        cl, E, Cd = st
        if not E:
            continue
        tr = all_transversals(beta, cl, r, t)
        idx = {tuple(sorted(x[0].items())): x for x in tr}
        pos = {k: {v: i + 1 for i, v in enumerate(sorted(Cd[k], reverse=True))} for k in E}
        for (sel, T, w, deg) in tr:
            P = sorted(sel.values())
            for k in E:
                for v2 in Cd[k]:
                    if v2 == sel[k]:
                        continue
                    u, v = max(sel[k], v2), min(sel[k], v2)
                    B = sum(1 for x in beta if v < x < u)
                    M = sum(1 for x in P if v < x < u)
                    s2 = dict(sel); s2[k] = v2
                    y = idx[tuple(sorted(s2.items()))]
                    n += 1
                    a += (y[2] == (-1) ** (1 + B + M) * w)
                    b_ += (y[2] == (-1) ** (pos[k][sel[k]] - pos[k][v2]) * w)
    TA += a; TB += b_; TOT += n
    print("  %2d %2d %9d %12d %9d" % (t, r, n, a, b_))
print()
print("  candidata %d/%d      senuelo %d/%d" % (TA, TOT, TB, TOT))
