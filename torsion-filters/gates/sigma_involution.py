# -*- coding: utf-8 -*-
"""LA INVOLUCION sigma -- el mecanismo de la direccion SUFICIENTE del criterio, para todo t y r.

QUE SE MIDE.  Sobre las beta con S simetrico bajo sigma(x) = C - x, separadas por si hay o no un
incremento igual a C, se comprueban las cuatro piezas del argumento:

    (a) sigma manda transversal a transversal        -- refleja SOLO las clases de exceso; las
                                                        singleton estan forzadas y no se tocan
    (b) sigma conserva el grado de cada transversal
    (c) w(sigma(g)) = -w(g)
    (d) sigma no tiene transversales fijas

y si Phi == 0.  Resultado: con incremento = C las cuatro se cumplen al 100% y Phi == 0 al 100%;
SIN incremento = C, sigma sigue actuando pero SIEMPRE tiene punto fijo, y Phi nunca se anula.

O sea las dos clausulas del criterio hacen dos trabajos distintos, y este guion los separa.

DOS ERRORES MIOS QUE ESTE GUION CAZO, y por eso existen sus dos controles:
  * la primera version reflejaba tambien las clases SINGLETON, que no estan en S y no tienen por
    que ser simetricas: sigma salia mal definida en 44 de 127 casos.  No era una refutacion, era
    un fallo de la definicion.
  * la primera version pasaba solo 3 puntos de evaluacion con r = 4, o sea evaluaba un alfabeto de
    t+6 letras contra una matriz de t+8: dos anomalias falsas en t=2, r=4.

Authors: Carles Marin, Claude (AI assistant).
"""
import sys
sys.path.insert(0, r"E:\proyectos\Curiosity\research\orbit-pair\gates")
from fractions import Fraction as F
from second_stratum import setup, all_transversals
from criterion_control import betas, lam_of, value

PTS = [[F(3, 2), F(5, 3), F(7, 4), F(11, 6)], [F(5, 2), F(7, 3), F(9, 4), F(13, 5)],
       [F(4, 3), F(9, 5), F(11, 6), F(15, 7)], [F(7, 5), F(11, 4), F(13, 7), F(17, 9)]]


def is_zero(b, t, r):
    lam = lam_of(list(b))
    for p in PTS:
        if value(lam, t, r, p[:r]) != 0:
            return False
    return True


def anat(b, t, r):
    st = setup(b, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted(x for k in E for x in Cd[k])
    C = S[0] + S[-1]
    inc = []
    for k in E:
        c = sorted(Cd[k], reverse=True)
        inc += [c[i] + c[i + 1] for i in range(len(c) - 1)]
    return cl, E, Cd, S, C, inc


CFG = [(4, 2, 15), (6, 2, 15), (8, 2, 15), (4, 3, 13), (6, 3, 13), (2, 2, 14), (2, 3, 13), (2, 4, 13)]
print("=" * 110)
print("sigma CORREGIDA: refleja SOLO las clases de exceso; las singleton estan forzadas y no se tocan.")
print("Poblacion: beta con S C-simetrico.  Se separa por si hay incremento = C, que es la otra mitad")
print("del criterio, para ver QUE hace exactamente esa condicion.")
print("=" * 110)
print("           |------------- CON incremento = C -------------|------ SIN incremento = C ------|")
print("  t  r      n   sigma ok  grado ok   w=-w   sin fijos  Phi=0 |    n   sigma ok  sin fijos  Phi=0")
for (t, r, W) in CFG:
    fila = []
    for tiene_incr in (True, False):
        n = ok1 = ok2 = ok3 = nofix = nz = 0
        for b in betas(t, r, W):
            d = anat(b, t, r)
            if d is None:
                continue
            cl, E, Cd, S, C, inc = d
            if sorted(C - x for x in S) != S:
                continue
            if (C in inc) != tiene_incr:
                continue
            n += 1
            tr = all_transversals(b, cl, r, t)
            idx = {tuple(sorted(x[0].items())): x for x in tr}
            o1 = o2 = o3 = True
            fijo = False
            for (sel, T, w, deg) in tr:
                ssel = {}
                for k, v in sel.items():
                    if k in E:
                        ssel[(C - k) % t] = C - v
                    else:
                        ssel[k] = v            # singleton: forzada, no se refleja
                key = tuple(sorted(ssel.items()))
                if key not in idx:
                    o1 = False
                    break
                y = idx[key]
                if y[3] != deg:
                    o2 = False
                if key == tuple(sorted(sel.items())):
                    fijo = True
                elif y[2] != -w:
                    o3 = False
            ok1 += o1; ok2 += o2; ok3 += o3; nofix += (not fijo)
            nz += is_zero(b, t, r)
        fila.append((n, ok1, ok2, ok3, nofix, nz))
    (n1, a1, b1, c1, d1, z1), (n0, a0, b0, c0, d0, z0) = fila
    print("  %2d %2d %6d %8d %9d %6d %10d %6d | %4d %9d %10d %6d"
          % (t, r, n1, a1, b1, c1, d1, z1, n0, a0, d0, z0))
