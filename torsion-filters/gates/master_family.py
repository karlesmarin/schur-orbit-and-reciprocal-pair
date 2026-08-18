# -*- coding: utf-8 -*-
# LA FAMILIA MAESTRA: ¿es CANCEL eventualmente CONSTANTE o solo eventualmente PERIODICA?
# 15 de agosto de 2026.
#
# LA CORRECCION QUE ORIGINA ESTE GUION, y no es nuestra.  Habiamos medido la saturacion con pasos
# s=26 y s=28 en t=4 y concluido "satura en 4 en los dos regimenes".  La consulta externa observa que
# esas NO son dos familias distintas: son dos MUESTREOS de la misma semirrecta beta(n) = beta0 + n*v.
#
#     28j = 0 (mod 4)        -> solo ve el estado residual rho = 0
#     26j = 0,2,0,2 (mod 4)  -> ve rho = 0 y rho = 2
#
# Asi que lo que teniamos demostrado era c_0 = c_2 = 4, y NO habiamos mirado c_1 ni c_3.  Tiene razon.
#
# LO QUE SE MIDE.  La familia maestra con paso 1, que visita TODOS los estados residuales:
#
#     beta(n): el maximo de S sube n y el minimo baja n,  n = 0,1,2,3,...
#
# y para cada n se registra si sigue siendo superviviente y cuanto vale CANCEL.  Dos desenlaces:
#
#   EVENTUALMENTE PERIODICA con valores distintos por clase  ->  es el marco EQP estandar y la
#       constancia que creiamos ver era un artefacto de nuestro muestreo;
#   EVENTUALMENTE CONSTANTE en TODAS las clases              ->  "rigidez residual", que NO sale del
#       marco EQP y seria el contenido nuevo.
#
# CONTROL.  Se imprime la tabla entera por n, no un resumen, y se agrupa por n mod t al final.  Si
# alguna clase no tiene ni un solo superviviente, se dice: una clase vacia no es una clase que
# coincida.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run: python master_family.py

import os
from collections import defaultdict

from second_stratum import setup

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "k_vs_m.py")
_src = open(SRC, encoding="utf-8").read()
_head = _src.split("# ===================================================================== C0")[0]
_ns = {"__name__": "k_vs_m_preamble"}
exec(compile(_head, SRC, "exec"), _ns)
probe = _ns["probe"]


def extremos_S(beta, t):
    st = setup(beta, t)
    if st is None:
        return None
    cl, E, Cd = st
    if not E:
        return None
    S = sorted({v for k in E for v in Cd[k]})
    return S[-1], S[0]


CASOS = [(4, 2, (18, 17, 11, 8, 7, 6, 1, 0), 60),
         (6, 3, (16, 15, 14, 13, 11, 6, 5, 4, 3, 2, 1, 0), 48)]

print("=" * 100)
print("FAMILIA MAESTRA paso 1: ¿CANCEL eventualmente CONSTANTE, o solo periodica por clase?")
print("=" * 100)
for (t, r, seed, NMAX) in CASOS:
    hi, lo = extremos_S(seed, t)
    print("\n--- t=%d r=%d  semilla %s   (hi=%d, lo=%d)" % (t, r, seed, hi, lo))
    porclase = defaultdict(list)
    nosurv = defaultdict(int)
    for n in range(NMAX + 1):
        b = tuple(sorted([(x + n if x == hi else (x - n if x == lo else x)) for x in seed],
                         reverse=True))
        # OJO -- aqui estuvo el fallo del primer intento y se deja escrito para no repetirlo: puse
        # una guarda `min(b) < 0` y descarto TODOS los n >= 1, porque el minimo de S baja a negativo.
        # La familia original de saturation.py tampoco traslada y funciona: un beta con entradas
        # negativas sigue siendo un conjunto de enteros distintos y la maquinaria lo admite.  Con la
        # guarda puesta sobrevivia solo n=0 y el guion IMPRIMIA IGUAL un veredicto de rigidez sobre
        # una poblacion de uno.
        if len(set(b)) != len(b):
            continue
        rec = probe(b, t, r)
        rho = n % t
        if rec is None or not rec['surv'] or rec['prof'] is None:
            nosurv[rho] += 1
            continue
        porclase[rho].append((n, rec['vac_cancelan'], rec['prof']))
    print("    rho | n con superviviente | CANCEL observados (n: valor)            | cola estable")
    valores_cola = {}
    for rho in range(t):
        datos = porclase[rho]
        if not datos:
            print("    %3d | %19d | (ninguno; %d descartados)                  | --"
                  % (rho, 0, nosurv[rho]))
            continue
        muestra = '  '.join('%d:%d' % (n, c) for n, c, _ in datos[:8])
        cola = sorted({c for _, c, _ in datos[len(datos) // 2:]})
        valores_cola[rho] = cola
        print("    %3d | %19d | %-38s | %s"
              % (rho, len(datos), muestra[:38], cola))
    POBLACION_MINIMA = 5
    vacias = [rho for rho in range(t) if not porclase[rho]]
    flacas = [rho for rho in valores_cola if len(porclase[rho]) < POBLACION_MINIMA]
    if vacias:
        print("    *** las clases %s NO TIENEN NI UN SUPERVIVIENTE en el rango: c_rho no esta"
              " definido ahi." % vacias)
        print("        No se dice 'constante en todas las clases' con clases vacias: una clase"
              " vacia no coincide, falta.")
    if flacas:
        print("    *** NO SE DECLARA VEREDICTO: las clases %s tienen menos de %d supervivientes."
              % (flacas, POBLACION_MINIMA))
        print("        Una clase con dos puntos no es una clase estabilizada.")
    elif valores_cola:
        estables = all(len(v) == 1 for v in valores_cola.values())
        vals = {rho: v[0] for rho, v in valores_cola.items() if len(v) == 1}
        print("    -> clases con cola de un solo valor: %d de %d"
              % (sum(1 for v in valores_cola.values() if len(v) == 1), len(valores_cola)))
        if estables:
            if len(set(vals.values())) == 1:
                print("       VEREDICTO: CANCEL eventualmente CONSTANTE, valor %d, en las %d"
                      " clases QUE TIENEN supervivientes (%s)"
                      % (list(vals.values())[0], len(vals), sorted(vals)))
            else:
                print("       VEREDICTO: eventualmente PERIODICA, vector c = %s"
                      "  ->  es el marco EQP, no hay rigidez" % vals)
        else:
            print("       VEREDICTO: alguna clase no se ha estabilizado en el rango barrido")

print()
print("=" * 100)
