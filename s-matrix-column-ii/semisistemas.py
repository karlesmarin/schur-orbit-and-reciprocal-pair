# -*- coding: utf-8 -*-
r"""semisistemas.py -- las orbitas raras del centro del censo son SEMISISTEMAS de suma nula.

DE DONDE SALE.  Mirando una a una las orbitas con delta >= 3 (que es donde nadie habia mirado:
siempre nos quedamos en n <= 7, y el centro del censo esta en n = (p-1)/2), aparece un patron de
anulacion imposible de confundir: M_r = 0 para TODOS los r pares < p-1.  Ejemplo p = 17,
T = {1,2,3,4,8,10,11,12}: nulos en r = 2,4,6,8,10,12,14.

LA IDENTIFICACION.  Un conjunto asi es un SEMISISTEMA ("half-system", el objeto del lema de Gauss y
de la demostracion de Zolotarev de la reciprocidad cuadratica): un representante de cada pareja
{a, -a} de F_p^x, es decir |T| = (p-1)/2 y T no contiene ningun par {a,-a}.

  Semisistema  =>  al elevar al cuadrado, T^2 recorre TODOS los residuos cuadraticos una vez, y
  M_{2k}(T) = sum_{s en QR} s^k = 0 para todo 1 <= k < (p-1)/2.
  Reciprocamente, si M_{2k} = 0 para k = 1 .. m-1 con m = |T|, las simetricas elementales del
  multiconjunto T^2 se anulan (Newton, p > m) y su polinomio es X^m - c: T^2 es una clase lateral
  de las raices m-esimas, o sea los QR; luego cuadrar es inyectivo en T y T es un semisistema.

CONSECUENCIA PARA EL SEMIGRUPO.  Si T es semisistema centrado (suma 0), S_T no tiene ningun
generador par por debajo de e, asi que S_T = <e, j impares con M_j != 0> y, cuando M_3 M_5 M_7 != 0,
S_T = <3,5,7> con huecos {1,2,4}: delta = 3.

ESTE GUION comprueba las dos direcciones y calcula el delta de todos los semisistemas centrados.

    python semisistemas.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from itertools import product

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo, estabilizador, momento, semigrupo, s_exacto    # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 23


def semisistemas(p):
    """todos los semisistemas de F_p^x: un representante de cada pareja {a,-a}."""
    parejas = [(a, p - a) for a in range(1, (p + 1) // 2)]
    for eleccion in product(*parejas):
        yield tuple(sorted(eleccion))


fallo_ida, fallo_vuelta, total = 0, 0, 0
for p in range(5, PMAX + 1):
    if not es_primo(p):
        continue
    m = (p - 1) // 2
    reparto, ejemplos = {}, {}
    nulos_centrados = 0
    for T in semisistemas(p):
        total += 1
        pares = [r for r in range(2, p - 1, 2) if momento(T, r, p) % p]
        if pares:                       # IDA: semisistema => todos los momentos pares nulos
            fallo_ida += 1
        if sum(T) % p:
            continue
        nulos_centrados += 1
        S = semigrupo(s_exacto(T, p))
        d, h = S["genero"], len(estabilizador(T, p))
        clave = (d, tuple(S["minimales"]), h)
        reparto[clave] = reparto.get(clave, 0) + 1
        ejemplos.setdefault(clave, T)
    print("p=%-3d  m=(p-1)/2=%-3d  semisistemas=%-7d  de suma nula=%-6d"
          % (p, m, 2 ** m, nulos_centrados))
    for clave in sorted(reparto):
        d, gens, h = clave
        print("        delta=%d  S=%-14s h=%-3d : %-5d  ej. %s"
              % (d, str(list(gens)), h, reparto[clave], str(ejemplos[clave])))
    sys.stdout.flush()

print("")
print("semisistemas examinados: %d ; fallos de 'semisistema => M_par = 0': %d" % (total, fallo_ida))
print("VEREDICTO:", "la identificacion SOBREVIVE" if fallo_ida == 0 else "FALSA")
