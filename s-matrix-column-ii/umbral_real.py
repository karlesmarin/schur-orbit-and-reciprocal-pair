# -*- coding: utf-8 -*-
r"""umbral_real.py -- ¿en que primo aparece de verdad la primera columna de defecto maximo?

LA PREGUNTA (Carles, 20-sep): ¿por que el caso raro se degrada factorialmente con n?

HAY DOS FACTORIALES Y NO SON EL MISMO:

  * el de la existencia: el estrato extremo pide que el Frobenius sea TRIVIAL en S_n --- que el
    trinomio X^n + aX + b se escinda del todo --- lo que tiene densidad 1/n!.  El numero ESPERADO
    de columnas de defecto maximo es del orden de p/n!, que llega a 1 cuando p ~ n!;

  * el de la cota: Hasse-Weil enfrenta ese termino principal contra 2 g_n sqrt(p), con
    g_n = 1 + [(n-2)(n-3)-4](n-2)!/4 ~ n!/4, y exige p > (2 g_n)^2 ~ (n!/2)^2, el CUADRADO.

=================================================================================================
POR QUE ESTE GUION ESTA REESCRITO (20-sep-2026, tarde).  LEER ANTES DE TOCARLO.
=================================================================================================
La version anterior tenia PMIN = 601 por defecto, porque se escribio para comprobar la prediccion
"deberia aparecer cerca de 720 = 6!".  Encontro el 607 --- el primer primo de SU VENTANA --- y la
nota lo publico como "el primer primo con un sextinomio totalmente descompuesto".  Falso: el
primero es 163.  La ventana la habia elegido la hipotesis que el guion iba a comprobar, y el
resultado no podia sino confirmarla.

Lo mismo, por otra via, con n = 5: el 233 se leyo de la tabla de extremo_bring.py, que imprime
con el filtro "p <= 60 or not ok or p > PMAX - 20" y esconde p = 61..229.  La compuerta SI habia
calculado N(67) = 1; el filtro de impresion se convirtio en una afirmacion.  Y la PARTE I, ya
publicada, dice en su propio texto que en p = 67 el polinomio de T = {9,27,50,52,63} es
X^5 + 2X + 27 con delta = 6.

DOS REMEDIOS, y el segundo es el que importa:
  1. se barre desde p = n + 1, no desde donde convenga;
  2. el guion SE NIEGA a imprimir "primera aparicion" si el barrido no empezo en n + 1.  Un
     minimo de una ventana no es un minimo.

COMO SE CUENTA, sin factorizar: para cada raiz x y cada a, el b que hace a x raiz de X^n + aX + b
esta forzado, b = -(x^n + ax).  Una tabla (a,b) -> numero de raices DISTINTAS se llena en O(p^2)
recorriendo (x, a); la casilla vale n exactamente cuando el trinomio se escinde del todo con
raices distintas.  Se exige ab != 0, que es la condicion de la Parte I.

    python umbral_real.py [N] [PMAX] [PMIN]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import math
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo                                              # noqa: E402

NS = [int(sys.argv[1])] if len(sys.argv) > 1 else [5, 6, 7]
PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 700
PMIN = int(sys.argv[3]) if len(sys.argv) > 3 else None       # solo para depurar; ver la guarda


def pares_escindidos(n, p):
    """numero de pares (a,b) con ab != 0 tales que X^n + aX + b tiene n raices distintas."""
    x = np.arange(p, dtype=np.int64)
    xn = np.array([pow(int(t), n, p) for t in x], dtype=np.int64)
    tabla = np.zeros((p, p), dtype=np.int16)
    for a in range(1, p):
        np.add.at(tabla, (a, (-(xn + a * x)) % p), 1)
    tabla[:, 0] = 0
    return int(np.count_nonzero(tabla == n))


for N in NS:
    arranque = PMIN if PMIN is not None else N + 1
    fact = math.factorial(N)
    g = 1 + ((N - 2) * (N - 3) - 4) * math.factorial(N - 2) // 4
    print("=" * 92)
    print("n = %d ; n! = %d ; genero g_n = %d ; umbral de la COTA (2g)^2 = %d"
          % (N, fact, g, (2 * g) ** 2))
    print("barrido: primos de %d a %d%s" % (arranque, PMAX,
                                            "" if PMIN is None else "   *** ARRANQUE FORZADO ***"))
    print("=" * 92)
    print("%-6s %-12s %-10s %-16s %s" % ("p", "pares (a,b)", "orbitas", "p/n!", "hay?"))
    primero = None
    for p in range(arranque, PMAX + 1):
        if not es_primo(p) or p <= N:
            continue
        pares = pares_escindidos(N, p)
        orb = pares // (p - 1) if pares else 0
        if orb and primero is None:
            primero = p
        # se imprimen los que tienen algo, y los diez primeros en cualquier caso
        if orb or p < arranque + 30:
            print("%-6d %-12d %-10d %-16.3f %s" % (p, pares, orb, p / fact, "SI" if orb else "no"))
            sys.stdout.flush()

    print("")
    if primero is None:
        print("  ninguna en el rango barrido; ampliar antes de concluir nada.")
    elif arranque > N + 1:
        print("  PRIMERO DE LA VENTANA: p = %d.  NO se afirma que sea el primero en absoluto:" % primero)
        print("  el barrido empezo en %d y no en %d.  Un minimo de una ventana no es un minimo."
              % (arranque, N + 1))
    else:
        print("  PRIMERA APARICION: p = %d,  que es %.3f veces n! = %d   (la cota exigia p > %d)"
              % (primero, primero / fact, fact, (2 * g) ** 2))
        print("  El barrido empezo en %d, o sea es una primera aparicion de verdad." % arranque)
    print("")

print("=" * 92)
print("QUE SE PUEDE DECIR Y QUE NO")
print("=" * 92)
print("  SI: la cota de Hasse-Weil es cuadraticamente pesimista, porque (2g_n)^2 ~ (n!/2)^2")
print("      mientras que el numero esperado de columnas extremas es p/n!.  Eso es aritmetica.")
print("  NO: que el fenomeno 'empieza' en p ~ n!.  p/n! >= 1 es cuando se espera UNA en promedio,")
print("      y la primera aparicion medida cae MUY por debajo --- una heuristica sobre el numero")
print("      esperado no fija el primer primo, y aqui el primero llega bastante antes.")
