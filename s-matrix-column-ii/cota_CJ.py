# -*- coding: utf-8 -*-
r"""cota_CJ.py -- ¿aguanta la cota explicita que propone la revision?

El Paso 2 del Teorema 5 solo da  C_J = binom(p,n)/p^{|J|} + O_n(p^{n/2}), sin constante.  La
revision propone una cota EXPLICITA, obtenida con la misma estrategia:

    | C_J(n,p) - p^{-|J|} binom(p,n) |  <=  (1 - p^{-|J|}) * (1/n!) * prod_{r=0}^{n-1} (2 sqrt p + r)

De donde sale: en una frecuencia no nula el polinomio de fase tiene grado entre 1 y 3, cada suma
de caracteres vale a lo sumo 2 sqrt p por Weil, y la identidad  sum_{sigma en S_n} x^{c(sigma)}
= x(x+1)...(x+n-1)  controla la suma de las contribuciones del tamiz por tipo de ciclo.

AQUI NO SE CONCEDE NADA: se coge el error MEDIDO de nuestros propios censos y se mira si cabe
dentro.  Una cota que no se comprueba es una cota que se cree.

Se mira ademas lo que de verdad importa para la nota: cuanto MARGEN sobra.  Una cota correcta
pero mil veces holgada no sirve para decir nada cuantitativo, y eso hay que decirlo si es el caso.

    python cota_CJ.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import io
import os
import sys
from math import comb, factorial, sqrt

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
N = 8


def cota(n, p, j):
    pr = 1.0
    for r in range(n):
        pr *= (2 * sqrt(p) + r)
    return (1 - p ** (-j)) * pr / factorial(n)


# los censos ya archivados: p, C_{1,2}, C_{1,3}, C_{1,2,3}
datos = []
texto = io.open(os.path.join(AQUI, "ley_2p_largo_OUT.txt"), encoding="utf-8",
                errors="replace").read()
for linea in texto.splitlines():
    f = linea.split()
    if len(f) == 7 and f[0].isdigit():
        datos.append((int(f[0]), int(f[1]), int(f[2]), int(f[3])))

print("n = %d.  Cota propuesta:  |C_J - binom(p,n)/p^|J||  <=  (1-p^-|J|) prod_{r<n}(2 sqrt p + r)/n!"
      % N)
print("")
print("%-5s %-9s %-14s %-14s %-10s %s"
      % ("p", "J", "error medido", "cota", "margen", "cabe?"))
casos = fallos = 0
margenes = []
for p, c12, c13, c123 in datos:
    tot = comb(p, N)
    for etiqueta, C, j in (("{1,2}", c12, 2), ("{1,3}", c13, 2), ("{1,2,3}", c123, 3)):
        err = abs(C - tot / p ** j)
        K = cota(N, p, j)
        ok = err <= K
        casos += 1
        fallos += 0 if ok else 1
        margenes.append(K / err if err else float("inf"))
        print("%-5d %-9s %-14.1f %-14.1f %-10.1f %s"
              % (p, etiqueta, err, K, K / err if err else float("inf"), "si" if ok else "NO"))
        sys.stdout.flush()

print("")
print("casos: %d ; fallos: %d" % (casos, fallos))
if margenes:
    m = sorted(margenes)
    print("margen de la cota sobre el error medido: minimo %.1f, mediana %.1f, maximo %.1f"
          % (m[0], m[len(m) // 2], m[-1]))
print("VEREDICTO:", "la cota explicita se sostiene en todo el rango medido" if not fallos
      else "LA COTA FALLA: no puede entrar en la nota")
print("")
print("LO QUE HAY QUE DECIR AL USARLA: es una cota de caso peor, no una estimacion del error.")
print("Sirve para hacer VERIFICABLE la dependencia en n --- que es lo que el O_n(p^{n/2}) esconde")
print("--- y para distinguir una cota de verdad de una envolvente ajustada al dato.  No sirve")
print("para predecir el tamano del error, y el margen medido lo deja claro.")
