# -*- coding: utf-8 -*-
r"""sato_tate.py -- la fluctuacion del recuento de columnas de defecto maximo es Sato-Tate.

LA CADENA.  Para n = 5, el numero de columnas con delta maximo esta dado (extremo_bring.py,
51 primos, 0 fallos) por

    N_max(5,p) = [ p + 1 - 4 a_p - 60 eps_R - 30 eps_4 - 24 eps_5 ] / 120,

con a_p la traza de la curva eliptica 50.a3.  Esa curva NO tiene multiplicacion compleja (su
j-invariante es -25/2, que no es uno de los trece j de CM), asi que el teorema de Sato-Tate
--- demostrado para curvas sin CM --- dice que los angulos theta_p definidos por
a_p = 2 sqrt(p) cos(theta_p) se equidistribuyen con densidad (2/pi) sin^2(theta).

CONSECUENCIA MEDIBLE.  El recuento de columnas de defecto maximo fluctua alrededor de p/120 con
una desviacion de tamano sqrt(p)/30 y la forma del semicirculo.  Esto no es una analogia: es un
test.  Aqui se calculan los theta_p para todos los primos de un rango y se comparan los momentos
y el histograma con la prediccion.

    python sato_tate.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import math
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo                                            # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5000


def ap(p):
    """traza de 50.a3 : y^2 + xy + y = x^3 - x - 2, por fuerza bruta."""
    a1, a2, a3, a4, a6 = 1, 0, 1, -1, -2
    n = 1
    cuad = [0] * p
    for y in range(p):
        cuad[(y * y) % p] = 1
    for x in range(p):
        c = (x * x * x + a2 * x * x + a4 * x + a6) % p
        b = (a1 * x + a3) % p
        d = (b * b + 4 * c) % p
        if d == 0:
            n += 1
        elif cuad[d]:
            n += 2
    return p + 1 - n


thetas, hasse = [], 0
for p in range(7, PMAX + 1):
    if not es_primo(p) or p in (2, 5):
        continue
    a = ap(p)
    x = a / (2 * math.sqrt(p))
    if abs(x) > 1:
        hasse += 1
        continue
    thetas.append(math.acos(x))

n = len(thetas)
print("primos usados: %d  (hasta %d) ; violaciones de Hasse: %d" % (n, PMAX, hasse))

# momentos de a_p/(2 sqrt p) = cos(theta): los de Sato-Tate son los numeros de Catalan alternados:
# E[cos^{2k}] = C_k / 4^k  con C_k el k-esimo Catalan; E[cos^{impar}] = 0.
def catalan(k):
    c = 1
    for i in range(k):
        c = c * 2 * (2 * i + 1) // (i + 2)
    return c


print("")
print("%-6s %-14s %-14s %s" % ("k", "E[x^k] medido", "Sato-Tate", "diferencia"))
for k in range(1, 9):
    med = sum(math.cos(t) ** k for t in thetas) / n
    if k % 2:
        pred = 0.0
    else:
        pred = catalan(k // 2) / (4 ** (k // 2))
    print("%-6d %-14.6f %-14.6f %+.6f" % (k, med, pred, med - pred))

print("")
print("histograma de theta en 6 cajas, contra (2/pi) sin^2:")
cajas = 6
for i in range(cajas):
    a0, a1_ = i * math.pi / cajas, (i + 1) * math.pi / cajas
    obs = sum(1 for t in thetas if a0 <= t < a1_) / n
    # integral de (2/pi) sin^2 entre a0 y a1
    F = lambda t: (t - math.sin(2 * t) / 2) / math.pi
    pred = F(a1_) - F(a0)
    barra = "#" * int(round(obs * 120))
    print("   [%.2f, %.2f)  obs %.4f  pred %.4f  %s" % (a0, a1_, obs, pred, barra))
