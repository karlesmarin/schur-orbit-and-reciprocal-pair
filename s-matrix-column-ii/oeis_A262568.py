# -*- coding: utf-8 -*-
r"""oeis_A262568.py -- la forma cerrada contra la b-file REAL de la OEIS, con su offset.

POR QUE EXISTE.  La auditoria adversaria del 20-sep encontro que dos documentos nuestros decian
"comprobado contra la b-file" y NINGUN guion la abria: se habia mirado a ojo. Eso es una compuerta
declarada que no existe.  Aqui se descarga, se guarda y se compara, incluido el OFFSET --- que es
el riesgo que la comparacion consigo mismo no puede ver: si la sucesion empieza en n = 1 y nosotros
indexamos desde n = 3, la formula enviada estaria desplazada y ningun control interno lo notaria.

    python oeis_A262568.py            (usa la copia local si existe; si no, la descarga)

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
import urllib.request

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
COPIA = os.path.join(AQUI, "b262568.txt")
URL = "https://oeis.org/A262568/b262568.txt"


def jacobi(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def phi(n):
    r, x, d = n, n, 2
    while d * d <= x:
        if x % d == 0:
            while x % d == 0:
                x //= d
            r -= r // d
        d += 1
    if x > 1:
        r -= r // x
    return r


def cerrada(n):
    """nuestra formula, con m = 2n+1."""
    m = 2 * n + 1
    tot = sum(phi(d) * jacobi(2, d) * 2 ** ((m // d - 1) // 2) for d in range(1, m + 1) if m % d == 0)
    assert tot % m == 0, (n, m, tot)
    return tot // m


if not os.path.exists(COPIA):
    print("descargando %s ..." % URL)
    req = urllib.request.Request(URL, headers={"User-Agent": "curiosity-research/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        open(COPIA, "wb").write(r.read())
print("b-file local: %s (%d bytes)" % (COPIA, os.path.getsize(COPIA)))

pares = []
for linea in open(COPIA, encoding="utf-8"):
    linea = linea.strip()
    if not linea or linea.startswith("#"):
        continue
    a, b = linea.split()
    pares.append((int(a), int(b)))

print("terminos en la b-file: %d ; primer indice: %d ; ultimo: %d"
      % (len(pares), pares[0][0], pares[-1][0]))
print("OFFSET de la sucesion = %d ; a(%d) = %d" % (pares[0][0], pares[0][0], pares[0][1]))
print("")
print("%-6s %-6s %-24s %-24s %s" % ("n", "m=2n+1", "b-file", "nuestra formula", "coincide"))
fallos = 0
for n, v in pares:
    c = cerrada(n)
    ok = (c == v)
    fallos += 0 if ok else 1
    if n <= pares[0][0] + 6 or not ok or n >= pares[-1][0] - 2:
        print("%-6d %-6d %-24d %-24d %s" % (n, 2 * n + 1, v, c, "si" if ok else "NO"))

print("")
print("terminos comparados: %d ; discrepancias: %d" % (len(pares), fallos))
print("VEREDICTO:", "la formula reproduce la b-file ENTERA con el offset de la OEIS"
      if fallos == 0 else "HAY DESFASE O ERROR: no enviar nada")
