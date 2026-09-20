# -*- coding: utf-8 -*-
r"""censo_delta.py -- la tabla de Delta(n, p) = suma de delta(T) sobre las orbitas
multiplicativas de espectros centrados T subset F_p con |T| = n.

POR QUE.  La formula por orbitas del pegado del anillo de fusion,
    G(n, p) = (n/2)[ D(a-1) + b - 2*Delta ],   D = C(p-1, n-1)/n,  a = (n-1)(p-n-1)/(p-1),
con b = numero de orbitas, deja UN solo termino sin forma cerrada: Delta = sum_orbitas delta(T).
Todo lo demas es combinatoria elemental.  Esto mide Delta para buscarle la forma.

Se usa el instrumento publicado de P3a (release_tool/P3a.py): `centrados` enumera un representante
por orbita y `s_exacto` da el semigrupo, cuyo genero es delta (Teorema del semigrupo exacto, altura
prima).  Para p <= 31 y n <= 8 el calculo es de segundos.

    python censo_delta.py [PMAX] [NMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys
from fractions import Fraction
from math import comb

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import centrados, es_primo, semigrupo, s_exacto            # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 31
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8


def delta_de(T, p):
    """delta(T) = genero de S_T; s_exacto devuelve los generadores."""
    return semigrupo(s_exacto(T, p))["genero"]


def fila(n, p):
    """Delta(n, p), el numero de orbitas b, y el reparto de delta por valor."""
    orbs = centrados(p, n)
    reparto, total = {}, 0
    for T in orbs:
        d = delta_de(T, p)
        reparto[d] = reparto.get(d, 0) + 1
        total += d
    return total, len(orbs), reparto


print("=" * 96)
print("Delta(n, p) = suma de delta sobre orbitas de espectros centrados,  p <= %d, n <= %d" % (PMAX, NMAX))
print("=" * 96)
for p in range(5, PMAX + 1):
    if not es_primo(p):
        continue
    filas = []
    for n in range(2, min(NMAX, p - 2) + 1):
        D, b, rep = fila(n, p)
        filas.append((n, D, b, rep))
    print("p = %d" % p)
    for n, D, b, rep in filas:
        C = comb(p - 1, n - 1)
        print("   n=%-2d  Delta=%-6d  orbitas=%-5d  columnas=C(p-1,n-1)=%-6d  Delta/orbitas=%-8s  reparto=%s"
              % (n, D, b, C, str(Fraction(D, b)) if b else "-",
                 " ".join("%d:%d" % (k, rep[k]) for k in sorted(rep))))
    sys.stdout.flush()
