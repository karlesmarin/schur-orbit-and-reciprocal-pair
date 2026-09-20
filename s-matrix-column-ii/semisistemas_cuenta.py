# -*- coding: utf-8 -*-
r"""semisistemas_cuenta.py -- cuantos semisistemas de F_p tienen suma nula, y la formula.

EL DATO.  Contando semisistemas centrados (semisistemas.py): p=7: 2, 11: 2, 13: 4, 17: 16,
19: 26, 23: 90.  Con m = (p-1)/2 y 2^m semisistemas en total, esos numeros son exactamente

        N(p) = ( 2^m + (2|p) * (p-1) ) / p ,        (2|p) = simbolo de Legendre de 2,

es decir, N(p) = (2^m - (2|p))/p + (2|p): la parte entera es el COCIENTE DE FERMAT de 2 a mitad de
exponente, que es entero por el criterio de Euler (2^m = (2|p) mod p).

POR QUE.  N(p) = (1/p) sum_{a mod p} prod_{t=1}^{m} (zeta^{at} + zeta^{-at}).  El termino a = 0 da
2^m.  Para a != 0 el producto recorre un semisistema y vale (2|p) por el lema de Gauss (el mismo
argumento que da el caracter cuadratico de 2).  De ahi los (p-1) terminos iguales.

Aqui se cuenta N(p) por programacion dinamica -- O(m p), no 2^m -- y se enfrenta a la formula.

    python semisistemas_cuenta.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo                                            # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 199


def cuenta(p):
    """semisistemas de suma 0 mod p, por DP sobre las parejas {t, -t}."""
    m = (p - 1) // 2
    d = [0] * p
    d[0] = 1
    for t in range(1, m + 1):
        nd = [0] * p
        for s in range(p):
            v = d[s]
            if v:
                nd[(s + t) % p] += v
                nd[(s - t) % p] += v
        d = nd
    return d[0]


def legendre(a, p):
    r = pow(a % p, (p - 1) // 2, p)
    return 0 if r == 0 else (1 if r == 1 else -1)


print("%-5s %-6s %-22s %-22s %s" % ("p", "(2|p)", "N(p) contado", "formula", "coinciden"))
fallos = 0
for p in range(5, PMAX + 1):
    if not es_primo(p):
        continue
    m = (p - 1) // 2
    N = cuenta(p)
    eps = legendre(2, p)
    F, resto = divmod(2 ** m + eps * (p - 1), p)
    ok = (resto == 0 and F == N)
    fallos += 0 if ok else 1
    if p <= 60 or not ok or p > PMAX - 20:
        print("%-5d %-6d %-22d %-22d %s" % (p, eps, N, F, "si" if ok else "NO"))
print("")
print("primos comprobados hasta %d ; fallos: %d" % (PMAX, fallos))
print("VEREDICTO:", "la formula SOBREVIVE" if fallos == 0 else "FALSA")
