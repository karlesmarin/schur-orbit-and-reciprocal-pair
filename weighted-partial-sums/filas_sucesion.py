#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""filas_sucesion.py -- que fila es cada valor de la sucesion 1,1,1,1,2,2,2,4,1,9,5,14.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

La nota dice, entre las fronteras: "el recuento de soluciones por fila, cuyos primeros valores
1,1,1,1,2,2,2,4,1,9,5,14 no sabemos situar".  NO DICE DE QUE FILAS.  Si un lector supone que son
las filas m = 1, 2, 3, ..., la sucesion es falsa: a(1) = 0, porque la suma vale 1 y no cero.

Lo que a(m) cuenta: las funciones eps de {1..m} en {+-1}, COMPLETAMENTE MULTIPLICATIVAS -- luego
determinadas por sus valores en los primos <= m -- con sum_{c<=m} c eps(c) = 0.

Uso:  python filas_sucesion.py          (escribe filas_sucesion_OUT.txt)
"""
import io
import itertools
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SAL = []


def di(t=""):
    print(t, flush=True)
    SAL.append(t)


def primos(m):
    return [p for p in range(2, m + 1)
            if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def a(m):
    """Cuantas eps multiplicativas anulan la suma."""
    ps = primos(m)
    # eps(c) = producto de eps(p)^{v_p(c)}: solo importa la paridad del exponente
    par = []
    for c in range(1, m + 1):
        e = []
        x = c
        for i, p in enumerate(ps):
            k = 0
            while x % p == 0:
                x //= p
                k += 1
            if k % 2:
                e.append(i)
        par.append((c, tuple(e)))
    n = 0
    sols = []
    for x in itertools.product((1, -1), repeat=len(ps)):
        s = 0
        for c, e in par:
            sg = 1
            for i in e:
                sg *= x[i]
            s += c * sg
        if s == 0:
            n += 1
            if len(sols) < 2:
                sols.append(dict(zip(ps, x)))
    return n, sols


def main():
    di("=" * 92)
    di("QUE FILA ES CADA VALOR de 1,1,1,1,2,2,2,4,1,9,5,14")
    di("=" * 92)
    di("%-5s %-6s %-14s %-9s %s" % ("m", "n=2m+1", "(n^2-1)/8 par?", "(2|n)", "a(m)"))
    vals = []
    for m in range(1, 29):
        n = 2 * m + 1
        t = (n * n - 1) // 8
        cnt, sols = a(m)
        vals.append((m, cnt))
        di("%-5d %-6d %-14s %-9s %d%s" % (
            m, n, "si" if t % 2 == 0 else "NO (impar)",
            "+1" if t % 2 == 0 else "-1", cnt,
            "   p.ej. " + str(sols[0]) if cnt and m <= 8 else ""))

    di("")
    nz = [(m, c) for m, c in vals if c]
    di("  filas NO VACIAS y su valor:")
    di("    m   = " + ", ".join(str(m) for m, _ in nz))
    di("    a(m)= " + ", ".join(str(c) for _, c in nz))
    di("")
    PUB = [1, 1, 1, 1, 2, 2, 2, 4, 1, 9, 5, 14]
    obt = [c for _, c in nz][:len(PUB)]
    coincide = (obt == PUB)
    di("  la sucesion publicada 1,1,1,1,2,2,2,4,1,9,5,14 es:")
    di("    la de TODAS las filas m=1,2,3,...      -> %s" %
       ("SI" if [c for _, c in vals][:len(PUB)] == PUB else "NO"))
    di("    la de las filas NO VACIAS, en orden    -> %s" % ("SI" if coincide else "NO"))
    di("")
    if coincide:
        di("  Es decir: son las filas m = %s y no m = 1,2,3,..." %
           ", ".join(str(m) for m, _ in nz[:len(PUB)]))
        di("  Esas m son exactamente las que cumplen m = 0 o 3 (mod 4), que es donde")
        di("  (n^2-1)/8 es par -- la Proposicion de las filas.")
    di("")
    di("=" * 92)
    di("VEREDICTO: la nota %s decir de que filas habla." %
       ("TIENE que" if coincide else "dice mal la sucesion y hay que"))
    di("=" * 92)
    io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "filas_sucesion_OUT.txt"), "w",
            encoding="utf-8", newline="\n").write("\n".join(SAL) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
