# -*- coding: utf-8 -*-
r"""soporte_por_momentos.py -- la conjetura del soporte reducida a momentos: barrido aritmetico.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Reduccion (tangente_momentos.py: 2974 celdas, 0 fallos): con W_0 = d y E >= 2, A es maximal en p
sii existe una clase c (con 2c != 0 si el segmento es simetrico) con p no | S(c).  conj:sop sigue si,
para todo segmento del papel FUERA de las familias (r no divide n-1, n, n+1) y no estable, alguna
S(c) es no nula mod p.  Esto es pura aritmetica: aqui se barre sin capas.
Segmentos del papel (convencion del texto, lineas 261-274 y 611): SU(n) {n-1, n-3, ..., 1-n} mod
M = 2q (n par) o q (n impar, exponentes divididos por 2, conjugado de Galois); Sp(2m) {+-1..+-m}
mod 2q? -- aqui SU y GL, que son los dos casos del barrido del papel.
Uso: python soporte_por_momentos.py QMAX
"""
import sys
from math import gcd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
# CONTROL: con "familias" se barre DENTRO de las familias, donde el indice si es divisible por p:
# alli la prueba tiene que encontrar casos con todas las S(c) nulas, o no ve nada.
DENTRO = len(sys.argv) > 2 and sys.argv[2] == "familias"


def fact(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def vp(n, p):
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v


def segmento(tipo, n, q):
    if tipo == "GL":
        return list(range(n)), q
    if n % 2 == 0:
        return [n + 1 - 2 * j for j in range(1, n + 1)], 2 * q
    return [((n + 1 - 2 * j) // 2) for j in range(1, n + 1)], q  # (n+1-2j)/2 entero para n impar


tot = 0; malos = []
for q in range(6, QMAX + 1):
    for p in fact(q):
        if p == 2:
            continue
        for tipo in ("SU", "GL"):
            for n in range(2, q):
                a, M = segmento(tipo, n, q)
                r = q // p ** vp(q, p)
                qp = M // p ** vp(M, p)
                en_familia = any(x % r == 0 for x in (n - 1, n, n + 1))
                if r < 3 or en_familia != DENTRO:
                    continue
                sim = sorted(x % M for x in a) == sorted((-x) % M for x in a)
                # estabilidad en q': el multiconjunto de clases invariante por las unidades
                cls = sorted(x % qp for x in a)
                estable = all(sorted(u * x % qp for x in cls) == cls for u in range(1, qp) if gcd(u, qp) == 1)
                if estable:
                    continue
                Sc = {}
                for x in a:
                    Sc[x % qp] = Sc.get(x % qp, 0) + x
                ok = any(Sc[c] % p for c in Sc if not (sim and (2 * c) % qp == 0))
                tot += 1
                if not ok:
                    malos.append((tipo, n, q, p, dict(sorted(Sc.items()))))
print("segmentos fuera de las familias, inestables, q <= %d, p impar: %d" % (QMAX, tot))
print("con TODAS las S(c) divisibles por p (contraejemplo potencial): %d" % len(malos))
for m in malos[:10]:
    print("  ", m)
