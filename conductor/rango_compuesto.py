# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# Observacion 12 (q' compuesto) de la nota del conductor:
# para n = q' >= 3 (n no = 2 mod 4), p primo impar con p no | n, y S(c) = 2c - n (N = 1) sobre
# TODOS los residuos c in [1, n-1] (tambien no unidades), compara
#     rango_{F_p} [S(u^{-1} c)]_{u unidad, c residuo} == phi(n)/2      con      p no | h^-(Q(zeta_n)).
# h^- se calcula con Washington, Thm 4.17: h^- = Q w prod_{chi impar} (-1/2 B_{1,chi*}), chi* primitivo,
# Q = 1 si n es potencia de primo y 2 si no, w = 2n (n impar) o n (n = 0 mod 4).  Se contrasta con
# valores conocidos (n primo de bernoulli_W1 y 39 -> 2, 65 -> 64).
#
# Uso: python rango_compuesto.py [NMAX PMAX]
import cmath
import os
import sys
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bernoulli_W1 import rango_mod_p, es_primo  # noqa: E402


def factoriza(n):
    out, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def caracteres(n):
    """Todos los caracteres de (Z/n)^* como dict a -> complejo."""
    units = [a for a in range(1, n) if gcd(a, n) == 1]
    # generadores por componentes primarias y CRT
    comps = []
    for p, k in factoriza(n).items():
        pk = p ** k
        grupo = [a for a in range(1, pk) if gcd(a, pk) == 1]
        if p == 2 and k >= 3:
            gens = [(pk - 1, 2), (5, 2 ** (k - 2))]
        else:
            orden = len(grupo)
            g = next(a for a in grupo if all(pow(a, orden // r, pk) != 1 for r in factoriza(orden)))
            gens = [(g, orden)]
        for g, o in gens:
            comps.append((pk, g, o))
    # logaritmo discreto de cada unidad en la base de generadores (fuerza bruta, n pequeno)
    import itertools
    rangos = [range(o) for (_, _, o) in comps]
    log = {}
    for exps in itertools.product(*rangos):
        # construir el residuo por CRT
        x = 1
        for (pk, g, o), e in zip(comps, exps):
            # elemento que es g^e mod pk y 1 mod n/pk
            y = pow(g, e, pk)
            m = n // pk
            # CRT
            t = (y - 1) * pow(m, -1, pk) % pk
            z = (1 + m * t) % n
            x = x * z % n
        log[x] = exps
    assert len(log) == len(units)
    chars = []
    for js in itertools.product(*rangos):
        ch = {a: cmath.exp(2j * cmath.pi * sum(j * e / o for j, e, (_, _, o) in zip(js, log[a], comps)))
              for a in units}
        chars.append(ch)
    return chars


def conductor(ch, n):
    for d in sorted(x for x in range(1, n + 1) if n % x == 0):
        if all(abs(ch[a] - 1) < 1e-9 for a in ch if a % d == 1 % d):
            return d
    return n


def B1_primitivo(ch, n):
    f = conductor(ch, n)
    if f == 1:
        return None
    tot = 0
    for a in range(1, f + 1):
        if gcd(a, f) != 1:
            continue
        # levantar a a unidad mod n con a' = a mod f
        lift = next(b for b in range(a, n * f + 1, f) if gcd(b, n) == 1) % n
        tot += ch[lift] * a
    return tot / f


def h_menos(n):
    chars = caracteres(n)
    prod = 1
    for ch in chars:
        if abs(ch[n - 1] + 1) > 1e-9:
            continue  # par
        B = B1_primitivo(ch, n)
        prod *= -B / 2
    Q = 1 if len(factoriza(n)) == 1 else 2
    w = 2 * n if n % 2 else n
    val = Q * w * prod
    assert abs(val.imag) < 1e-6 * max(1, abs(val.real)), (n, val)
    return round(val.real)


def rango_todas(n, p):
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    filas = [[(2 * ((u * c) % n) - n) % p for c in range(1, n)] for u in units]
    return rango_mod_p(filas, p)


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 90
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    print("control de h^-:", {n: h_menos(n) for n in (23, 31, 39, 41, 65)}, " (esperado 23:3 31:9 39:2 41:121 65:64)")
    pares = coinciden = 0
    fallos = []
    positivos = 0
    for n in range(3, NMAX):
        if n % 4 == 2:
            continue
        full = sum(1 for u in range(1, n) if gcd(u, n) == 1) // 2
        h = h_menos(n)
        for p in range(3, PMAX + 1):
            if not es_primo(p) or n % p == 0:
                continue
            pares += 1
            total = rango_todas(n, p) == full
            pred = h % p != 0
            positivos += (not pred)
            if total == pred:
                coinciden += 1
            else:
                fallos.append((n, p, h))
    print("n in [3,%d), n no = 2 mod 4;  p impar <= %d, p no | n:  pares %d   coinciden %d   (p | h^-: %d)   fallos %s"
          % (NMAX, PMAX, pares, coinciden, positivos, fallos[:10]))


if __name__ == "__main__":
    main()
