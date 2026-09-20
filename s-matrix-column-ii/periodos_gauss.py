# -*- coding: utf-8 -*-
r"""periodos_gauss.py -- los semisistemas H-invariantes se cuentan con PERIODOS DE GAUSS.

LA CADENA.  Sea H un subgrupo de orden IMPAR h de (Z/m)^x (impar porque -1 no puede estar en el
estabilizador de un semisistema).  Un semisistema T invariante por H es una union de H-orbitas que
toma exactamente una de cada pareja {O, -O}.  La condicion de suma nula se vuelve

        suma de los sigma(O) elegidos = 0 (mod m),      sigma(O) = suma de los elementos de O,

y esos sigma(O) son exactamente los PERIODOS DE GAUSS de orden h modulo m.  Es decir:

  #{semisistemas H-invariantes de suma nula} = #{ elecciones de signo sobre los periodos de Gauss
                                                  de orden h cuya suma con signo es 0 mod m }.

El caso h = 1 es el problema original (los periodos son los propios elementos).  O sea que el
reticulo de orbitas del recuento se organiza como el MISMO problema, un piso mas arriba, con los
periodos de Gauss en lugar de los residuos.  Los periodos de Gauss son ademas los generadores del
orden C_T de P3a: los dos hilos se tocan aqui.

ESTE GUION comprueba la identificacion contando de las dos maneras.

    python periodos_gauss.py [MMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import sys
from itertools import product
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 31


def es_primo(x):
    return x > 1 and all(x % d for d in range(2, int(x ** 0.5) + 1))


def subgrupo(g, h, p):
    """el subgrupo de orden h generado por g^{(p-1)/h}."""
    b = pow(g, (p - 1) // h, p)
    S, x = set(), 1
    for _ in range(h):
        S.add(x)
        x = x * b % p
    return sorted(S)


def raiz_primitiva(p):
    for g in range(2, p):
        vistos, x = set(), 1
        for _ in range(p - 1):
            x = x * g % p
            vistos.add(x)
        if len(vistos) == p - 1:
            return g
    raise ValueError


print("%-5s %-4s %-30s %-16s %-16s %s"
      % ("p", "h", "periodos de Gauss (sumas)", "invariantes", "por periodos", "coincide"))
fallos = casos = 0
for p in range(7, MMAX + 1, 2):
    if not es_primo(p):
        continue
    n = (p - 1) // 2
    g = raiz_primitiva(p)
    impares = [h for h in range(1, p) if (p - 1) % h == 0 and h % 2 == 1]
    for h in impares:
        H = subgrupo(g, h, p)
        if -1 % p in H or (p - 1) in H:
            continue                               # -1 en H: no hay semisistemas invariantes
        # las H-orbitas y sus periodos
        orbitas, vistos = [], set()
        for x in range(1, p):
            if x in vistos:
                continue
            O = sorted(x * a % p for a in H)
            vistos |= set(O)
            orbitas.append(O)
        # emparejar O con -O
        pares, usados = [], set()
        for i, O in enumerate(orbitas):
            if i in usados:
                continue
            mO = sorted((p - y) % p for y in O)
            j = next(k for k, Q in enumerate(orbitas) if Q == mO)
            usados |= {i, j}
            pares.append((O, orbitas[j]))
        # (a) recuento directo: semisistemas H-invariantes de suma nula
        directo = 0
        for elec in product(*[(A, B) for A, B in pares]):
            T = [x for O in elec for x in O]
            if len(set(T)) == n and sum(T) % p == 0:
                directo += 1
        # (b) por periodos: elecciones de signo sobre los periodos de Gauss
        sigmas = [sum(A) % p for A, B in pares]
        porper = sum(1 for e in product([1, -1], repeat=len(pares))
                     if sum(s * ei for s, ei in zip(sigmas, e)) % p == 0)
        casos += 1
        ok = (directo == porper)
        fallos += 0 if ok else 1
        print("%-5d %-4d %-30s %-16d %-16d %s"
              % (p, h, str(sigmas)[:29], directo, porper, "si" if ok else "NO"))
        sys.stdout.flush()

print("")
print("casos: %d ; fallos: %d" % (casos, fallos))
print("VEREDICTO:", "los semisistemas H-invariantes SON las elecciones de signo sobre periodos de Gauss"
      if not fallos else "la identificacion FALLA")

print("")
print("=" * 96)
print("EL REMATE: para h > 1 TODOS los periodos son cero, y el recuento se cierra")
print("=" * 96)
print("La suma de un subgrupo H de orden h > 1 de F_p^x es 0: H es el conjunto de raices de")
print("X^h - 1, cuya suma es el coeficiente de X^{h-1}, que es 0 en cuanto h >= 2.  Luego cada")
print("periodo sigma(aH) = a * sum(H) = 0 y TODA eleccion de signos cumple la condicion:")
print("")
print("        A(h) := #{semisistemas invariantes por el subgrupo de orden h} = 2^{(p-1)/(2h)},")
print("")
print("para todo h > 1 impar, mientras que A(1) = N(p).  Con eso el RETICULO DE ORBITAS queda")
print("determinado por inversion de Mobius sobre los divisores impares de p-1.")
print("")
print("%-5s %-5s %-14s %-16s %-16s %s" % ("p", "h", "A(h) formula", "estab. exacto h", "orbitas", "medido"))
fallos2 = casos2 = 0
for p in range(7, MMAX + 1, 2):
    if not es_primo(p):
        continue
    n = (p - 1) // 2
    g = raiz_primitiva(p)
    impares = [h for h in range(1, p) if (p - 1) % h == 0 and h % 2 == 1]
    # medicion directa del reticulo
    U = list(range(1, p))
    sols = []
    for elec in product(*[(t, p - t) for t in range(1, n + 1)]):
        T = frozenset(elec)
        if sum(T) % p == 0:
            sols.append(T)
    conj = set(sols)
    medido = {}
    for T in sols:
        h = sum(1 for a in U if frozenset(a * t % p for t in T) == T)
        medido[h] = medido.get(h, 0) + 1
    for h in impares:
        A = (2 ** ((p - 1) // (2 * h))) if h > 1 else len(sols)
        # estabilizador exactamente h: Mobius sobre los multiplos impares de h que dividen p-1
        exacto = 0
        for hp in impares:
            if hp % h == 0:
                k = hp // h
                mu = 1
                x, primos = k, set()
                d = 2
                while d * d <= x:
                    if x % d == 0:
                        c = 0
                        while x % d == 0:
                            x //= d
                            c += 1
                        if c > 1:
                            mu = 0
                        primos.add(d)
                    d += 1
                if x > 1:
                    primos.add(x)
                if mu:
                    mu = (-1) ** len(primos)
                Ahp = (2 ** ((p - 1) // (2 * hp))) if hp > 1 else len(sols)
                exacto += mu * Ahp
        orb = exacto // ((p - 1) // h) if exacto else 0
        m_ = medido.get(h, 0)
        ok = (exacto == m_)
        casos2 += 1
        fallos2 += 0 if ok else 1
        print("%-5d %-5d %-14d %-16d %-16s %-8d %s"
              % (p, h, A, exacto, str(orb), m_, "si" if ok else "NO"))
print("")
print("casos: %d ; fallos: %d" % (casos2, fallos2))
print("VEREDICTO:", "el reticulo de orbitas queda descrito en forma cerrada"
      if not fallos2 else "la descripcion FALLA")
