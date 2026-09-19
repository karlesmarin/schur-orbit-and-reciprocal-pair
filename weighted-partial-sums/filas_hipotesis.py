# -*- coding: utf-8 -*-
r"""filas_hipotesis.py -- las dos hipotesis que le faltan a la frase "In particular" de prop:filas.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

DOS ACUSACIONES DEL INFORME ADVERSARIO del 19-sep-2026 (F02 y F03), contra la FRASE EN PROSA que
sigue al recuadro.  El recuadro dice, correctamente, que para epsilon con valores +-1

    suma_{c<=m} c epsilon(c) = (n^2-1)/8  (mod 2),   luego solo se anula si (2/n) = 1,

y la frase anade: "no component T(m,chi) with chi quadratic vanishes unless 2 is a quadratic
residue modulo the length n".  Las dos acusaciones son sobre esa frase, no sobre el recuadro:

  F02  un caracter cuadratico NO toma valores +-1: vale 0 en los enteros no coprimos con su
       conductor.  Con un cero en el rango la identidad de paridad no aplica.
       Testigo propuesto: f=35, m=6, n=13.

  F03  (2/n) para n compuesto es el SIMBOLO DE JACOBI, que vale 1 sin que 2 sea un cuadrado.
       "quadratic residue modulo n" es por tanto una afirmacion distinta y mas fuerte.
       Testigo propuesto: m=7, n=15.

Se comprueban los dos testigos y se barre alrededor para saber si son aislados o sistematicos.
Salida: gates/capas_tipo/filas_hipotesis_OUT.txt
"""
import sys
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FALLOS = 0


def jacobi(a, n):
    """(a/n) de Jacobi, n impar positivo."""
    assert n % 2 == 1 and n > 0
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


def es_cuadrado(a, n):
    """2 es un cuadrado modulo n de verdad?  Se enumera; n es pequeno."""
    return any((x * x - a) % n == 0 for x in range(n))


def caracter_cuadratico(f):
    """chi(c) = (c/f) de Jacobi: el caracter cuadratico modulo f, que vale 0 si gcd(c,f)>1."""
    return lambda c: 0 if gcd(c, f) > 1 else jacobi(c, f)


def linea(txt):
    print("     " + txt)


print("=" * 92)
print("LAS DOS HIPOTESIS QUE LE FALTAN A LA FRASE 'In particular' DE prop:filas")
print("=" * 92)

# ---------------------------------------------------------------------------- F02
print("\n  F02  un caracter cuadratico vale CERO fuera de su conductor, no +-1")
print("       el testigo del informe: f=35, m=6, n=13")
f, m = 35, 6
n = 2 * m + 1
chi = caracter_cuadratico(f)
vals = [chi(c) for c in range(1, m + 1)]
S = sum(c * chi(c) for c in range(1, m + 1))
j = jacobi(2, n)
linea("chi(c) para c=1..6 : %s" % vals)
linea("suma c chi(c)      : %d" % S)
linea("(2/%d) de Jacobi    : %+d" % (n, j))
linea("chi es impar       : %s   (chi(-1) = %+d)" % (chi(f - 1) == -1, chi(f - 1)))
linea("chi es primitivo   : conductor 35 = 5*7, ambos componentes no triviales")
roto = (S == 0 and j == -1)
linea("LA FRASE DICE que S=0 obliga a (2/n)=+1.  Aqui S=0 y (2/n)=%+d  ->  %s"
      % (j, "REFUTADA" if roto else "no refutada"))
if not roto:
    FALLOS += 1
    linea("*** el testigo no reproduce: hay que mirarlo ***")

print("\n       barrido: cuantos pares (f,m) con un cero en el rango la rompen")
rotos = []
for f in range(5, 200, 2):
    chi = caracter_cuadratico(f)
    for m in range(2, 41):
        n = 2 * m + 1
        if all(gcd(c, f) == 1 for c in range(1, m + 1)):
            continue                      # sin ceros en el rango: la identidad si aplica
        S = sum(c * chi(c) for c in range(1, m + 1))
        if S == 0 and jacobi(2, n) == -1:
            rotos.append((f, m, n))
linea("pares (f,m) con f<200, m<=40 que rompen la frase: %d" % len(rotos))
linea("los seis primeros: %s" % rotos[:6])
if not rotos:
    FALLOS += 1
    linea("*** ninguno: la acusacion no se sostiene ***")

print("\n       CONTROL: sin ceros en el rango la frase tiene que aguantar siempre")
malos = 0
for f in range(5, 200, 2):
    chi = caracter_cuadratico(f)
    for m in range(2, 41):
        n = 2 * m + 1
        if any(gcd(c, f) > 1 for c in range(1, m + 1)):
            continue
        S = sum(c * chi(c) for c in range(1, m + 1))
        if S == 0 and jacobi(2, n) == -1:
            malos += 1
linea("contraejemplos con gcd(f, m!) = 1: %d  (tienen que ser 0)" % malos)
if malos:
    FALLOS += 1
    linea("*** entonces el fallo no es solo el cero, y la reparacion propuesta no basta ***")

# ---------------------------------------------------------------------------- F03
print("\n  F03  (2/n) de Jacobi no dice que 2 sea un cuadrado modulo n")
print("       el testigo del informe: m=7, n=15, epsilon = (1,1,-1,1,-1,-1,1)")
m, n = 7, 15
eps = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1, 7: 1}
mult = all(eps.get(a * b) == eps[a] * eps[b]
           for a in eps for b in eps if a * b in eps)
S = sum(c * eps[c] for c in range(1, m + 1))
linea("epsilon es multiplicativa en el intervalo : %s" % mult)
linea("suma c epsilon(c)                         : %d" % S)
linea("(2/15) de Jacobi                          : %+d" % jacobi(2, 15))
linea("2 es un cuadrado modulo 15                : %s  (los cuadrados son %s)"
      % (es_cuadrado(2, 15), sorted({x * x % 15 for x in range(15)})))
roto3 = (S == 0 and jacobi(2, 15) == 1 and not es_cuadrado(2, 15))
linea("LA FRASE DICE 'unless 2 is a quadratic residue mod n'  ->  %s"
      % ("REFUTADA" if roto3 else "no refutada"))
if not roto3:
    FALLOS += 1

print("\n       cuantos n impares compuestos tienen (2/n)=+1 sin que 2 sea cuadrado")
disc = [n for n in range(9, 200, 2) if jacobi(2, n) == 1 and not es_cuadrado(2, n)]
linea("n < 200 con (2/n)=+1 y 2 no cuadrado: %d   los seis primeros %s" % (len(disc), disc[:6]))
linea("en cambio, para n PRIMO las dos cosas coinciden siempre: %s"
      % all(jacobi(2, p) == (1 if es_cuadrado(2, p) else -1)
            for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)))
if not disc:
    FALLOS += 1

print("\n" + "=" * 92)
print("TOTAL DE FALLOS DE ESTA COMPUERTA: %d   (0 = las dos acusaciones se sostienen)" % FALLOS)
print("=" * 92)
sys.exit(1 if FALLOS else 0)
