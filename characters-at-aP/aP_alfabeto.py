# -*- coding: utf-8 -*-
# EL ALFABETO DE a_P^k ES EL DE PRASAD 2016 MENOS DOS PUNTOS.   25 de agosto de 2026.
#
# POR QUE.  La nota decia que los dos enunciados "tienen la misma forma" y que no habiamos derivado
# uno del otro.  Pero si la forma coincide entrada por entrada, lo natural es que coincidan los
# OBJETOS.  Esta gate lo comprueba: el alfabeto en el que se evalua el caracter simplectico en
# a_P^k es, como MULTICONJUNTO, el alfabeto del Teorema 2 de Prasad con (m,n) = (d,q), quitando
# exactamente dos elementos.
#
# EL ENUNCIADO QUE SE MIDE.  Sp(2m), t = 2m+2, zeta de orden t, d = gcd(k,t), q = t/d.  El alfabeto
# simplectico de a_P^k es  A = { zeta^{+- j k} : j = 1..m }.  Entonces:
#
#      d PAR    ->   A = d copias de mu_q,  menos DOS copias de  +1
#      d IMPAR  ->   A = d copias de mu_q,  menos una de +1 y una de -1
#
# LA RAZON, que es de tres lineas y va escrita antes de medir:
#   {+-j : j=1..m} = Z/t \ {0, t/2}, porque t = 2m+2.
#   zeta^k tiene orden q, luego zeta^{k e} solo depende de e mod q, y cada residuo mod q lo alcanzan
#   exactamente d valores de e en Z/t.  Asi que { zeta^{k e} : e in Z/t } = d copias de mu_q.
#   Falta quitar los dos e excluidos:  e = 0    -> residuo 0     -> el elemento +1;
#                                      e = t/2  -> t/2 = dq/2, que mod q vale 0 si d es par
#                                                  (elemento +1) y q/2 si d es impar (elemento -1).
#
# CONSECUENCIA, y es lo que la nota necesitaba: el alfabeto de Prasad en su Teorema 2 con t=1 es
# "mu_n con multiplicidad m".  El nuestro es mu_q con multiplicidad d, menos dos puntos.  Luego su
# m ES nuestro d y su n ES nuestro q --- no por parecido, por identidad de alfabetos --- y los dos
# puntos que sobran son exactamente las dos CLASES FIJAS del plegado.
#
# QUE SE MIDE
#   A0  CONTROL: el multiconjunto A calculado a pelo tiene 2m elementos y es cerrado por inverso.
#   A1  la identidad de multiconjuntos, para todo (m,k) del rango.
#   A2  CONTROL: el enunciado con la paridad de d CAMBIADA debe fallar.  Si no fallara, A1 no
#       estaria midiendo la paridad.
#   A3  la consecuencia sobre H(z):  d par -> (1-z)^2/(1-z^q)^d ;  d impar -> (1-z^2)/(1-z^q)^d.
#       Se contrasta contra la serie calculada desde los autovalores.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python aP_alfabeto.py > aP_alfabeto_OUT.txt 2>&1

import cmath
import sys
from collections import Counter
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []


def ok(cond, etiqueta, detalle=""):
    print(("   OK   " if cond else "  FALLA ") + etiqueta + ("   " + detalle if detalle else ""))
    if not cond:
        fallos.append(etiqueta)


def alfabeto(m, k):
    """A = { zeta^{+- j k} }, como multiconjunto de EXPONENTES mod t."""
    t = 2 * m + 2
    return Counter([(j * k) % t for j in range(1, m + 1)] +
                   [(-j * k) % t for j in range(1, m + 1)])


def como_residuos_mod_q(expon, t, q):
    """zeta^e con zeta de orden t; zeta^k tiene orden q.  El elemento es zeta_q^{e/d * ...}:
    lo que identifica al elemento es e mod t, y su imagen en mu_q es e*(q/t)... se hace numerico."""
    return None


print(SEP)
print("A0 -- CONTROL: el alfabeto tiene 2m elementos y es cerrado por inverso")
for (m, k) in [(3, 1), (5, 2), (7, 4), (6, 3)]:
    t = 2 * m + 2
    A = alfabeto(m, k)
    cerrado = all(A[(-e) % t] == c for e, c in A.items())
    ok(sum(A.values()) == 2 * m and cerrado, f"m={m} k={k}",
       f"{sum(A.values())} elementos, cerrado por inverso: {cerrado}")

print(SEP)
print("A1 -- LA IDENTIDAD DE MULTICONJUNTOS")
print("     d par   ->  A = d*mu_q  -  2*{+1}")
print("     d impar ->  A = d*mu_q  -  {+1} - {-1}")
print()
print(f"{'m':>3} {'k':>3} {'t':>4} {'d':>3} {'q':>4} {'paridad':>8} {'coincide':>9}   {'diferencia si no'}")
for m in range(2, 16):
    t = 2 * m + 2
    for k in range(1, 13):
        d = gcd(k, t)
        q = t // d
        if q < 2:
            continue
        # el alfabeto, como multiconjunto de elementos de mu_q identificados por su exponente mod q
        z = cmath.exp(2j * cmath.pi / t)
        obs = Counter()
        for e in alfabeto(m, k).elements():
            w = z ** e                                   # elemento del circulo
            # su exponente en mu_q:  w = exp(2 pi i r / q)
            r = round(cmath.phase(w) / (2 * cmath.pi) * q) % q
            assert abs(w - cmath.exp(2j * cmath.pi * r / q)) < 1e-9, (m, k, e)
            obs[r] += 1
        pred = Counter({r: d for r in range(q)})
        if d % 2 == 0:
            pred[0] -= 2
        else:
            pred[0] -= 1
            pred[q // 2] -= 1                            # q par garantizado si d impar
        pred = Counter({r: c for r, c in pred.items() if c})
        igual = obs == pred
        dif = "" if igual else f"obs-pred={sorted((obs-pred).items())} pred-obs={sorted((pred-obs).items())}"
        if not igual or (m <= 5 and k <= 4):
            print(f"{m:>3} {k:>3} {t:>4} {d:>3} {q:>4} {'par' if d%2==0 else 'impar':>8} "
                  f"{'SI' if igual else 'NO':>9}   {dif}")
        if not igual:
            fallos.append(f"A1 m={m} k={k}")
ok(not [f for f in fallos if f.startswith("A1")], "A1: la identidad se cumple en todo el rango barrido")

print(SEP)
print("A2 -- CONTROL: con la paridad de d CAMBIADA, el enunciado tiene que FALLAR")
mal = tot = 0
for m in range(2, 12):
    t = 2 * m + 2
    for k in range(1, 9):
        d = gcd(k, t)
        q = t // d
        if q < 3:
            continue
        z = cmath.exp(2j * cmath.pi / t)
        obs = Counter()
        for e in alfabeto(m, k).elements():
            r = round(cmath.phase(z ** e) / (2 * cmath.pi) * q) % q
            obs[r] += 1
        malo = Counter({r: d for r in range(q)})
        if d % 2 == 1:                                   # <- A PROPOSITO al reves
            malo[0] -= 2
        else:
            malo[0] -= 1
            malo[q // 2] -= 1
        malo = Counter({r: c for r, c in malo.items() if c})
        tot += 1
        mal += (obs == malo)
ok(mal == 0, "el enunciado con la paridad invertida no acierta nunca",
   f"{mal} aciertos de {tot} --- si fuera >0, A1 no estaria midiendo la paridad")

print(SEP)
print("A3 -- LA CONSECUENCIA SOBRE H(z)")
print("     d par:  H(z) = (1-z)^2 / (1-z^q)^d        d impar:  H(z) = (1-z^2) / (1-z^q)^d")
N = 26


def h_serie(m, k):
    t = 2 * m + 2
    z = cmath.exp(2j * cmath.pi / t)
    h = [0j] * (N + 1)
    h[0] = 1 + 0j
    for e in alfabeto(m, k).elements():
        y = z ** e
        acc, nue = 0j, [0j] * (N + 1)
        for n in range(N + 1):
            acc = acc * y + h[n]
            nue[n] = acc
        h = nue
    return h


def h_cerrada(d, q):
    from math import comb
    c = [0] * (N + 1)
    for j in range(N // q + 1):
        c[q * j] = comb(j + d - 1, d - 1)
    out = [0] * (N + 1)
    for n in range(N + 1):
        if d % 2 == 0:                                   # (1-z)^2 = 1 - 2z + z^2
            out[n] = c[n] - 2 * (c[n - 1] if n >= 1 else 0) + (c[n - 2] if n >= 2 else 0)
        else:                                            # (1-z^2)
            out[n] = c[n] - (c[n - 2] if n >= 2 else 0)
    return out


print(f"{'m':>3} {'k':>3} {'d':>3} {'q':>4} {'max|serie-cerrada|':>20} {'h_0..h_6':>28}")
for (m, k) in [(2, 1), (3, 2), (4, 1), (5, 2), (5, 3), (5, 4), (7, 4), (6, 3), (9, 5), (11, 6)]:
    t = 2 * m + 2
    d = gcd(k, t)
    q = t // d
    if q < 2:
        continue
    hs, hc = h_serie(m, k), h_cerrada(d, q)
    peor = max(abs(hs[n] - hc[n]) for n in range(N + 1))
    ok_fila = peor < 1e-7
    if not ok_fila:
        fallos.append(f"A3 m={m} k={k}")
    print(f"{m:>3} {k:>3} {d:>3} {q:>4} {peor:>20.2e} {str(hc[:7]):>28}"
          + ("" if ok_fila else "   <-- FALLA"))
ok(not [f for f in fallos if f.startswith("A3")], "A3: la forma cerrada reproduce la serie")

print(SEP)
print("A4 -- LO QUE ESO SIGNIFICA: es EL MISMO ELEMENTO, visto en Sp(t) y en Sp(t-2)")
print("     d*mu_q tiene dq = t elementos y es cerrado por inverso, luego es el multiconjunto de")
print("     autovalores de un elemento de Sp(t) --- y ese elemento es el del Teorema 2 de Prasad")
print("     con (m,n) = (d,q) evaluado en t=1.  Nuestro a_P^k es ESE MISMO con un plano")
print("     simplectico quitado:  {1,1} si d es par,  {1,-1} si d es impar.")
print()
print(f"{'m':>3} {'k':>3} {'t':>4} {'d':>3} {'q':>4} {'|d*mu_q|':>9} {'= t':>5} "
      f"{'inv-cerrado':>12} {'plano':>9}")
malA4 = 0
for m in range(2, 14):
    t = 2 * m + 2
    for k in range(1, 9):
        d = gcd(k, t)
        q = t // d
        if q < 2:
            continue
        big = Counter({r: d for r in range(q)})
        n_big = sum(big.values())
        inv = all(big[(-r) % q] == c for r, c in big.items())
        bien = (n_big == t) and inv
        malA4 += (not bien)
        if m <= 6 and k <= 4:
            print(f"{m:>3} {k:>3} {t:>4} {d:>3} {q:>4} {n_big:>9} {str(n_big == t):>5} "
                  f"{str(inv):>12} {('{1,1}' if d % 2 == 0 else '{1,-1}'):>9}")
ok(malA4 == 0, "d*mu_q tiene t elementos y es inverso-cerrado en todo el rango")

print()
print("A5 -- Y AL NIVEL DE LOS h, quitar el plano es una DIFERENCIA FINITA explicita:")
print("        d par   ->  h^A = (1 - S)^2 h'        d impar ->  h^A = (1 - S^2) h'")
print("     con h'_n = C(n/q + d - 1, d - 1) los del alfabeto grande y S el desplazamiento.")
Nh = 24
malA5 = 0
for (m, k) in [(3, 2), (5, 2), (5, 3), (5, 4), (7, 4), (9, 5), (11, 6), (6, 3)]:
    from math import comb
    t = 2 * m + 2
    d = gcd(k, t)
    q = t // d
    hs = h_serie(m, k)
    hp = [comb(n // q + d - 1, d - 1) if n % q == 0 else 0 for n in range(N + 1)]
    if d % 2 == 0:
        pred = [hp[n] - 2 * (hp[n - 1] if n >= 1 else 0) + (hp[n - 2] if n >= 2 else 0)
                for n in range(N + 1)]
    else:
        pred = [hp[n] - (hp[n - 2] if n >= 2 else 0) for n in range(N + 1)]
    peor = max(abs(hs[n] - pred[n]) for n in range(min(Nh, N) + 1))
    malA5 += (peor > 1e-7)
    print(f"     m={m:>2} k={k:>2} d={d} q={q:>2}:  max|h_nuestro - (diferencia)(h_grande)| = {peor:.1e}")
ok(malA5 == 0, "la diferencia finita reproduce nuestros h a partir de los suyos")

print(SEP)
print("LECTURA.")
print("  El alfabeto del Teorema 2 de Prasad en t=1 es mu_n con multiplicidad m.  El de a_P^k es")
print("  mu_q con multiplicidad d, MENOS DOS PUNTOS.  Luego su m ES nuestro d y su n ES nuestro q,")
print("  y no por parecido: por identidad de multiconjuntos.  Y los dos puntos que se quitan son")
print("  +1 y +1 (d par) o +1 y -1 (d impar) --- que son exactamente los dos elementos fijados por")
print("  la involucion w -> w^{-1}, o sea LAS DOS CLASES FIJAS DEL PLEGADO.")
print()
print("  Eso convierte la tabla de la nota de parecido en correspondencia: las tres primeras filas")
print("  quedan demostradas, y la cuarta --- la naranja --- queda explicada, porque la diferencia")
print("  entre los dos enunciados es exactamente la diferencia entre los dos alfabetos.")
print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {sorted(set(fallos))[:6]}")
    sys.exit(1)
print("RESULTADO: la identidad de alfabetos se cumple, y su control con la paridad invertida falla.")
