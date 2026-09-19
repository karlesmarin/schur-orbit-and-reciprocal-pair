# -*- coding: utf-8 -*-
r"""inter_autocorrelacion.py -- la prueba de lem:inter por autocorrelacion, medida paso a paso.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

LA PRUEBA NUEVA (19-sep-2026).  Sustituye a la de seis casos (cotas phi(q') <= 2 sqrt(q') - 2 e
inspeccion de q' en {4,6,10,12,18,30}).  Un intervalo que es union de estratos es estable por toda
unidad u.  Sea J el mas corto de I y su complemento, de longitud b con 2 <= b <= m/2, y
C_J(t) = |J cap (J+t)|.  Si uJ = J entonces C_J(ut) = C_J(t), y en particular C_J(u) = C_J(1) =
b-1.  Pero, para 1 <= t < m,

    C_J(t) = max(0, b-t) + max(0, b-(m-t)),

y con 2b <= m los dos sumandos no son positivos a la vez; luego C_J(t) = b-1 sii t = 1 o t = m-1.
Asi u = +-1.  Es MAS FUERTE que el enunciado: basta UNA unidad u != +-1 para que ningun intervalo
de longitud 2..m-2 sea u-estable.

QUE SE MIDE.
  P1  la formula de C_J(t), contra el recuento directo, m <= 120, todo b en [1, m-1], todo t.
  P2  la forma fuerte: para m <= 120, todo intervalo con 2 <= B <= m-2 y toda unidad u,
      uI = I  =>  u = +-1 (mod m).
  P3  la clasificacion del lema: intervalos union de estratos con 2 <= B <= m-2, m <= 150;
      debe salir exactamente (6,3) con {5,0,1} y {2,3,4}.
  N1  CASI-ACIERTO (la compuerta sabe fallar): en los extremos B = 1, m-1 SI hay intervalos
      estables por u != +-1 -- el singleton {0}.  Si P2 se midiera sin la cota 2 <= B, fallaria.
  N2  CASI-ACIERTO: un conjunto que NO es intervalo, {1,5} con m = 24 y u = 5, es u-estable.
      La hipotesis "intervalo" es la que trabaja.

Salida: gates/capas_tipo/inter_autocorrelacion_OUT.txt
"""
import sys
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FALLOS = 0


def intervalo(a, B, m):
    return frozenset((a + i) % m for i in range(B))


def C(J, t, m):
    return len(J & frozenset((x + t) % m for x in J))


def unidades(m):
    return [u for u in range(1, m) if gcd(u, m) == 1]


def estable(I, u, m):
    return frozenset((u * x) % m for x in I) == I


# P1
n1 = 0
for m in range(3, 121):
    for b in range(1, m):
        J = intervalo(0, b, m)
        for t in range(1, m):
            directo = C(J, t, m)
            formula = max(0, b - t) + max(0, b - (m - t))
            if 2 * b <= m and directo != formula:
                FALLOS += 1
                print("P1 FALLO m=%d b=%d t=%d directo=%d formula=%d" % (m, b, t, directo, formula))
            n1 += 1
print("P1  C_J(t) = max(0,b-t)+max(0,b-(m-t)) con 2b<=m, m<=120: %d medidas" % n1)

# P2
n2 = 0
for m in range(3, 121):
    U = unidades(m)
    for B in range(2, m - 1):
        for a in range(m):
            I = intervalo(a, B, m)
            for u in U:
                n2 += 1
                if estable(I, u, m) and u not in (1, m - 1):
                    FALLOS += 1
                    print("P2 FALLO m=%d B=%d a=%d u=%d" % (m, B, a, u))
print("P2  uI=I => u=+-1, 2<=B<=m-2, m<=120: %d medidas" % n2)

# P3
hallados = set()
for m in range(3, 151):
    estr = {}
    for c in range(m):
        estr.setdefault(gcd(c, m), set()).add(c)
    for B in range(2, m - 1):
        for a in range(m):
            I = intervalo(a, B, m)
            if all(s <= I or not (s & I) for s in map(frozenset, estr.values())):
                hallados.add((m, B, tuple(sorted(I))))
esperado = {(6, 3, (0, 1, 5)), (6, 3, (2, 3, 4))}
if hallados != esperado:
    FALLOS += 1
    print("P3 FALLO: %s" % sorted(hallados))
print("P3  intervalos union de estratos, 2<=B<=m-2, m<=150: %s" % sorted(hallados))

# N1
testigos = [(m, u) for m in range(5, 40) for u in unidades(m)
            if u not in (1, m - 1) and estable(intervalo(0, 1, m), u, m)]
if not testigos:
    FALLOS += 1
    print("N1 FALLO: el singleton {0} no aparece como estable")
print("N1  extremo B=1: {0} es u-estable con u!=+-1 en %d pares (m,u), p.ej. %s"
      % (len(testigos), testigos[:3]))

# N2
if not (estable(frozenset({1, 5}), 5, 24)):
    FALLOS += 1
    print("N2 FALLO")
print("N2  {1,5} (no intervalo), m=24, u=5: estable = %s" % estable(frozenset({1, 5}), 5, 24))

print("\nFALLOS = %d" % FALLOS)
sys.exit(1 if FALLOS else 0)
