# -*- coding: utf-8 -*-
u"""P3c_tres_cubicas.py -- las tres cubicas de las tablas son las tres cubicas clasicas, y el tamiz las relaciona.

Las tablas de P3c nombran a Clebsch, Cayley y Segre.  Lo que se mide aqui es que las tres no
son tres cubicas cualesquiera con nombre:

  * el estrato (1^5) es la superficie de CLEBSCH, de la que se hizo un modelo con sus
    27 rectas reales (Gotinga, 1872);
  * el estrato (2,1^4) es la cubica de CAYLEY, la unica superficie cubica con CUATRO nodos ---
    el maximo posible para una superficie cubica;
  * el estrato (1^6) es la cubica de SEGRE, la unica hipersuperficie cubica de P^4 con DIEZ
    nodos --- el maximo posible para un threefold cubico.

Es decir: los dos estratos mas largos de n = 6, (1^6) y (2,1^4), realizan los dos objetos
extremales de la geometria cubica del XIX, uno en cada dimension; el mas largo de n = 5 es
la Clebsch.

Y HAY UN TERCER HECHO, QUE ES EL BUENO.  Clasicamente, la cubica de Segre cortada por el
hiperplano y_i = y_j ES la cubica nodal de Cayley.  En el tamiz de Li-Wan, pasar de (1^6) a
(2,1^4) es FUNDIR DOS PARTES de la particion --- y fundir dos partes es exactamente poner
y_i = y_j.  La operacion combinatoria del tamiz y la seccion hiperplana clasica son la misma.
Eso no es una observacion sobre nombres: es una identidad de variedades, y aqui se comprueba.

LO QUE SE MIDE:
  1. V_(1^6) = {sum y = sum y^3 = 0} en P^5, con 10 puntos singulares.  (Segre)
  2. V_(2,1^4) = V_(1^6) cortada por y_1 = y_2, como IDENTIDAD de ecuaciones; y 4 nodos. (Cayley)
  3. V_(1^5) = {sum x = sum x^3 = 0} en P^4, lisa.  (Clebsch)
  4. Y el 5 de chi_p(5) del recuento de Clebsch es el cuerpo de definicion de sus 27 rectas:
     se cuentan las rectas sobre F_p y salen 27 cuando chi_p(5) = +1 y 15 cuando vale -1.
     Ese es el origen aritmetico del unico caracter que aparece en esa fila.

Todo por enumeracion exacta sobre F_p, sin algebra simbolica: la fuerza bruta aqui es barata y
no hay que creerle nada a nadie.  Controles: los nodos se identifican uno a uno (no solo se
cuentan), y los primos de las rectas cubren las dos ramas chi_p(5) = +1 y -1.

Uso:  python P3c_tres_cubicas.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import sys

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import io
import os
from itertools import product

ok = 0
mal = 0
L = []


def comprueba(que, cond, detalle=u''):
    global ok, mal
    if cond:
        ok += 1
        L.append(u'  ok   %-58s %s' % (que, detalle))
    else:
        mal += 1
        L.append(u'  MAL  %-58s %s' % (que, detalle))


def chi(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def proyectivos(l, p):
    u"""Un representante por punto de P^(l-1)(F_p): primera coordenada no nula igual a 1."""
    for i in range(l):
        for cola in product(range(p), repeat=l - 1 - i):
            yield (0,) * i + (1,) + cola


def estrato(w, p):
    u"""Los puntos de V_w(F_p) = {sum w_i x_i = sum w_i x_i^3 = 0} en P^(l-1)."""
    l = len(w)
    pts = []
    for x in proyectivos(l, p):
        if sum(w[i] * x[i] for i in range(l)) % p:
            continue
        if sum(w[i] * pow(x[i], 3, p) for i in range(l)) % p:
            continue
        pts.append(x)
    return pts


def singulares(w, p):
    u"""Los puntos donde los dos gradientes son proporcionales: la matriz jacobiana 2 x l
    tiene rango <= 1.  Los gradientes son (w_i) y (3 w_i x_i^2)."""
    l = len(w)
    out = []
    for x in estrato(w, p):
        f1 = [w[i] % p for i in range(l)]
        f2 = [(3 * w[i] * pow(x[i], 2, p)) % p for i in range(l)]
        rango1 = all((f1[i] * f2[j] - f1[j] * f2[i]) % p == 0
                     for i in range(l) for j in range(i + 1, l))
        if rango1:
            out.append(x)
    return out


# --------------------------------------------------------------------------
L.append(u'1.  SEGRE: el estrato (1^6) y sus DIEZ nodos')
L.append(u'-' * 78)
# --------------------------------------------------------------------------
W6 = (1, 1, 1, 1, 1, 1)
for p in (11, 13, 17, 19):
    s = singulares(W6, p)
    comprueba(u'p=%2d: V_(1^6) tiene 10 puntos singulares' % p, len(s) == 10,
          u'salen %d' % len(s))
# Y QUE SON, no solo cuantos: los nodos tienen que ser los repartos de signos tres y tres ---
# que es el lema del lugar singular, y la razon de que 10 = C(6,3)/2.  Volver a contar 10 no
# bastaria: seria un control que no puede fallar.
for p in (13, 19, 23):
    s = singulares(W6, p)
    bien = []
    for x in s:
        # normalizar el primer signo a +1 y mirar las seis coordenadas
        inv = pow(x[0], p - 2, p) if x[0] else None
        if inv is None:
            continue
        y = [(c * inv) % p for c in x]
        signos = [1 if c == 1 else (-1 if c == p - 1 else 0) for c in y]
        bien.append(0 not in signos and sum(signos) == 0)
    comprueba(u'p=%2d: los 10 nodos son repartos +-1 de suma cero (tres y tres)' % p,
          len(s) == 10 and len(bien) == 10 and all(bien),
          u'%d nodos, %d de la forma buena' % (len(s), sum(bien)))
comprueba(u'y por eso son 10: los repartos 3+3 de seis pesos, C(6,3)/2',
      (6 * 5 * 4) // 6 // 2 == 10, u'C(6,3)/2 = %d' % ((6 * 5 * 4) // 6 // 2))

# --------------------------------------------------------------------------
L.append(u'')
L.append(u'2.  CAYLEY: el estrato (2,1^4) es la SECCION y_1 = y_2 de Segre')
L.append(u'-' * 78)
# --------------------------------------------------------------------------
W5 = (2, 1, 1, 1, 1)
for p in (11, 13, 17, 19):
    s = singulares(W5, p)
    comprueba(u'p=%2d: V_(2,1^4) tiene 4 puntos singulares' % p, len(s) == 4,
          u'salen %d' % len(s))

# la identidad de ecuaciones, que es lo que de verdad se afirma
L.append(u'')
L.append(u'   la identidad, no la coincidencia numerica:')
L.append(u'   poner y_1 = y_2 = x_1 en  sum_{i=1}^{6} y_i = sum y_i^3 = 0  da')
L.append(u'        2x_1 + x_2 + x_3 + x_4 + x_5 = 0,   2x_1^3 + x_2^3 + ... + x_5^3 = 0,')
L.append(u'   que es exactamente V_(2,1^4).  FUNDIR DOS PARTES = CORTAR POR y_1 = y_2.')
for p in (11, 13, 17, 19, 23):
    corte = [x for x in estrato(W6, p) if x[0] == x[1]]
    # los puntos del corte, reescritos en las coordenadas del estrato fundido
    reescritos = set()
    for x in corte:
        y = (x[0], x[2], x[3], x[4], x[5])
        # normalizar a representante proyectivo
        for i in range(5):
            if y[i]:
                inv = pow(y[i], p - 2, p)
                y = tuple((c * inv) % p for c in y)
                break
        reescritos.add(y)
    directo = set(estrato(W5, p))
    comprueba(u'p=%2d: el corte y_1=y_2 de Segre ES V_(2,1^4), punto por punto' % p,
          reescritos == directo,
          u'%d vs %d' % (len(reescritos), len(directo)))

# --------------------------------------------------------------------------
L.append(u'')
L.append(u'3.  CLEBSCH: el estrato (1^5) es LISO')
L.append(u'-' * 78)
# --------------------------------------------------------------------------
W5c = (1, 1, 1, 1, 1)
for p in (11, 13, 17, 19, 23, 29):
    s = singulares(W5c, p)
    comprueba(u'p=%2d: V_(1^5) no tiene puntos singulares' % p, not s,
          u'salen %d' % len(s))

# --------------------------------------------------------------------------
L.append(u'')
L.append(u'4.  EL 5 DE chi_p(5) ES EL CUERPO DE DEFINICION DE LAS 27 RECTAS')
L.append(u'-' * 78)
L.append(u'    15 de las 27 rectas de la Clebsch son racionales y 12 se conjugan sobre Q(raiz 5);')
L.append(u'    asi que sobre F_p se ven 27 si 5 es resto cuadratico y 15 si no.  Ese es el unico')
L.append(u'    caracter de la fila, y no es una casualidad de ajuste.')
# --------------------------------------------------------------------------


def rectas(p):
    u"""Cuenta las rectas contenidas en V_(1^5) sobre F_p, por fuerza bruta sobre pares de
    puntos: una recta que une dos puntos de la superficie esta contenida en ella si y solo si
    todos sus puntos lo estan.  Se cuentan como conjuntos de puntos, no como pares."""
    pts = estrato(W5c, p)
    S = set(pts)
    vistas = set()
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            a, b = pts[i], pts[j]
            linea = set()
            dentro = True
            for t in range(p):
                q = tuple((a[k] + t * b[k]) % p for k in range(5))
                for k in range(5):
                    if q[k]:
                        inv = pow(q[k], p - 2, p)
                        q = tuple((c * inv) % p for c in q)
                        break
                if q not in S:
                    dentro = False
                    break
                linea.add(q)
            if not dentro:
                continue
            for k in range(5):
                if b[k]:
                    inv = pow(b[k], p - 2, p)
                    bb = tuple((c * inv) % p for c in b)
                    break
            if bb not in S:
                continue
            linea.add(bb)
            vistas.add(frozenset(linea))
    return len(vistas)


# LOS PRIMOS SE ELIGEN POR LAS DOS RAMAS, NO POR COMODIDAD.  5 es resto cuadratico
# exactamente cuando p = +-1 (mod 5); unos primos todos de esa clase no tocarian nunca la rama
# de las 15 rectas, que es justo la mitad que hace falta para que el caracter signifique algo.
PRIMOS_RECTAS = (11, 19, 29, 31,     # p = +-1 mod 5:  chi = +1, 27 rectas
                 7, 13, 17, 23)      # p = +-2 mod 5:  chi = -1, 15 rectas
for p in PRIMOS_RECTAS:
    n = rectas(p)
    esperado = 27 if chi(5, p) == 1 else 15
    comprueba(u'p=%2d (chi_p(5)=%+d): %d rectas sobre F_p' % (p, chi(5, p), n),
          n == esperado, u'esperadas %d' % esperado)
comprueba(u'las dos ramas del caracter se han probado, no solo una',
      len({chi(5, p) for p in PRIMOS_RECTAS}) == 2,
      u'chi vistos: %s' % sorted({chi(5, p) for p in PRIMOS_RECTAS}))

# --------------------------------------------------------------------------
L.append(u'')
L.append(u'VEREDICTO.  Los dos estratos mas largos de n = 6 realizan los dos')
L.append(u'objetos extremales de la geometria cubica del siglo XIX --- Cayley con el maximo de')
L.append(u'nodos para una superficie, Segre con el maximo para un threefold ---, el mas largo de')
L.append(u'n = 5 es la Clebsch, y la operacion del tamiz que lleva de Segre a Cayley es la seccion hiperplana')
L.append(u'clasica que los relaciona.  El unico caracter de la fila de Clebsch es el cuerpo de')
L.append(u'definicion de sus 27 rectas.')
L.append(u'')
L.append(u'=' * 70)
L.append(u'VERIFICACION LAS TRES CUBICAS:  %d ok, %d MAL' % (ok, mal))
L.append(u'=' * 70)

print(u'\n'.join(L))
io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     u'P3c_tres_cubicas_OUT.txt'), 'w',
        encoding='utf-8', newline='').write(u'\n'.join(L) + u'\n')
sys.exit(1 if mal else 0)
