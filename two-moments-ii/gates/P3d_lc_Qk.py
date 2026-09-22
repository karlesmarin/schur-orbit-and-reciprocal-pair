# -*- coding: utf-8 -*-
u"""P3d_lc_Qk.py -- lc(Q_k) | k para todo k >= 3, que es lo que hace valer prop:s3 para todo k.

La divisibilidad A = lc(Q_k) | k se demuestra con el coeficiente de Waring del monomio
extremo, sin computo; la comprobacion numerica de abajo la contrasta por otra ruta (la
recursion de Newton) en 3 <= k <= 400, y lleva un control negativo que tiene que fallar.

LA DEMOSTRACION, entera.

  Con e_1 = 0, la formula de Waring da para  s_k = sum_{2a+3b=k} c_{a,b} e2^a e3^b

        |c_{a,b}|  =  k (a+b-1)! / (a! b!).

  Los monomios de s_k son una progresion:  (a, b) = (a_0 + 3i, b_max - 2i),  i = 0..d,
  porque 2a + 3b = k.  Como  u = e2^3 / e3^2,  el grado en u ES i, luego el coeficiente
  principal crudo es el del monomio con b MINIMO, y  b_min = k mod 2  (3b = k - 2a fuerza
  b = k mod 2).  Ahi el factorial de arriba se cancela contra el de abajo:

    k par:    b = 0,  a = k/2,      a+b = k/2
              |c| = k (k/2 - 1)! / ((k/2)! 0!)      = k / (k/2)  =  2

    k impar:  b = 1,  a = (k-3)/2,  a+b = (k-1)/2
              |c| = k ((k-3)/2)! / (((k-3)/2)! 1!)  =  k

  Q_k es ese polinomio PRIMITIVADO, asi que  lc(Q_k) = |crudo| / cont(Q_k)  divide al
  crudo.  Luego lc(Q_k) divide a 2 si k es par --- y 2 | k --- y divide a k si k es impar.
  En los dos casos  lc(Q_k) | k.   []

  No se usa nada de Q_k salvo su monomio extremo, asi que vale para todo k >= 3 sin computo.

Y DE PROPINA, LA FORMA CERRADA --- pero esta es MEDIDA, no demostrada, y se dice asi.
El contenido resulta valer

        cont(Q_k) = 2  si k es potencia de 2
                  = p  si k = p^a con p primo impar
                  = 1  en el resto

  y por tanto   lc(Q_k) = 1,  2,  k/p,  k   en las cuatro ramas respectivamente.
  Medido 3 <= k <= 400.  Lo que le falta para ser teorema: el contenido es el gcd de
  (k/n) C(n, a) con n = n_0 + i variando CON i, asi que NO es el teorema clasico del gcd
  de los binomiales de un n fijo.  Se dice asi y no se escribe como demostrado.

LO QUE MIDE ESTA COMPUERTA (y cada cosa por una ruta distinta de la que quiere confirmar):
  1. b_min = k mod 2.
  2. Waring, calculada aparte, reproduce el coeficiente que da la RECURSION en ese monomio.
  3. El crudo vale 2 (k par) / k (k impar).
  4. CONTROL NEGATIVO: "el crudo vale k tambien para k par" tiene que FALLAR, y fallar exactamente
     en los k pares.  Sin el, el punto 3 no distingue las dos ramas.
  5. cont(Q_k) por ramas, con las cuatro ramas no vacias.
  6. El corolario lc(Q_k) | k, medido aparte de la formula.
  7. Control del control: la recursion reproduce s_3 = 3 e_3 y s_4 = 2 e_2^2.

Uso:  python P3d_lc_Qk.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys
from collections import Counter
from math import factorial

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

KMAX = 400

try:
    from math import gcd
except ImportError:                                    # py2
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a


# ---------------------------------------------------------------- s_k, sin sympy
# s_k es isobarico de peso k con peso(e2) = 2, peso(e3) = 3, asi que un dict
# {(a, b): coef} es representacion exacta y la recursion es aritmetica entera.
S = {0: {(0, 0): 3}, 1: {}, 2: {(1, 0): -2}, 3: {(0, 1): 3}}


def mul_shift(d, da, db, mult):
    return dict(((a + da, b + db), c * mult) for (a, b), c in d.items())


def add(d1, d2):
    r = dict(d1)
    for key, v in d2.items():
        r[key] = r.get(key, 0) + v
        if r[key] == 0:
            del r[key]
    return r


for k in range(4, KMAX + 1):
    S[k] = add(mul_shift(S[k - 2], 1, 0, -1), mul_shift(S[k - 3], 0, 1, 1))


def waring(k, a, b):
    u"""|coef de e2^a e3^b en s_k| segun Waring:  k (a+b-1)! / (a! b!)."""
    return k * factorial(a + b - 1) // (factorial(a) * factorial(b))


def factoriza(n):
    f = Counter()
    m, p = n, 2
    while p * p <= m:
        while m % p == 0:
            f[p] += 1
            m //= p
        p += 1
    if m > 1:
        f[m] += 1
    return f


def rama_de(k):
    primos = sorted(factoriza(k))
    if k % 2 == 0:
        return ('potencia de 2', 2) if primos == [2] else ('par no 2-adico', 1)
    if len(primos) == 1:
        return ('p^a impar', primos[0])
    return ('impar compuesto', 1)


# ---------------------------------------------------------------- comprobaciones
OK = MAL = 0
L = []


def comprueba(et, cond, det=u''):
    global OK, MAL
    if cond:
        OK += 1
        L.append(u'  ok   %s%s' % (et, (u'   ' + det) if det else u''))
    else:
        MAL += 1
        L.append(u'  MAL  %s%s' % (et, (u'   ' + det) if det else u''))


# 7. control del control, primero: si la recursion esta mal, todo lo de abajo miente
comprueba(u'control: la recursion da s_3 = 3 e_3', S[3] == {(0, 1): 3})
comprueba(u'control: la recursion da s_4 = 2 e_2^2', S[4] == {(2, 0): 2})

mal_bmin, mal_waring, mal_crudo, mal_cont, mal_div = [], [], [], [], []
control_neg_pares = []
ramas = Counter()
tabla = []

for k in range(3, KMAX + 1):
    d = S[k]
    b_min = min(b for (a, b) in d)
    a_top = max(a for (a, b) in d if b == b_min)

    if b_min != k % 2:
        mal_bmin.append((k, b_min))

    crudo = abs(d[(a_top, b_min)])
    if crudo != waring(k, a_top, b_min):
        mal_waring.append((k, crudo, waring(k, a_top, b_min)))

    esperado = 2 if k % 2 == 0 else k
    if crudo != esperado:
        mal_crudo.append((k, crudo, esperado))
    if crudo != k:
        control_neg_pares.append(k)

    # Q_k primitivo: el indice en u es i = (a - a_min)/3
    a_min = min(a for (a, b) in d)
    co = {}
    for (a, b), c in d.items():
        assert (a - a_min) % 3 == 0, (k, a, b)
        co[(a - a_min) // 3] = c
    grado = max(co)
    coef = [co.get(i, 0) for i in range(grado + 1)]
    g = 0
    for c in coef:
        g = gcd(g, abs(c))
    A = abs(coef[-1]) // g

    if A == 0 or k % A != 0:
        mal_div.append((k, A))

    nombre, cont_esperado = rama_de(k)
    ramas[nombre] += 1
    if g != cont_esperado:
        mal_cont.append((k, g, cont_esperado, nombre))
    if k <= 30:
        tabla.append((k, grado, crudo, g, A, nombre))

comprueba(u'1. b_min = k mod 2 (el monomio extremo es b=0 si k par, b=1 si k impar)',
      not mal_bmin, u'fallos: %s' % (mal_bmin[:5],))
comprueba(u'2. Waring reproduce el coeficiente de la RECURSION en el monomio extremo',
      not mal_waring, u'fallos: %s' % (mal_waring[:5],))
comprueba(u'3. el coeficiente principal CRUDO vale 2 (k par) / k (k impar)',
      not mal_crudo, u'fallos: %s' % (mal_crudo[:5],))
comprueba(u'4. CONTROL NEGATIVO: "el crudo vale k tambien para k par" FALLA',
      len(control_neg_pares) > 0, u'%d casos' % len(control_neg_pares))
comprueba(u'   ... y falla EXACTAMENTE en los k pares (si no, el punto 3 no separa ramas)',
      sorted(control_neg_pares) == [k for k in range(3, KMAX + 1) if k % 2 == 0])
comprueba(u'5. cont(Q_k) = 2 / p / 1 por ramas  (MEDIDO, no demostrado)',
      not mal_cont, u'fallos: %s' % (mal_cont[:5],))
for r in (u'potencia de 2', u'par no 2-adico', u'p^a impar', u'impar compuesto'):
    comprueba(u'   ... rama "%s" no vacia' % r, ramas[r] > 0, u'%d casos' % ramas[r])
comprueba(u'6. COROLARIO lc(Q_k) | k, medido aparte de la formula, 3 <= k <= %d' % KMAX,
      not mal_div, u'fallos: %s' % (mal_div[:5],))


# ---------------------------------------------------------------- el pago: prop:s3
# Con A | k demostrada, A <= k, luego  k < 2^d  ==>  A < 2^d.  La pregunta deja de
# necesitar el valor de A y pasa a ser aritmetica pura sobre d = floor(k/6) - [k=1 mod 6].
def d_de(k):
    return k // 6 - (1 if k % 6 == 1 else 0)


excepciones = [k for k in range(3, 5000) if d_de(k) > 0 and k >= 2 ** d_de(k)]

comprueba(u'7. con A <= k, el conjunto excepcional {k >= 2^d} es finito y su maximo es 37',
      excepciones and max(excepciones) == 37, u'%d valores, max %s'
      % (len(excepciones), max(excepciones) if excepciones else None))
comprueba(u'   ... y esta CONTENIDO en 6 <= k <= 46, que es donde hay certificados exactos',
      all(6 <= k <= 46 for k in excepciones),
      u'fuera del rango: %s' % ([k for k in excepciones if not 6 <= k <= 46],))
comprueba(u'   ... k < 2^d para todo 38 <= k <= 4999 (la cola no reabre)',
      all(k < 2 ** d_de(k) for k in range(38, 5000) if d_de(k) > 0))
comprueba(u'   ... y el margen 2^d/k en la cola no es apretado (>= 1.4 ya en 38 <= k < 200)',
      min((2 ** d_de(k)) / float(k) for k in range(38, 200) if d_de(k) > 0) >= 1.4,
      u'minimo %.3f' % min((2 ** d_de(k)) / float(k)
                           for k in range(38, 200) if d_de(k) > 0))
# Control de que la cota demostrada es ESTRICTAMENTE mas floja que la medida, y aun asi
# basta.  Si saliera igual, el trabajo de arriba no habria comprado nada.
def A_medida(k):
    d = S[k]
    a_min = min(a for (a, b) in d)
    co = {}
    for (a, b), c in d.items():
        co[(a - a_min) // 3] = c
    coef = [co.get(i, 0) for i in range(max(co) + 1)]
    g = 0
    for c in coef:
        g = gcd(g, abs(c))
    return abs(coef[-1]) // g


exc_medida = [k for k in range(3, KMAX + 1)
              if d_de(k) > 0 and A_medida(k) >= 2 ** d_de(k)]

comprueba(u'8. con A MEDIDA el conjunto excepcional es {6,9,10,15,21,33,35} (el de prop:s3)',
      exc_medida == [6, 9, 10, 15, 21, 33, 35], u'%s' % (exc_medida,))
comprueba(u'   ... la cota demostrada A <= k es ESTRICTAMENTE mas floja (29 casos frente a 7)',
      set(exc_medida) < set(excepciones),
      u'medida %d, demostrada %d' % (len(exc_medida), len(excepciones)))
comprueba(u'   ... y aun asi basta, porque los 29 caben en el rango certificado',
      max(excepciones) <= 46)

# ap:patron: cuando es monico Q_k (con deg Q_k > 0).  Helou da una direccion (k primo => monico);
# la otra es MEDIDA, no demostrada, y para k par ni siquiera es cierta.
es_primo = lambda n: n > 1 and all(n % q for q in range(2, int(n ** 0.5) + 1))
mon_impar = [k for k in range(3, KMAX + 1, 2) if d_de(k) > 0 and A_medida(k) == 1]
pri_impar = [k for k in range(3, KMAX + 1, 2) if d_de(k) > 0 and es_primo(k)]
comprueba(u'9. k impar, 3 <= k <= %d: Q_k monico EXACTAMENTE para k primo (MEDIDO)' % KMAX,
          mon_impar == pri_impar, u'%d casos' % len(mon_impar))
mon_par = [k for k in range(4, KMAX + 1, 2) if d_de(k) > 0 and A_medida(k) == 1]
comprueba(u'   ... y para k par los monicos son las potencias de 2: el "exactamente primo" FALLA',
          mon_par == [k for k in range(4, KMAX + 1, 2) if d_de(k) > 0 and k & (k - 1) == 0]
          and 8 in mon_par and 16 in mon_par, u'%s' % (mon_par,))

L.append(u'')
L.append(u'  k    deg  crudo   cont   lc(Q_k)   rama')
L.append(u'  ' + u'-' * 58)
for (k, grado, crudo, g, A, nombre) in tabla:
    L.append(u'  %-4d %-4d %-7d %-6d %-9d %s' % (k, grado, crudo, g, A, nombre))

L.append(u'')
L.append(u'VEREDICTO')
L.append(u'  lc(Q_k) | k queda DEMOSTRADA para todo k >= 3: el coeficiente principal crudo')
L.append(u'  es 2 (k par) o k (k impar) por Waring en el monomio b = k mod 2, y primitivar')
L.append(u'  solo puede dividirlo.  prop:s3 vale ahora para todo k sin el computo finito.')
L.append(u'  La FORMA CERRADA de lc(Q_k) --- 1, 2, k/p, k --- sigue siendo MEDIDA:')
L.append(u'  le falta el contenido, que es el gcd de (k/n)C(n,a) con n variando con i, y')
L.append(u'  eso NO es el teorema clasico del gcd de binomiales de un n fijo.')

print(u'\n'.join(L))
print(u'')
print(u'=' * 72)
print(u'VERIFICACION lc(Q_k) | k:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 72)
sys.exit(1 if MAL else 0)
