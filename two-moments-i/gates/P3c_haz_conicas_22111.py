# -*- coding: utf-8 -*-
u"""P3c_haz_conicas_22111.py -- el haz de conicas de V_(2,2,1,1,1), derivado y verificado.

La derivacion por haz de conicas de la Proposicion 9 de P3c (etiqueta prop:22111) da la formula
    #V_(2,2,1,1,1)(F_p) = p^2 + 1 + p(3 + chi_p(5) + sum_{P(r)=0} chi_p(-2r)).
Aqui se rederiva entera, con P explicito, y se mide.

LA DERIVACION.  V_w = {sum w_i x_i = 0} cap {sum w_i x_i^3 = 0} en P^4, w = (2,2,1,1,1).
Con  s = x1+x2,  u = x1-x2,  t = x3+x4,  v = x3-x4  y  x5 = -(2s+t)  de la lineal, la cubica
queda, en P^3 de coordenadas (s:u:t:v),

    F = 2s^3 + 6su^2 + t^3 + 3tv^2 - 4(2s+t)^3 = 0.

  1. La recta  L : s = t = 0  esta en la superficie (F se anula identicamente).
  2. El plano  t = r s  del haz por L corta a F en L mas una conica residual:
         F|_{t=rs} = s * ( 6u^2 + 3r v^2 + B(r) s^2 ),   B(r) = 2 + r^3 - 4(2+r)^3 = -3 P(r),
     con   P(r) = r^3 + 8r^2 + 16r + 10.
     Dividiendo por 3:   C_r :  2u^2 + r v^2 - P(r) s^2 = 0   en (s:u:v).
  3. disc C_r = -2 r P(r).  Luego las fibras degeneradas son r = 0, las raices de P, y r = oo.
  4. Recuento fibra a fibra (p impar, conica lisa ternaria = p+1 puntos siempre):
       r generico   : p+1 puntos, corta a L en 1 + chi(-2r)  ->  aporta p - chi(-2r)
       r = 0        : 2u^2 - 10 s^2, rango 2, par de rectas racional sii 20 = cuadrado
                      ->  p+1+p*chi(5) puntos, corta a L en 1  ->  aporta p + p*chi(5)
       P(r) = 0     : 2u^2 + r v^2, rango 2  ->  p+1+p*chi(-2r), corta a L en 1+chi(-2r)
                      ->  aporta p + (p-1)*chi(-2r)
       r = oo       : el plano s=0 da 3t(v^2-t^2) = 0, tres rectas  ->  aporta 2p
     y sum_{r in F_p} chi(-2r) = 0 reabsorbe los genericos.  Sumando con #L = p+1 sale
         #V = p^2 + 1 + p(3 + chi_p(5) + sum_{P(r)=0, r in F_p} chi_p(-2r)).

LO QUE SE MIDE AQUI:
  0. CONTROL POSITIVO: el mismo contador sobre w = (1,1,1,1,1) tiene que dar la Clebsch
     p^2 + (6+chi_p(5))p + 1, que es un numero PUBLICADO en la nota.  Sin esto el contador no
     acredita nada.
  1. F se anula en L, identicamente en (u:v)  --- simbolico.
  2. La factorizacion  F|_{t=rs} = s*(6u^2+3rv^2+B(r)s^2)  y  B = -3P  --- simbolico.
  3. La formula contra el recuento por fuerza bruta, primo a primo.
  4. CONTROL NEGATIVO: quitar el termino chi_p(5), o cambiar chi(-2r) por chi(2r), tiene que FALLAR.
  5. Los primos malos se nombran y se explican por disc(P) y los coeficientes, no se esconden.

Uso:  python P3c_haz_conicas_22111.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys

import numpy as np
import sympy as sp

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

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


# ---------------------------------------------------------------- 1-2. lo simbolico
s, u, t, v, r = sp.symbols('s u t v r')

F = sp.expand(2 * s ** 3 + 6 * s * u ** 2 + t ** 3 + 3 * t * v ** 2 - 4 * (2 * s + t) ** 3)

# control del control: F tiene que ser la cubica de partida, reconstruida desde las x
x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
x5 = -(2 * x1 + 2 * x2 + x3 + x4)
cubica = sp.expand(2 * x1 ** 3 + 2 * x2 ** 3 + x3 ** 3 + x4 ** 3 + x5 ** 3)
sust = {x1: (s + u) / 2, x2: (s - u) / 2, x3: (t + v) / 2, x4: (t - v) / 2}
comprueba(u'control: F es 4x la cubica original en las coordenadas (s,u,t,v)',
      sp.simplify(sp.expand(cubica.subs(sust) * 4) - F) == 0)

comprueba(u'1. la recta L : s = t = 0 esta en la superficie',
      sp.simplify(F.subs({s: 0, t: 0})) == 0)

Frs = sp.expand(F.subs(t, r * s))
coc = sp.simplify(sp.cancel(Frs / s))
B = sp.expand(sp.simplify(coc - 6 * u ** 2 - 3 * r * v ** 2)) / s ** 2
B = sp.expand(sp.simplify(B))
P = sp.Poly(r ** 3 + 8 * r ** 2 + 16 * r + 10, r)

comprueba(u'2. F|_{t=rs} = s * (6u^2 + 3r v^2 + B(r) s^2)',
      sp.simplify(Frs - s * (6 * u ** 2 + 3 * r * v ** 2 + B * s ** 2)) == 0)
comprueba(u'   ... y B(r) = -3 P(r) con P = r^3 + 8r^2 + 16r + 10',
      sp.expand(B + 3 * P.as_expr()) == 0, u'B = %s' % sp.factor(B))
comprueba(u'   ... disc de la conica = -2 r P(r), luego degenera solo en r=0 y P(r)=0',
      sp.expand(sp.det(sp.diag(-P.as_expr(), 2, r)) + 2 * r * P.as_expr()) == 0)
comprueba(u'   ... P(0) = 10 != 0, luego r = 0 NO es raiz de P (las ramas no se solapan)',
      P.eval(0) == 10)

discP = sp.discriminant(P.as_expr(), r)
comprueba(u'   ... P es irreducible sobre Q y su grupo es S_3 (disc no cuadrado)',
      sp.Poly(P, r).is_irreducible and not sp.sqrt(sp.Integer(discP)).is_Integer,
      u'disc P = %s' % discP)

# ---------------------------------------------------------------- 3. la medida
def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def contar(w, p):
    u"""#V_w(F_p) para l = 5, eliminando x5 con la lineal.

    Se itera sobre x1 y se vectoriza el resto: el array vivo es p^3 y no p^4, que es lo que
    permite llegar a p ~ 150 sin agotar memoria.
    """
    w1, w2, w3, w4, w5 = [x % p for x in w]
    inv5 = pow(w5, p - 2, p)
    g = np.arange(p, dtype=np.int64)
    X2 = g.reshape(p, 1, 1)
    X3 = g.reshape(1, p, 1)
    X4 = g.reshape(1, 1, p)
    cub = (g ** 3) % p                      # tabla de cubos, una vez
    base = (w2 * cub[X2] + w3 * cub[X3] + w4 * cub[X4]) % p
    lin = (w2 * X2 + w3 * X3 + w4 * X4) % p
    ceros = 0
    for a in range(p):
        X5 = (-(w1 * a + lin) * inv5) % p
        Fv = (w1 * cub[a] + base + w5 * cub[X5]) % p
        ceros += int(np.count_nonzero(Fv == 0))
    return (ceros - 1) // (p - 1)


def raices_P(p):
    return [rr for rr in range(p) if (rr ** 3 + 8 * rr ** 2 + 16 * rr + 10) % p == 0]


def formula(p, signo=-2, con_chi5=True):
    S = sum(chi(signo * rr, p) for rr in raices_P(p))
    return p * p + 1 + p * (3 + (chi(5, p) if con_chi5 else 0) + S)


PRIMOS = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
          79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149]
# p = 3 y p = 5 van aparte: la derivacion divide por 3, y en p = 5 la rama r = 0 y la rama
# P(r) = 0 se SOLAPAN porque P(0) = 10.  Hay que mirarlos, no darlos por excluidos.
PRIMOS_PEQ = [3, 5]

# control positivo sobre un numero PUBLICADO
malos_clebsch = [p for p in PRIMOS
                 if contar((1, 1, 1, 1, 1), p) != p * p + (6 + chi(5, p)) * p + 1]
comprueba(u'0. CONTROL POSITIVO: el contador reproduce la Clebsch publicada, p^2+(6+chi_p(5))p+1',
      not malos_clebsch, u'falla en: %s' % malos_clebsch)

medido = {}
for p in PRIMOS:
    medido[p] = contar((2, 2, 1, 1, 1), p)

malos = [(p, medido[p], formula(p)) for p in PRIMOS if medido[p] != formula(p)]
buenos = [p for p in PRIMOS if medido[p] == formula(p)]

comprueba(u'3. la formula acierta en TODOS los primos probados, sin excepciones',
      not malos, u'discrepan: %s' % ([p for (p, _, _) in malos],))
comprueba(u'   ... y son %d primos, 7 <= p <= %d' % (len(buenos), max(buenos)),
      len(buenos) >= 18)

# 5. el primo de mala reduccion NO rompe la formula: disc(P) = -140 se anula en p = 7, P tiene
# ahi raiz doble, y la formula AGUANTA igual.  Por eso p = 7 se mide y no se excluye.
comprueba(u'5. 7 divide a disc(P) = %s y P tiene raiz DOBLE mod 7' % discP,
      discP % 7 == 0 and len(raices_P(7)) == 2,
      u'raices de P mod 7: %s' % raices_P(7))
comprueba(u'   ... y AUN ASI la formula acierta en p = 7 (la suma sobre las 2 raices da -2)',
      medido[7] == formula(7),
      u'medido %d, formula %d' % (medido[7], formula(7)))
comprueba(u'   ... la rama de TRES raices esta ejercitada, no es teorica',
      any(len(raices_P(p)) == 3 for p in buenos),
      u'p con P totalmente escindido: %s'
      % [p for p in buenos if len(raices_P(p)) == 3])

# p = 3 y p = 5: se miden, no se declaran excluidos de oidas.
peq = [(p, contar((2, 2, 1, 1, 1), p), formula(p)) for p in PRIMOS_PEQ]
# En p = 5, P(0) = 10 = 0 y las dos ramas degeneradas de la derivacion se SOLAPAN (r = 0 es
# raiz de P).  No falla: el solape colapsa dos aportaciones en una y la cuenta sale igual.
comprueba(u'6. p = 5 AGUANTA, aunque ahi r = 0 ES raiz de P y las dos ramas se solapan',
      [(m == f) for (p, m, f) in peq if p == 5] == [True],
      u'p=5: medido %d, formula %d' % (peq[1][1], peq[1][2]))
comprueba(u'   ... p = 3 es el UNICO que falla, y es donde la derivacion divide por 3',
      [(m == f) for (p, m, f) in peq if p == 3] == [False],
      u'p=3: medido %d, formula %d' % (peq[0][1], peq[0][2]))
comprueba(u'   ... luego la hipotesis correcta es p > 3, no p > 5 ni "p de buena reduccion"',
      not malos and all(m == f for (p, m, f) in peq if p > 3))

# ---------------------------------------------------------------- 4. controles negativos
d1 = [p for p in buenos if medido[p] != formula(p, con_chi5=False)]
comprueba(u'4. CONTROL NEGATIVO: quitar el termino chi_p(5) FALLA', len(d1) > 0,
      u'falla en %d de %d primos buenos' % (len(d1), len(buenos)))
d2 = [p for p in buenos if medido[p] != formula(p, signo=2)]
comprueba(u'   CONTROL NEGATIVO: chi(+2r) en vez de chi(-2r) FALLA', len(d2) > 0,
      u'falla en %d de %d primos buenos' % (len(d2), len(buenos)))
# y una variante que NO puede distinguirse: -2r y -8r son el mismo caracter (difieren en un cuadrado)
d3 = [p for p in buenos if medido[p] != formula(p, signo=-8)]
comprueba(u'   CONTROL DE EQUIVALENCIA: chi(-8r) NO puede fallar --- -8r y -2r difieren en el cuadrado 4',
      len(d3) == 0)

L.append(u'')
L.append(u'  p     #V medido   formula   raices de P mod p    sum chi(-2r)')
L.append(u'  ' + u'-' * 68)
for p in PRIMOS:
    rr = raices_P(p)
    L.append(u'  %-5d %-11d %-9d %-20s %d'
             % (p, medido[p], formula(p), rr if rr else u'ninguna',
                sum(chi(-2 * z, p) for z in rr)))

L.append(u'')
L.append(u'VEREDICTO')
L.append(u'  #V_(2,2,1,1,1)(F_p) = p^2 + 1 + p( 3 + chi_p(5) + sum_{P(r)=0} chi_p(-2r) )')
L.append(u'  con P(r) = r^3 + 8r^2 + 16r + 10, irreducible sobre Q, grupo S_3, disc = %s.' % discP)
L.append(u'  Derivada por el haz de conicas por la recta s = t = 0, no ajustada: los cuatro')
L.append(u'  pasos simbolicos estan arriba y la medida es un control, no la fuente.')
L.append(u'  Hipotesis: p > 3, y NADA MAS.  Los dos primos sospechosos no hay que excluirlos:')
L.append(u'  p = 7, donde P tiene raiz doble, y p = 5, donde r = 0 es raiz de P y las dos')
L.append(u'  ramas degeneradas se solapan.  Los dos aguantan; el unico que cae es p = 3.')

print(u'\n'.join(L))
print(u'')
print(u'=' * 72)
print(u'VERIFICACION HAZ DE CONICAS DE V_(2,2,1,1,1):  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 72)
sys.exit(1 if MAL else 0)
