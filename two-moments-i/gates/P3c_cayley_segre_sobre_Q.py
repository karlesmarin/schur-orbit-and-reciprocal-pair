# -*- coding: utf-8 -*-
u"""P3c_cayley_segre_sobre_Q.py -- los modelos de Segre y Cayley, sobre el cuerpo base y no solo sobre la clausura.

QUE SE ESTABLECE.  El apendice ap:extremales de P3c dice: "a cubic surface can have at most
four nodes, and there is exactly one that attains it: Cayley's nodal cubic ... lambda=(2,1^4)
gives Cayley, lambda=(1^6) gives Segre".  Las dos unicidades son GEOMETRICAS --- sobre k-barra y
salvo equivalencia proyectiva.  Sobre Q puede haber formas torcidas del mismo tipo geometrico, y
una torcida CUENTA DISTINTO sobre F_p.  Como el articulo entero va de contar sobre F_p, la
identificacion tiene que ser sobre el cuerpo base o no vale para lo que se usa.

LO QUE SE ESTABLECE AQUI:

  (a) SEGRE no necesita equivalencia ninguna.  El modelo clasico de la cubica de Segre es
      {sum y_i = 0} cap {sum y_i^3 = 0} en P^5, que es V_(1^6) LITERALMENTE, sobre Z.  Sus diez
      nodos son los vectores +-1 con tres +1, todos racionales.  No hay torcida posible: es el
      mismo esquema, no una forma de el.

  (b) CAYLEY si la necesita, y la hay sobre Q.  Los cuatro nodos de V_(2,1^4) son RACIONALES y
      estan en posicion general (det = -8), luego la matriz M cuyas COLUMNAS son los nodos esta
      en GL_4(Q) --- de hecho en GL_4(Z[1/2]) --- y

          F(M u)  =  24 * ( u0u1u2 + u0u1u3 + u0u2u3 + u1u2u3 ),

      que es la forma estandar de Cayley EXACTA, sin reescalar variables.  Luego V_(2,1^4) es
      proyectivamente equivalente a la cubica nodal de Cayley sobre el CUERPO PRIMO, siempre que
      la caracteristica no divida a 6.

  (c) Y la hipotesis char no divide 6 es AFILADA, no de conveniencia: det M = -8 se anula en
      p = 2 y la constante 24 se anula en p = 2 y en p = 3.  Se mide, no se supone.

  La racionalidad de los nodos es lo que DECIDE (b): cuatro nodos racionales en posicion general
  se mandan al simplex coordenado por una matriz del cuerpo base, y una cubica con nodos en los
  cuatro puntos coordenados esta forzada a la forma  a u0u1u2 + b u0u1u3 + c u0u2u3 + d u1u2u3.

CONTROLES QUE PUEDEN FALLAR, y por eso estan:
  - los nodos se CALCULAN del jacobiano, no se copian del apendice;
  - la identidad F(Mu) = 24*Cayley es simbolica, no numerica;
  - CONTROL NEGATIVO: sustituir un nodo por otro punto racional cualquiera tiene que ROMPER la forma;
  - y el control de la ruta de ida y vuelta: la seccion y_1 = y_2 de Segre ES V_(2,1^4).

Uso:  python P3c_cayley_segre_sobre_Q.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import sys
from itertools import combinations

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


# ---------------------------------------------------------------- (a) SEGRE
y = sp.symbols('y0:6')
Lin = sum(y)
Cub = sum(t ** 3 for t in y)

pm = set()
for S in combinations(range(6), 3):
    v = tuple(1 if i in S else -1 for i in range(6))
    pm.add(min(v, tuple(-c for c in v)))          # clases proyectivas v ~ -v

comprueba(u'a1. los vectores +-1 de suma cero, salvo signo, son DIEZ', len(pm) == 10)

sing = []
for v in sorted(pm):
    d = dict(zip(y, v))
    if sp.simplify(Lin.subs(d)) != 0 or sp.simplify(Cub.subs(d)) != 0:
        continue
    J = sp.Matrix([[sp.diff(Lin, t).subs(d) for t in y],
                   [sp.diff(Cub, t).subs(d) for t in y]])
    if J.rank() < 2:
        sing.append(v)

comprueba(u'a2. los diez estan en V_(1^6) y son SINGULARES de la interseccion completa',
      len(sing) == 10)
comprueba(u'a3. y son racionales (entradas +-1), luego V_(1^6) ES el modelo estandar de Segre'
      u' sobre Z, no una forma de el',
      all(all(c in (1, -1) for c in v) for v in sing))

# ---------------------------------------------------------------- (b) CAYLEY
x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
u0, u1, u2, u3 = sp.symbols('u0 u1 u2 u3')
VX = [x1, x2, x3, x4]
x5 = -(2 * x1 + x2 + x3 + x4)
F = sp.expand(2 * x1 ** 3 + x2 ** 3 + x3 ** 3 + x4 ** 3 + x5 ** 3)

# control del control: F es la cubica de V_(2,1^4), reconstruida de la definicion
comprueba(u'b0. control: F se anula donde la definicion de V_(2,1^4) dice',
      sp.simplify(F.subs({x1: 1, x2: -1, x3: -1, x4: -1})) == 0)

# los nodos, CALCULADOS
J = [sp.expand(sp.diff(F, v)) for v in VX]
nodos = set()
for v in VX:
    otras = [w for w in VX if w is not v]
    ss = sp.solve([j.subs(v, 1) for j in J], otras, dict=True)
    for s in ss:
        if any(getattr(val, 'free_symbols', set()) for val in s.values()):
            continue
        pt = tuple(sp.Integer(1) if w is v else sp.nsimplify(s.get(w, 0)) for w in VX)
        for c in pt:
            if c != 0:
                pt = tuple(sp.nsimplify(c2 / c) for c2 in pt)
                break
        nodos.add(pt)

nodos = sorted(nodos, key=str)
comprueba(u'b1. V_(2,1^4) tiene EXACTAMENTE cuatro nodos, calculados del jacobiano',
      len(nodos) == 4, u'%s' % (nodos,))
comprueba(u'b2. y los cuatro son RACIONALES --- que es lo que decide toda la cuestion',
      all(c.is_rational for P in nodos for c in P))

M = sp.Matrix([[nodos[j][i] for j in range(4)] for i in range(4)])
comprueba(u'b3. la matriz de los nodos esta en GL_4(Q): det = %s' % M.det(), M.det() != 0)
comprueba(u'   ... y sus entradas son +-1, luego esta en GL_4(Z[1/2])',
      all(abs(c) == 1 for c in M) and abs(M.det()) == 8)

G = sp.expand(F.subs(dict(zip(VX, list(M * sp.Matrix([u0, u1, u2, u3]))))))
CAY = u0 * u1 * u2 + u0 * u1 * u3 + u0 * u2 * u3 + u1 * u2 * u3
lam = sp.simplify(sp.cancel(G / CAY))
comprueba(u'b4. F(Mu) = 24 * (u0u1u2 + u0u1u3 + u0u2u3 + u1u2u3), IDENTIDAD simbolica',
      sp.expand(G - 24 * CAY) == 0, u'cociente = %s' % lam)
comprueba(u'   ... y no hace falta reescalar ninguna variable: el factor es una constante',
      lam.is_number and lam == 24)

coord_nodos = 0
for e in sp.eye(4).tolist():
    d = dict(zip((u0, u1, u2, u3), e))
    if (sp.simplify(G.subs(d)) == 0
            and all(sp.simplify(sp.diff(G, v).subs(d)) == 0 for v in (u0, u1, u2, u3))):
        coord_nodos += 1
comprueba(u'b5. en las coordenadas nuevas los cuatro puntos COORDENADOS son los nodos',
      coord_nodos == 4)

# CONTROL NEGATIVO: cambiar un nodo por otro punto racional tiene que romperlo
Mmal = M.copy()
Mmal[:, 3] = sp.Matrix([1, 1, 1, 1])
roto = True
if Mmal.det() != 0:
    Gm = sp.expand(F.subs(dict(zip(VX, list(Mmal * sp.Matrix([u0, u1, u2, u3]))))))
    roto = sp.simplify(sp.cancel(Gm / CAY)).is_number is not True
comprueba(u'b6. CONTROL NEGATIVO: sustituir un nodo por (1,1,1,1) ROMPE la forma de Cayley', roto)

# ---------------------------------------------------------------- (c) la hipotesis, afilada
comprueba(u'c1. la hipotesis char no divide 6 es AFILADA: det M = -8 se anula en p = 2',
      M.det() % 2 == 0 and M.det() % 3 != 0)
comprueba(u'   ... y la constante 24 = 2^3 * 3 se anula en p = 2 y en p = 3, y en ningun otro',
      sorted(sp.factorint(24)) == [2, 3])

# ---------------------------------------------------------------- el puente del tamiz
z = sp.symbols('z1:6')
seg_cortada = sp.expand(Cub.subs({y[0]: z[0], y[1]: z[0], y[2]: z[1],
                                  y[3]: z[2], y[4]: z[3], y[5]: z[4]}))
v22 = sp.expand(2 * z[0] ** 3 + z[1] ** 3 + z[2] ** 3 + z[3] ** 3 + z[4] ** 3)
comprueba(u'd1. la seccion y_1 = y_2 de Segre ES V_(2,1^4): las dos ecuaciones, identicas',
      sp.expand(seg_cortada - v22) == 0
      and sp.expand(Lin.subs({y[0]: z[0], y[1]: z[0], y[2]: z[1], y[3]: z[2],
                              y[4]: z[3], y[5]: z[4]})
                    - (2 * z[0] + z[1] + z[2] + z[3] + z[4])) == 0)

L.append(u'')
L.append(u'  Los cuatro nodos de V_(2,1^4):')
for P in nodos:
    L.append(u'    %s' % (tuple(int(c) for c in P),))
L.append(u'')
L.append(u'  M = columnas de los nodos, det = %s' % M.det())
for fila in M.tolist():
    L.append(u'    [%s]' % u'  '.join(u'%2d' % int(c) for c in fila))

L.append(u'')
L.append(u'VEREDICTO')
L.append(u'  SEGRE: V_(1^6) es el modelo estandar sobre Z, no una forma de el. Nada que probar.')
L.append(u'  CAYLEY: equivalencia proyectiva SOBRE EL CUERPO PRIMO, con la matriz explicita de')
L.append(u'  entradas +-1 arriba, valida en toda caracteristica que no divida a 6.  Lo que la')
L.append(u'  hace posible es que los cuatro nodos sean racionales; el control negativo comprueba que no')
L.append(u'  vale cualquier cuarta columna.  Luego los recuentos transfieren sin torcer.')

print(u'\n'.join(L))
print(u'')
print(u'=' * 72)
print(u'VERIFICACION SEGRE Y CAYLEY SOBRE EL CUERPO BASE:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 72)
sys.exit(1 if MAL else 0)
