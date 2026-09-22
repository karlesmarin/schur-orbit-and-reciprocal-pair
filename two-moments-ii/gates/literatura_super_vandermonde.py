# -*- coding: utf-8 -*-
u"""literatura_super_vandermonde.py -- el caso k = 4, n = 3 es un super-Vandermonde set de tamanyo 3.

LA DEFINICION (Sziklai-Takats, Finite Fields Appl. 14 (2008) 1056-1067; verbatim del
capitulo 1 de la tesis de Takats, ELTE 2014, que ES ese articulo):

  "Let w = w_S be the smallest positive integer k such that pi_k != 0 if such a k exists,
   otherwise w = infinity."
  "Definition 1.5.  Let 1 < t < q.  We say that T = {y_1, ..., y_t} in GF(q) is a
   super-Vandermonde set, if pi_k = sum_i y_i^k = 0 for all 1 <= k <= t-1."

Y EL TEOREMA 1.12: "Suppose that T in GF(q) is a super-Vandermonde set of size |T| < p.
Then T is a (transform of a) multiplicative subgroup."

PARA t = 3 Y q = p PRIMO eso dice: f(Y) = Y^3 - b_0 y 3 | p-1, o sea T = a{1, w, w^2} con w
raiz cubica primitiva.  Que es EXACTAMENTE la familia del caso k = 4.

La equivalencia con el caso k = 4 pasa por Newton: con e_1 = 0, p_4 = 2 e_2^2, luego
{p_1 = p_4 = 0} <=> {e_1 = e_2 = 0} <=> {p_1 = p_2 = 0}.  Hace falta p distinto de 2 --- en
caracteristica 2 la identidad es vacua --- y p distinto de 3, donde no hay sV-sets de tamanyo
3 porque p no puede dividir a t.

QUE CUBRE ESA LITERATURA Y QUE NO
  CUBIERTO: la estructura del caso k = 4, n = 3.  Y el recuento (p-1)/3 es corolario de una
            linea del Teorema 1.12.
  NO CUBIERTO: los momentos prescritos NO CONSECUTIVOS.  w_T es por definicion el primer
            momento no nulo, asi que la condicion es siempre el segmento inicial completo.
            La familia de la nota --- p_1 = 0 y p_k = 0 con los intermedios LIBRES --- no cabe.
  NO CUBIERTO: los recuentos C_k(3,p) como funcion de k y p.  Esa literatura es estructural.

Se comprueba por fuerza bruta sobre todas las ternas de F_p, 5 <= p < 120.  El testigo
negativo (seccion 4) son ternas con p_1 = p_6 = 0 y p_2 != 0, que ninguna definicion de
Vandermonde set admite; si no existieran, la seccion 4 fallaria.

chi (el simbolo de Legendre) y es_primo van definidos aqui mismo: el guion no importa nada
fuera de la biblioteca estandar.

Uso:  python literatura_super_vandermonde.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
from __future__ import print_function

import itertools
import sys


def es_primo(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def chi(a, p):
    """El simbolo de Legendre (a|p)."""
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

OK = 0
MAL = 0
L = []


def comprueba(et, cond, det=u''):
    global OK, MAL
    if cond:
        OK += 1
        L.append(u'  ok   %s%s' % (et, (u'   ' + det) if det else u''))
    else:
        MAL += 1
        L.append(u'  MAL  %s%s' % (et, (u'   ' + det) if det else u''))


def sec(t):
    L.append(u'')
    L.append(t)
    L.append(u'-' * len(t))


PS = [q for q in range(5, 120) if es_primo(q)]


def ternas(p, cond):
    u"""Las ternas de F_p distintas que cumplen cond(p1, p2, p3, p4)."""
    out = []
    for T in itertools.combinations(range(p), 3):
        p1 = sum(T) % p
        p2 = sum(t * t for t in T) % p
        p3 = sum(pow(t, 3, p) for t in T) % p
        p4 = sum(pow(t, 4, p) for t in T) % p
        if cond(p1, p2, p3, p4):
            out.append(T)
    return out


# ---------------------------------------------------------------------------
sec(u'1.  El caso k = 4 es el super-Vandermonde set de tamanyo 3')
# ---------------------------------------------------------------------------
malos = []
for p in PS:
    fam_k4 = set(ternas(p, lambda a, b, c, d: a == 0 and d == 0))
    sv = set(ternas(p, lambda a, b, c, d: a == 0 and b == 0 and c != 0))
    # la familia k = 4 incluye la terna degenerada con p_3 = 0 si la hubiera
    fam_k4_no_deg = set(T for T in fam_k4
                         if sum(pow(t, 3, p) for t in T) % p != 0)
    if fam_k4_no_deg != sv:
        malos.append((p, len(fam_k4_no_deg), len(sv)))
comprueba(u'{p_1 = p_4 = 0, p_3 != 0} = {p_1 = p_2 = 0, p_3 != 0} en %d primos' % len(PS),
      not malos, u'fallan: %s' % (malos[:3],))
L.append(u'      (por Newton con e_1 = 0: p_4 = 2 e_2^2, luego p_4 = 0 <=> e_2 = 0 <=> p_2 = 0)')

# ---------------------------------------------------------------------------
sec(u'2.  Y el Teorema 1.12 las describe: cosets de las raices cubicas')
# ---------------------------------------------------------------------------
malos = []
for p in PS:
    sv = ternas(p, lambda a, b, c, d: a == 0 and b == 0 and c != 0)
    # las raices cubicas de la unidad, si existen
    omegas = [w for w in range(2, p) if pow(w, 3, p) == 1]
    esperado = set()
    if omegas:
        w = omegas[0]
        for a in range(1, p):
            esperado.add(tuple(sorted({a % p, (a * w) % p, (a * w * w) % p})))
    if set(sv) != esperado:
        malos.append((p, len(sv), len(esperado)))
comprueba(u'toda terna sV es a{1, w, w^2} con w raiz cubica primitiva', not malos,
      u'fallan: %s' % (malos[:3],))

# ---------------------------------------------------------------------------
sec(u'3.  Y el recuento es (p-1)/3 si 3 | p-1, cero si no')
# ---------------------------------------------------------------------------
malos = []
for p in PS:
    n = len(ternas(p, lambda a, b, c, d: a == 0 and b == 0 and c != 0))
    esp = (p - 1) // 3 if p % 3 == 1 else 0
    if n != esp:
        malos.append((p, n, esp))
comprueba(u'#sV(3, p) = (p-1)/3 si p = 1 mod 3, y 0 si no', not malos,
      u'fallan: %s' % (malos[:3],))
comprueba(u'y eso es la formula C_4(3,p) = (p-1)(1+chi_p(-3))/6',
      all(((p - 1) * (1 + chi(-3, p)) // 6) == ((p - 1) // 3 if p % 3 == 1 else 0)
          for p in PS))
L.append(u'      => el recuento es corolario de UNA LINEA del Teorema 1.12.')

# ---------------------------------------------------------------------------
sec(u'4.  Lo que esa literatura NO cubre: los momentos no consecutivos')
# ---------------------------------------------------------------------------
# w_T es el PRIMER momento no nulo: la condicion es el segmento inicial completo.
# La familia de la nota deja los intermedios LIBRES.  Testigo: k = 6, donde p_2 no se anula.
p = 19
sueltas = ternas(p, lambda a, b, c, d: a == 0
                 and (sum(pow(t, 6, p) for t in (0, 0, 0)) == 0 or True))
con6 = ternas(p, lambda a, b, c, d: a == 0)
con6 = [T for T in con6 if sum(pow(t, 6, p) for t in T) % p == 0]
libres = [T for T in con6 if sum(t * t for t in T) % p != 0]
comprueba(u'en p = 19 hay ternas con p_1 = p_6 = 0 y p_2 DISTINTO de 0', len(libres) > 0,
      u'%d de %d' % (len(libres), len(con6)))
L.append(u'      => esas NO son Vandermonde sets de ningun tipo: su w_T es 2, no 3.')
L.append(u'         La familia de momentos no consecutivos no esta en esa literatura.')

L.append(u'')
L.append(u'VEREDICTO')
L.append(u'  CONOCIDO: la estructura del caso k = 4, n = 3 --- es el super-Vandermonde set de')
L.append(u'          tamanyo 3 de Sziklai-Takats 2008, Teorema 1.12, y el recuento es')
L.append(u'          corolario inmediato.  La nota lo cita como recuperacion.')
L.append(u'  NO CUBIERTO: los momentos prescritos NO CONSECUTIVOS, que es el objeto de la')
L.append(u'          nota; y los recuentos C_k(3,p) como funcion de k.')

print(u'\n'.join(L))
print(u'')
print(u'=' * 70)
print(u'VERIFICACION SUPER-VANDERMONDE:  %d ok, %d MAL' % (OK, MAL))
print(u'=' * 70)
sys.exit(1 if MAL else 0)
