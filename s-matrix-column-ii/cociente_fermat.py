# -*- coding: utf-8 -*-
r"""cociente_fermat.py -- el recuento de semisistemas detecta los primos de Wieferich.

DE DONDE SALE.  N(p) = #{semisistemas de F_p de suma nula} = (2^m + eps (p-1))/p con m = (p-1)/2 y
eps = (2|p) (semisistemas_cuenta.py, comprobado hasta 199).  Como eps^2 = 1,

        2^{p-1} - 1 = (2^m - eps)(2^m + eps),

y modulo p se tiene 2^m = eps, luego 2^m + eps = 2 eps.  Dividiendo,

        N(p) = eps * ( 1 + q_2(p) / 2 )   (mod p),      q_2(p) = (2^{p-1} - 1)/p.

CONSECUENCIA.  p es de WIEFERICH (p^2 divide 2^{p-1} - 1, o sea q_2(p) = 0 mod p) si y solo si

        N(p) = (2|p)   (mod p).

Es decir: un recuento combinatorio de semisistemas -- el objeto del lema de Gauss -- reconoce los
primos de Wieferich.  Los dos conocidos por debajo de 10^15 son 1093 y 3511.

ESTE GUION comprueba la identidad en todos los primos de un rango y, aparte, en 1093 y 3511.

    python cociente_fermat.py [PMAX]

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(AQUI), "note_P3a", "release_tool"))
from P3a import es_primo                                            # noqa: E402

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 500


def eps(p):
    return 1 if pow(2, (p - 1) // 2, p) == 1 else -1


def N(p):
    """el recuento exacto, como entero grande."""
    e = eps(p)
    t = 2 ** ((p - 1) // 2) + e * (p - 1)
    assert t % p == 0, p
    return t // p


def q2(p):
    t = (pow(2, p - 1, p * p) - 1) % (p * p)
    assert t % p == 0, p
    return (t // p) % p


fallos = casos = 0
wieferich = []
for p in range(5, PMAX + 1):
    if not es_primo(p):
        continue
    casos += 1
    e, q = eps(p), q2(p)
    pred = (e * (1 + q * pow(2, p - 2, p))) % p
    if N(p) % p != pred:
        fallos += 1
        print("FALLA en p = %d" % p)
    if q == 0:
        wieferich.append(p)

print("primos con la identidad N(p) = (2|p)(1 + q_2(p)/2) mod p comprobada: %d ; fallos: %d"
      % (casos, fallos))
print("primos de Wieferich hallados por debajo de %d: %s" % (PMAX, wieferich or "ninguno"))

print("")
print("los dos de Wieferich conocidos, con N(p) mod p calculado aparte:")
for p in (1093, 3511):
    e = eps(p)
    t = (pow(2, (p - 1) // 2, p * p) - e) % (p * p)
    npmod = ((t // p) + e) % p
    print("   p=%-6d (2|p)=%-3d  N(p) mod p = %-6d  (2|p) mod p = %-6d  %s"
          % (p, e, npmod, e % p, "COINCIDE: criterio de Wieferich" if npmod == e % p else "NO"))

print("")
print("un control de no-Wieferich, donde el criterio debe FALLAR:")
for p in (101, 2003, 5003):
    e = eps(p)
    t = (pow(2, (p - 1) // 2, p * p) - e) % (p * p)
    npmod = ((t // p) + e) % p
    print("   p=%-6d N(p) mod p = %-6d  (2|p) mod p = %-6d  %s"
          % (p, npmod, e % p, "distintos, como debe ser" if npmod != e % p else "IGUALES: mal"))
