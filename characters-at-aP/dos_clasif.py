# -*- coding: utf-8 -*-
"""Las DOS clasificaciones no son la misma, y esa es la aportacion:

  (I)  valores ENTEROS      :  q|2m,  q|(2m+1),  q|(2m+2),  o  q=6 y m=1 (3)
  (U)  cancelacion, |U|=1   :  q|m,   q|(2m+1),  q|(2m+2),  o  q=6 y m=1 (6)

(U) esta contenida en (I).  La diferencia son los puntos donde TODOS los caracteres son enteros
y aun asi hace falta un corrector.  Ahi es donde el articulo actual no llega.
"""
import os
import sys
from math import gcd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # antes: un scratchpad temporal (9-sep-2026)
from galois import chi_lim, galois_estable, condicion_I


def condicion_U(m, q):
    a = (q - 1) // 2
    u, v = divmod(m, q)
    return (v in (0, a, q - 1)) or (q == 6 and v == 1)


print("=" * 92)
print("(U) esta contenida en (I)?  y cuantos puntos hay en (I) y no en (U)?")
print("=" * 92)
dentro = fuera = solo_I = solo_U = 0
ejem = []
for q in range(2, 41):
    for m in range(1, 121):
        I, U = condicion_I(m, q), condicion_U(m, q)
        if U and I:
            dentro += 1
        if U and not I:
            solo_U += 1
        if I and not U:
            solo_I += 1
            if len(ejem) < 10:
                ejem.append((m, q))
print("   en (U) y en (I): %d      en (U) pero NO en (I): %d   <- deberia ser 0" % (dentro, solo_U))
print("   en (I) pero NO en (U): %d   <- enteros, pero con corrector" % solo_I)
print("   ejemplos:", ", ".join("m=%d q=%d" % e for e in ejem))

print()
print("=" * 92)
print("EL EJEMPLO NUEVO:  Sp(8), q=6, lambda=(2,1,0,0)  ->  dice que vale -6")
print("=" * 92)
m, q, lam = 4, 6, (2, 1, 0, 0)
print("   2m=%d mod q=%d | 2m+1=%d mod q=%d | 2m+2=%d mod q=%d | q=6 y m mod 3=%d"
      % (2 * m, 2 * m % q, 2 * m + 1, (2 * m + 1) % q, 2 * m + 2, (2 * m + 2) % q, m % 3))
print("   en (I): %s     en (U): %s" % (condicion_I(m, q), condicion_U(m, q)))
v = chi_lim(lam, q)
print("   chi_(2,1,0,0)(g_{1/6}) = %s" % v)
print("   redondeado: %d" % round(float(v)))
print("   y el denominador de a_P aqui seria 2m+2 = %d, asi que q=6 NO es potencia de a_P" % (2 * m + 2))

print()
print("   y el que separa integridad de cancelacion:  Sp(4), q=4, lambda=(1,0)")
print("   en (I): %s (porque 4 | 2m=4)   en (U): %s"
      % (condicion_I(2, 4), condicion_U(2, 4)))
print("   chi_(1,0)(g_{1/4}) = %s   -> entero, y sin embargo |U| = 2" % chi_lim((1, 0), 4))
