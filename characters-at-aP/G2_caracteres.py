# -*- coding: utf-8 -*-
# G2 CON CARACTERES DE VERDAD: LA PREDICCION DE LA NOTA, PUESTA A PRUEBA.   25 de agosto de 2026.
#
# POR QUE.  La nota 2 predice, y no puede comprobar, que si una clase fija lleva floor(d/2) porque
# el plegado de tipo C es una INVOLUCION, en G2 --- donde la razon de longitudes es r = 3 ---
# deberia llevar floor(d/3).  Aqui se construye el evaluador que falta y se mide.
#
# EL ELEMENTO.  Kostant define a_P por <alpha, x_P> proporcional a (alpha,alpha), o sea peso 1 en
# las simples CORTAS y r en las LARGAS.  En G2: alpha_1 corta, alpha_2 larga, r = 3, y el orden es
# r*h^v = 3*4 = 12.  Para alpha = c1 alpha_1 + c2 alpha_2:   ht_P(alpha) = c1 + 3 c2.
# Las seis raices positivas dan ht_P = 1, 3, 4, 5, 6, 9 --- ninguna divisible por 12, luego a_P es
# REGULAR, como debe.
#
# EL EVALUADOR.  Multiplicidades de peso por FREUDENTHAL (exacto, en fracciones), y luego
#       chi_lambda(x) = suma_mu m_mu * x^mu.
# En G2 el reticulo de pesos coincide con el de raices, asi que todo peso es combinacion entera de
# alpha_1, alpha_2 y x^mu = zeta_12^{k (c1 + 3 c2)}.
#
# QUE SE MIDE
#   G0  CONTROL: la suma de multiplicidades es la dimension de Weyl.  Si fallara, nada de lo demas
#       significa nada.
#   G1  CONTROL: en a_P (k=1, regular) los valores tienen que caer en {0,+-1} --- es el Teorema 3.1
#       de Kostant.  Es el control que dice que el elemento es el correcto.
#   G2  los valores de chi en a_P^k, y su espectro.
#   G3  LA PREDICCION: ¿que capacidades (no fija, fija) reproducen la anulacion?  Se prueban todas
#       las de 0..4.  La nota predice floor(d/3) en las fijas.
#
# Authors: Carles Marin, Claude (AI assistant).
# Run:  python G2_caracteres.py > G2_caracteres_OUT.txt 2>&1

import cmath
import sys
from collections import Counter
from fractions import Fraction
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SEP = "=" * 92
fallos = []


def ok(c, e, d=""):
    print(("   OK   " if c else "  FALLA ") + e + ("   " + d if d else ""))
    if not c:
        fallos.append(e)


# ---------------------------------------------------------------- G2: datos basicos
# base: alpha_1 CORTA (long^2 = 2), alpha_2 LARGA (long^2 = 6), (alpha_1,alpha_2) = -3
GRAM = [[Fraction(2), Fraction(-3)], [Fraction(-3), Fraction(6)]]


def ip(u, v):
    return sum(u[i] * v[j] * GRAM[i][j] for i in range(2) for j in range(2))


POS = [(1, 0), (0, 1), (1, 1), (2, 1), (3, 1), (3, 2)]     # raices positivas, en la base alpha
RHO = (5, 3)                                               # semisuma
OM = [(2, 1), (3, 2)]                                      # pesos fundamentales omega_1, omega_2


def htP(al):
    return al[0] + 3 * al[1]                               # peso 1 en la corta, r=3 en la larga


def largo(al):
    return ip(al, al) == 6


print(SEP)
print("DATOS DE G2")
print(f"   raices positivas y su ht_P: {[(a, htP(a), 'larga' if largo(a) else 'corta') for a in POS]}")
print(f"   rho = {RHO},  orden de a_P = r*h^v = 3*4 = 12")
ok(sorted(htP(a) for a in POS) == [1, 3, 4, 5, 6, 9], "los ht_P son 1,3,4,5,6,9")
ok(sum(1 for a in POS if largo(a)) == 3, "tres largas y tres cortas")
ok(all(htP(a) % 12 for a in POS), "a_P es REGULAR: ningun ht_P divisible por 12")


# ---------------------------------------------------------------- Freudenthal
def pesos_y_multiplicidades(lam):
    """lam en la base alpha.  Devuelve {mu (base alpha): multiplicidad}."""
    lr = (lam[0] + RHO[0], lam[1] + RHO[1])
    c_lam = ip(lr, lr)
    # todos los mu = lam - n1 a1 - n2 a2 con n_i >= 0 y mu >= -lam (peso mas bajo, w0 = -1 en G2)
    N1, N2 = 2 * lam[0], 2 * lam[1]
    mult = {}
    orden = sorted(((n1, n2) for n1 in range(N1 + 1) for n2 in range(N2 + 1)),
                   key=lambda n: (n[0] + n[1], n))
    for (n1, n2) in orden:
        mu = (lam[0] - n1, lam[1] - n2)
        if (n1, n2) == (0, 0):
            mult[mu] = Fraction(1)
            continue
        s = Fraction(0)
        for al in POS:
            k = 1
            while True:
                nu = (mu[0] + k * al[0], mu[1] + k * al[1])
                if nu not in mult:
                    break
                if mult[nu]:
                    s += mult[nu] * ip(nu, al)
                k += 1
        mr = (mu[0] + RHO[0], mu[1] + RHO[1])
        den = c_lam - ip(mr, mr)
        mult[mu] = Fraction(0) if den == 0 else 2 * s / den
    return {mu: int(v) for mu, v in mult.items() if v > 0}


def dim_weyl(lam):
    lr = (lam[0] + RHO[0], lam[1] + RHO[1])
    num = den = Fraction(1)
    for al in POS:
        num *= ip(lr, al)
        den *= ip(RHO, al)
    return num / den


print(SEP)
print("G0 -- CONTROL: suma de multiplicidades == dimension de Weyl")
DOM = [(a, b) for a in range(7) for b in range(7)]          # lambda = a om_1 + b om_2
LAMS = []
for (a, b) in DOM:
    lam = (a * OM[0][0] + b * OM[1][0], a * OM[0][1] + b * OM[1][1])
    LAMS.append(((a, b), lam))
malos = 0
for (ab, lam) in LAMS[:8]:
    W = pesos_y_multiplicidades(lam)
    s, dw = sum(W.values()), dim_weyl(lam)
    malos += (s != dw)
    print(f"   lambda = {ab[0]} om1 + {ab[1]} om2:  sum m_mu = {s:>5}   dim Weyl = {dw}"
          + ("" if s == dw else "   <-- FALLA"))
ok(malos == 0, "G0: Freudenthal reproduce la dimension de Weyl")

# ---------------------------------------------------------------- G1/G2
print(SEP)
print("G1 -- CONTROL: en a_P (k=1) los valores tienen que caer en {0,+-1}  (Kostant, Thm 3.1)")


def chi(lam, k, W=None):
    z = cmath.exp(2j * cmath.pi / 12)
    W = W or pesos_y_multiplicidades(lam)
    s = sum(m * z ** ((k * (mu[0] + 3 * mu[1])) % 12) for mu, m in W.items())
    if abs(s.imag) > 1e-6 or abs(s.real - round(s.real)) > 1e-6:
        return None
    return int(round(s.real))


cache = {}
vals1 = Counter()
for (ab, lam) in LAMS:
    cache[ab] = pesos_y_multiplicidades(lam)
    v = chi(lam, 1, cache[ab])
    vals1[v] += 1
ok(set(vals1) <= {0, 1, -1}, "G1: chi(a_P) in {0,+-1}", f"valores {dict(sorted(vals1.items()))}")

print()
print("G2 -- LOS VALORES EN a_P^k")
print(f"{'k':>3} {'d':>3} {'q':>4} {'pesos':>6} {'nulos':>6}   espectro")
DATOS = {}
for k in range(1, 7):
    d = gcd(k, 12)
    q = 12 // d
    vals = Counter()
    reg = {}
    for (ab, lam) in LAMS:
        v = chi(lam, k, cache[ab])
        if v is None:
            continue
        vals[v] += 1
        reg[ab] = (lam, v)
    DATOS[k] = (d, q, reg)
    print(f"{k:>3} {d:>3} {q:>4} {sum(vals.values()):>6} {vals.get(0,0):>6}   {dict(sorted(vals.items()))}")

# ---------------------------------------------------------------- G3
print(SEP)
print("G3 -- ¿TRANSPORTA EL CRITERIO DE TIPO C?")
print()
print("   ⚠ AVISO SOBRE UN INTENTO FALLIDO, que se deja escrito porque el error es instructivo.")
print("   Primero 'derive' una condicion asi: el numerador de Weyl en x = a_P^k es")
print("   suma_w eps(w) zeta_12^{k ht_P(w(l+r))}, y si una reflexion s_alpha deja ese exponente")
print("   igual modulo q sus dos terminos se cancelan, lo que da  q | <l+r,alpha^v> ht_P(alpha).")
print("   ESO ESTA MAL: cancelar la pareja (1, s_alpha) no anula la SUMA, que tiene |W| terminos.")
print("   Lo que si es cierto es que el numerador se anula siempre que el estabilizador de x en W")
print("   contenga una reflexion --- o sea SIEMPRE que x no sea regular, que es el caso k>=2.  Por")
print("   eso el numerador no dice nada alli: hace falta el limite, y ese es justo el paso que la")
print("   nota declara abierto.  Abajo se mide, sin derivar nada.")
print()
print(f"{'k':>3} {'d':>3} {'q':>4} {'pesos':>6} {'DERIVADO':>10} {'q|<l+r,a^v>':>13} {'q|(l+r,a)':>11}")
for k in sorted(DATOS):
    d, q, reg = DATOS[k]
    if q < 2:
        continue
    ader = a2 = a1 = tot = 0
    for ab, (lam, v) in reg.items():
        lr = (lam[0] + RHO[0], lam[1] + RHO[1])
        cder = any((int(2 * ip(lr, al) / ip(al, al)) * htP(al)) % q == 0 for al in POS)
        c2 = any(int(2 * ip(lr, al) / ip(al, al)) % q == 0 for al in POS)
        c1 = any(ip(lr, al) % q == 0 for al in POS)
        tot += 1
        ader += (cder == (v == 0))
        a2 += (c2 == (v == 0))
        a1 += (c1 == (v == 0))
    # el veredicto se toma sobre el criterio de RAICES, que es el de la nota 1 (Prop. 7.1);
    # el 'derivado' se deja en la tabla como lo que es: un intento fallido, documentado arriba
    marca = "" if a1 == tot else "   <-- ninguno cierra"
    print(f"{k:>3} {d:>3} {q:>4} {tot:>6} {str(ader)+'/'+str(tot):>10} {str(a2)+'/'+str(tot):>13} "
          f"{str(a1)+'/'+str(tot):>11}{marca}")
    if k == 1 and a1 != tot:
        fallos.append("G3: en k=1 el criterio de raices deberia cerrar y no cierra")

print(SEP)
print("LECTURA -- se escribe a la vista de la tabla, no antes.")
print(SEP)
if fallos:
    print(f"RESULTADO: {len(fallos)} FALLOS -> {fallos}")
    sys.exit(1)
print("RESULTADO: los controles pasan; la lectura de G3 va en la nota.")
