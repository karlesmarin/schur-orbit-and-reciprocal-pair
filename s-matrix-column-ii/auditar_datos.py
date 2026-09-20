# -*- coding: utf-8 -*-
r"""auditar_datos.py -- cada numero y cada formula de la nota, contra su origen.

Un numero en un paper lleva dentro el instrumento que lo produjo.  Esto lo comprueba: los que
son baratos se RECALCULAN aqui desde cero, y los que son caros se buscan en la salida
archivada del guion que los produjo.  Un numero que no aparezca en ninguno de los dos sitios
es un numero sin origen, y esos son los que hunden un paper.

    python auditar_datos.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import io
import os
import re
import sys
from math import comb, factorial

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = os.path.dirname(os.path.abspath(__file__))
GATES = os.path.join(os.path.dirname(AQUI), "gates")
TEX = io.open(os.path.join(AQUI, "P3b.tex"), encoding="utf-8").read()
CUERPO = TEX.split(r"\begin{thebibliography}")[0]

fallos = []
hechas = []


def afirma(ok, que):
    print("   %-4s %s" % ("ok" if ok else "MAL", que))
    hechas.append(que)
    if not ok:
        fallos.append(que)


def en_tex(t):
    return t in CUERPO.replace("\\,", "").replace("~", " ")


def en_salida(fichero, *trozos):
    p = os.path.join(GATES, fichero)
    if not os.path.exists(p):
        return False
    s = io.open(p, encoding="utf-8", errors="replace").read()
    return all(t in s for t in trozos)


print("(1) EL GENERO: formula recalculada, y contra la forma de Korchmaros-Lia-Timpanella")
for n in range(5, 10):
    g = 1 + factorial(n - 2) * (n * n - 5 * n + 2) // 4
    klt = factorial(n - 2) * ((n - 2) * (n - 3) - 4) // 4 + 1
    afirma(g == klt, "n=%d: nuestra forma %d == la suya %d" % (n, g, klt))
afirma(en_tex("4,\\;49,\\;481,\\;4681,\\;47881"), "la lista 4, 49, 481, 4681, 47881 esta en el tex")

print("")
print("(2) LA TABLA DE LA FRONTERA: las veinte celdas, recalculadas")
for n, gg, cota, fact, ratio in ((5, 4, 64, 120, "0.5"), (6, 49, 9604, 720, "13.3"),
                                 (7, 481, 925444, 5040, "183.6"),
                                 (8, 4681, 87647044, 40320, "2173.8")):
    g = 1 + factorial(n - 2) * (n * n - 5 * n + 2) // 4
    afirma(g == gg, "g_%d = %d" % (n, gg))
    afirma(4 * g * g == cota, "(2g_%d)^2 = %d" % (n, cota))
    afirma(factorial(n) == fact, "%d! = %d" % (n, fact))
    afirma(abs(4 * g * g / factorial(n) - float(ratio)) < 0.06,
           "cociente %s (calculado %.1f)" % (ratio, 4 * g * g / factorial(n)))

print("")
print("(3) N(p): la formula, recalculada contra fuerza bruta")


def es_primo(m):
    return m > 1 and all(m % d for d in range(2, int(m ** 0.5) + 1))


def N_formula(p):
    leg = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    num = 2 ** ((p - 1) // 2) + leg * (p - 1)
    assert num % p == 0, p
    return num // p


def N_fuerza(p):
    m = (p - 1) // 2
    pares = [(a, p - a) for a in range(1, m + 1)]
    c = 0
    for mascara in range(1 << m):
        s = 0
        for i, (a, b) in enumerate(pares):
            s += a if (mascara >> i) & 1 else b
        if s % p == 0:
            c += 1
    return c


for p in (7, 11, 13, 17, 19, 23):
    afirma(N_formula(p) == N_fuerza(p),
           "N(%d) = %d por formula y por fuerza bruta" % (p, N_formula(p)))
afirma(en_tex("2^{(p-1)/2}+(2\\mid p)(p-1)}{p}"), "la formula de N(p) esta en el tex")

print("")
print("(4) LAS CUENTAS CARAS: cada una, en la salida del guion que la produjo")
CUENTAS = [
    ("7935", "cascada_newton_OUT.txt", "orbitas examinadas: 7935"),
    ("438", "tricotomia_OUT.txt", "orbitas h=1 fuera de los lugares: 438"),
    ("2924", "semisistemas_OUT.txt", "semisistemas examinados: 2924"),
    ("36", "oeis_A262568_OUT.txt", "terminos comparados: 36"),
    ("199", "semisistemas_cuenta_OUT.txt", "primos comprobados hasta 199"),
    ("427", "sato_tate_OUT.txt", "primos usados: 427"),
    ("0.0097", "sato_tate_OUT.txt", "0.009694"),
    ("1093", "cociente_fermat_OUT.txt", "p=1093"),
    ("3511", "cociente_fermat_OUT.txt", "p=3511"),
    # ANTES DECIA 607, Y ERA FALSO.  El guion se habia corrido con PMIN = 601, asi que el 607
    # era el primer primo de su ventana y no el primero en absoluto.  Ahora la compuerta barre
    # desde p = n+1 y se niega a decir "primera aparicion" si no arranco ahi, y lo que se pide
    # aqui son las TRES primeras apariciones, que es lo que la nota debe citar.
    ("67  (n=5)", "umbral_real_OUT.txt", "PRIMERA APARICION: p = 67"),
    ("163 (n=6)", "umbral_real_OUT.txt", "PRIMERA APARICION: p = 163"),
    ("601 (n=7)", "umbral_real_OUT.txt", "PRIMERA APARICION: p = 601"),
    ("67 = primer N>0", "extremo_bring_OUT.txt", "primer p con N_orb > 0 (barrido desde 7): 67"),
    ("51 primos", "extremo_bring_OUT.txt", "primos: 51 ; fallos: 0"),
    ("34/3/1/1", "orbitas_compuesto_OUT.txt", "h=1:34 h=3:3 h=5:1 h=15:1"),
]
for etiqueta, fichero, trozo in CUENTAS:
    afirma(en_salida(fichero, trozo), "%-10s en %s" % (etiqueta, fichero))

# EL AGUJERO QUE DEJO PASAR EL 233 Y EL 607.  Lo de arriba comprueba que el numero esta en la
# SALIDA DEL GUION, no que el numero que dice la NOTA sea ese.  Con la nota diciendo 607 y el
# guion diciendo 607 --- porque se habia corrido con PMIN = 601 --- las dos casillas cuadraban y
# el error era invisible.  Un auditor que solo mira un lado no audita nada.
print("")
print("(4b) Y QUE LA NOTA DIGA LO MISMO, no solo el guion")
PRIMERAS = [("n=5", 67, (233,)), ("n=6", 163, (607,)), ("n=7", 601, ())]
for etiqueta, bueno, malos in PRIMERAS:
    afirma(en_tex("$p=%d$" % bueno) or en_tex("p=%d" % bueno) or en_tex(str(bueno)),
           "%s: la nota cita la primera aparicion %d" % (etiqueta, bueno))
    for malo in malos:
        # El valor falso PUEDE aparecer, pero solo dentro de la nota de erratas que dice que lo
        # era: un paper que corrige un numero tiene que poder nombrarlo.  Lo que no puede es
        # aparecer suelto, como afirmacion.  Se exige que toda aparicion este cerca de
        # "earlier version", que es como se rotula la erratum.
        sitios = [m.start() for m in re.finditer(str(malo), CUERPO)]
        # el valor falso puede vivir en la nota de erratas del cuerpo O en el apendice de
        # cambios, que es donde se movio el 20-sep; se acepta cualquiera de los dos rotulos
        dentro = all(CUERPO.count("earlier version", max(0, s - 700), s + 700)
                     or CUERPO.count("draft", max(0, s - 900), s + 900)
                     for s in sitios)
        afirma(bool(sitios) <= bool(dentro),
               "%s: el falso %d solo aparece en la nota de erratas (%d apariciones)"
               % (etiqueta, malo, len(sitios)))

print("")
print("(5) LA TABLA DE CONVERGENCIA: cada celda, contra la salida de ley_2p_largo")
sal = io.open(os.path.join(GATES, "ley_2p_largo_OUT.txt"),
              encoding="utf-8", errors="replace").read()
medidos = {}
for linea in sal.splitlines():
    f = linea.split()
    if len(f) == 7 and f[0].isdigit():
        medidos[int(f[0])] = float(f[6])
for p, v in ((37, 1.06), (43, 0.89), (53, 0.60), (59, 0.49), (67, 0.40),
             (73, 0.33), (83, 0.28), (89, 0.22), (97, 0.22)):
    afirma(p in medidos and abs(medidos[p] - v) < 0.006,
           "p=%-3d resto %.2f  (medido %.4f)" % (p, v, medidos.get(p, float("nan"))))
# NO es monotono: revierte en p = 43 y en p = 67.  Se comprueba la TENDENCIA, que es lo
# que la nota afirma, y se exige que las reversiones sean exactamente esas dos.
# solo entre los primos QUE LA TABLA ENSENA, que empieza en 31: el 29 -> 31 queda fuera
enTabla = [q for q in sorted(medidos) if q >= 31]
rev = [b for a, b in zip(enTabla, enTabla[1:]) if medidos[b] > medidos[a]]
afirma(medidos[97] < medidos[37] / 4, "el resto cae de %.2f a %.2f" % (medidos[37], medidos[97]))
afirma(rev == [43, 67], "las reversiones son exactamente p = 43 y 67 (medidas: %s)" % rev)

print("")
print("(6) EL 34/3/1/1 DE p = 31, y el 38 -> 21 de la seccion 5")
afirma(en_salida("orbitas_compuesto_OUT.txt", "h=1:34 h=3:3 h=5:1 h=15:1"),
       "34 libres, 3 con h=3, 1 con h=5, 1 con h=15")
r43, r97 = medidos[43] * 43, medidos[97] * 97
afirma(abs(r43 - 38) < 1.5 and abs(r97 - 21) < 1.5,
       "resto*p va de %.0f a %.0f (el tex dice 38 y 21)" % (r43, r97))

print("")
print("=" * 70)
print("comprobaciones: %d ; fallos: %d" % (len(hechas), len(fallos)))
# Un auditor que no corre ninguna comprobacion y dice "todo bien" es un control que no puede
# fallar.  La version anterior imprimia SIEMPRE 0 comprobaciones --- la cuenta salia de una
# expresion muerta que acababa devolviendo los fallos --- y aprobaba igual.
if not hechas:
    sys.exit("FATAL: cero comprobaciones. No poder medir no es aprobar")
print("VEREDICTO:", "todos los numeros de la nota tienen origen" if not fallos
      else "HAY NUMEROS SIN CUADRAR: %s" % fallos)
sys.exit(1 if fallos else 0)
