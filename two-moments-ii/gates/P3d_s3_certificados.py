# -*- coding: utf-8 -*-
u"""P3d_s3_certificados.py -- prop:s3: los certificados de irreducibilidad de P_k y la tabla que los imprime.

La desigualdad A < 2^d cubre todo k salvo un conjunto finito, y esos k se resuelven con
certificados exactos de irreducibilidad para P_k.  Este guion hace ese computo finito
(6 <= k <= 46), mide A frente a 2^d hasta k = 160, y escribe la tabla de certificados de la
nota.  Todo es exacto (sympy: recursion entera, Sturm, factorizacion mod p); los controles
son la reproduccion de s_3, s_4, s_5, s_7 por la recursion y el conjunto excepcional, que
tiene que salir exactamente {6,9,10,15,21,33,35} con igualdad en {6,10}.

EL ARGUMENTO, ENTERO, PARA PODER MEDIR CADA PASO POR SEPARADO.

  Sobre la fibra u_0 (raiz de Q_k) vive la cubica  Y^3 + u_0 Y - u_0.  Su discriminante es
  -u_0^2 (4u_0 + 27).  Si u_0 > 0 en todo encaje real --- que es lo que da la realidad total de
  Q_k --- ese discriminante es negativo en todo encaje, luego NO es un cuadrado, luego el grupo
  de Galois no es ciclico: o S_3 o la cubica es reducible sobre Q(u_0).

  Reducible sobre Q(u_0) significa que la cubica tiene raiz Y en Q(u_0).  De Y^3 = u_0(1-Y) sale
  u_0 = Y^3/(1-Y), asi que las raices de todas las fibras son las raices de

      P_k(Y) = (1-Y)^d Q_k(Y^3/(1-Y)),     d = deg Q_k,  deg P_k = 3d.

  Si TODA fibra fuese reducible, P_k tendria un factor racional de grado d (una raiz por fibra).
  Luego:  P_k irreducible sobre Q  ==>  alguna fibra es S_3.   Eso es prop:s3.

  (Y de hecho da mas: si P_k es irreducible de grado 3d, entonces Q_k es irreducible de grado d
  y TODA fibra es S_3, porque [Q(Y):Q] = 3d se factoriza por Q(u_0) obligando [Q(Y):Q(u_0)] = 3.)

LO QUE SE MIDE AQUI, Y EL ORDEN IMPORTA:

  1. Q_k, exacta, por la recursion  s_k = -e_2 s_{k-2} + e_3 s_{k-3}  con e_1 = 0.
  2. Todas las raices de Q_k reales, positivas y simples  (Sturm, no numerico).
  3. Termino constante de Q_k no nulo  --- lo necesita el paso |h(0)h(1)| >= 1.
  4. El discriminante de la fibra es -u^2(4u+27), identidad simbolica, negativo para u > 0.
  5. P_k tiene EXACTAMENTE d raices reales, todas en (0,1), y 2d complejas.
  6. P_k irreducible sobre Q, con CERTIFICADO: un primo p donde P_k mod p sigue siendo
     irreducible.  Eso es verificable en un renglon por quien quiera, sin confiar en nuestro
     factorizador.
  7. A = lc(Q_k) divide a k  --- DEMOSTRADA en P3d_lc_Qk.py (Waring en
     el monomio extremo: el coeficiente crudo es 2 si k es par y k si es impar, y primitivar
     solo puede dividirlo).  Aqui se sigue MIDIENDO, como control independiente de la prueba.
  8. A < 2^d falla exactamente en un conjunto finito, y se imprime cual es.

OBSERVACION.  P_k NO tiene todas sus raices en (0,1): tiene grado 3d y solo d raices reales.
La cubica de cada fibra tiene discriminante negativo --- el paso 4 --- luego UNA raiz real y
dos complejas.  Lo cierto, y lo que el argumento necesita, es "sus d raices reales estan todas
en (0,1)": el factor hipotetico h es el de las raices reales.  El punto 5 mide exactamente eso.

SALIDA.  Ademas de P3d_s3_certificados_OUT.txt, escribe junto a si mismo la tabla de la nota,
tabla_certificados.tex, y su pie, tabla_certificados_pie.tex (solo si su contenido cambia).

Uso:  python P3d_s3_certificados.py

Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
"""
import sys

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass

import io
import os

import sympy as sp

e2, e3, u, Y = sp.symbols('e2 e3 u Y')

# hasta donde se certifica irreducibilidad (caro) y hasta donde se mide A vs 2^d (barato)
K_CERT = 46
K_MIDE = 160

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


# --------------------------------------------------------------------------
# 1. la familia Q_k, exacta
# --------------------------------------------------------------------------
S = {0: sp.Integer(3), 1: sp.Integer(0), 2: -2 * e2}
for k in range(3, K_MIDE + 1):
    S[k] = sp.expand(-e2 * S[k - 2] + e3 * S[k - 3])


def Qk(k):
    u"""El factor propio de s_k escrito en u = e_2^3/e_3^2, en forma entera primitiva."""
    P = sp.Poly(S[k], e2, e3)
    a0 = min(m[0] for m in P.monoms())
    b0 = min(m[1] for m in P.monoms())
    resto = sp.Poly(sp.expand(sp.cancel(S[k] / (e2 ** a0 * e3 ** b0))), e2, e3)
    terms = []
    for (a, b), c in zip(resto.monoms(), resto.coeffs()):
        if a % 3 != 0:
            raise AssertionError(u'k=%d: el exponente de e_2 no es multiplo de 3' % k)
        terms.append((a // 3, c))
    q = sp.Poly(sum(c * u ** i for i, c in terms), u)
    _, prim = q.primitive()
    if prim.LC() < 0:
        prim = -prim
    return prim


def Pdek(q, d):
    u"""P_k(Y) = (1-Y)^d Q_k(Y^3/(1-Y)), construido coeficiente a coeficiente.

    Sustituir la fraccion y cancelar despues deja a sympy adivinando generadores --- y no
    adivina.  Aqui no hay nada que adivinar: si Q_k = sum c_i u^i entonces
    P_k = sum c_i Y^(3i) (1-Y)^(d-i), que es un polinomio de entrada.
    """
    c = q.all_coeffs()[::-1]
    return sp.Poly(sp.expand(sum(c[i] * Y ** (3 * i) * (1 - Y) ** (d - i)
                                 for i in range(d + 1))), Y)


# control del control: la recursion tiene que reproducir los s_k que ya estan en la nota
comprueba(u'la recursion reproduce s_3 = 3e_3', sp.expand(S[3] - 3 * e3) == 0)
comprueba(u'la recursion reproduce s_4 = 2e_2^2', sp.expand(S[4] - 2 * e2 ** 2) == 0)
comprueba(u'la recursion reproduce s_5 = -5e_2e_3', sp.expand(S[5] + 5 * e2 * e3) == 0)
comprueba(u'la recursion reproduce s_7 = 7e_2^2e_3', sp.expand(S[7] - 7 * e2 ** 2 * e3) == 0)

Q = {}
for k in range(6, K_MIDE + 1):
    Q[k] = Qk(k)

comprueba(u'deg Q_k = 0 exactamente en k de {3,4,5,7}',
      [k for k in range(3, K_MIDE + 1) if Qk(k).degree() == 0] == [3, 4, 5, 7])

# la formula del grado, que es prop:grado, comprobada aqui tambien porque este argumento la usa
malgrado = [k for k in range(6, K_MIDE + 1)
            if Q[k].degree() != k // 6 - (1 if k % 6 == 1 else 0)]
comprueba(u'deg Q_k = floor(k/6) - [k = 1 mod 6] en todo el rango',
      not malgrado, u'fallan: %s' % malgrado[:8])

# --------------------------------------------------------------------------
# 2. el discriminante de la fibra, como identidad simbolica
# --------------------------------------------------------------------------
disc = sp.factor(sp.discriminant(sp.Poly(Y ** 3 + u * Y - u, Y), Y))
comprueba(u'disc(Y^3 + uY - u) = -u^2 (4u + 27)',
      sp.expand(disc + u ** 2 * (4 * u + 27)) == 0, u'sympy da: %s' % disc)
comprueba(u'ese discriminante es < 0 para todo u > 0',
      sp.ask(sp.Q.negative(-u ** 2 * (4 * u + 27)), sp.Q.positive(u)) is True)

# --------------------------------------------------------------------------
# 3. realidad total, positividad, simplicidad, y el termino constante
# --------------------------------------------------------------------------
L.append(u'')
L.append(u'LAS RAICES DE Q_k Y DE P_k')
L.append(u'-' * 78)
L.append(u'%4s %3s %8s %7s   %-24s %-22s' % (
    u'k', u'd', u'A=lc', u'A | k', u'Q_k: reales>0 simples', u'P_k: d reales en (0,1)'))

fallo_raices = []
fallo_const = []
fallo_P01 = []
for k in range(6, K_CERT + 1):
    q = Q[k]
    d = q.degree()
    if d == 0:
        continue
    npos = sp.polys.polytools.count_roots(q, 0, sp.oo)
    simples = (sp.gcd(q, q.diff(u)).degree() == 0)
    if not (npos == d and simples):
        fallo_raices.append(k)
    if q.eval(0) == 0:
        fallo_const.append(k)
    Pk = Pdek(q, d)
    n01 = sp.polys.polytools.count_roots(Pk, 0, 1)
    nreal = sp.polys.polytools.count_roots(Pk)
    if not (Pk.degree() == 3 * d and n01 == d and nreal == d):
        fallo_P01.append((k, Pk.degree(), n01, nreal))
    A = abs(int(q.LC()))
    L.append(u'%4d %3d %8d %7s   %-24s %-22s' % (
        k, d, A, u'si' if k % A == 0 else u'NO',
        u'%d de %d, %s' % (npos, d, u'simples' if simples else u'REPETIDA'),
        u'%d de %d reales, %d en (0,1)' % (nreal, 3 * d, n01)))

comprueba(u'Q_k: las d raices son reales, positivas y simples',
      not fallo_raices, u'fallan: %s' % fallo_raices)
comprueba(u'Q_k(0) != 0  (lo necesita |h(0)h(1)| >= 1)',
      not fallo_const, u'fallan: %s' % fallo_const)
comprueba(u'P_k tiene grado 3d y EXACTAMENTE d raices reales, todas en (0,1)',
      not fallo_P01, u'fallan: %s' % fallo_P01[:5])

# --------------------------------------------------------------------------
# 4. los certificados de irreducibilidad
# --------------------------------------------------------------------------
L.append(u'')
L.append(u'LOS CERTIFICADOS:  un primo donde P_k mod p SIGUE siendo irreducible')
L.append(u'-' * 78)
L.append(u'   (irreducible modulo un primo ==> irreducible sobre Q; se comprueba en un renglon)')
L.append(u'')

sin_cert = []
reducibles = []
certs = []
for k in range(6, K_CERT + 1):
    q = Q[k]
    d = q.degree()
    if d == 0:
        continue
    Pk = Pdek(q, d)
    if not Pk.is_irreducible:
        reducibles.append(k)
        continue
    cert = None
    for p in sp.primerange(2, 600):
        if int(Pk.LC()) % p == 0:
            continue
        Fp = sp.Poly(Pk, Y, modulus=p)
        if sp.gcd(Fp, Fp.diff(Y)).degree() != 0:
            continue
        if Fp.is_irreducible:
            cert = p
            break
    if cert is None:
        sin_cert.append(k)
    certs.append((k, 3 * d, cert))

for i in range(0, len(certs), 4):
    L.append(u'   ' + u'   '.join(
        u'k=%-3d deg %-3d p=%-4s' % (a, b, c if c else u'--') for a, b, c in certs[i:i + 4]))

comprueba(u'P_k es irreducible sobre Q para todo 6 <= k <= %d con d > 0' % K_CERT,
      not reducibles, u'reducibles: %s' % reducibles)
comprueba(u'y cada uno lleva certificado mod p (irreducible ya modulo un primo)',
      not sin_cert, u'sin certificado: %s' % sin_cert)
comprueba(u'son %d casos, no una muestra' % len(certs), len(certs) >= 38,
      u'k con d > 0 en 6..%d' % K_CERT)

# --------------------------------------------------------------------------
# 5. la desigualdad, y el conjunto finito que se le escapa
# --------------------------------------------------------------------------
L.append(u'')
L.append(u'LA DESIGUALDAD  A < 2^d,  Y LO QUE SE LE ESCAPA')
L.append(u'-' * 78)

escapan = []
nodivide = []
for k in range(6, K_MIDE + 1):
    q = Q[k]
    d = q.degree()
    if d == 0:
        continue
    A = abs(int(q.LC()))
    if k % A != 0:
        nodivide.append(k)
    if not A < 2 ** d:
        escapan.append(k)

L.append(u'   se le escapan, en 6 <= k <= %d:  %s' % (K_MIDE, escapan))
L.append(u'   todos por debajo de %d, luego dentro del rango certificado'
         % (K_CERT + 1))

comprueba(u'A = lc(Q_k) divide a k en todo 6 <= k <= %d  (DEMOSTRADA en P3d_lc_Qk.py;'
      u' aqui se mide como control independiente)' % K_MIDE,
      not nodivide, u'no divide: %s' % nodivide[:8])
comprueba(u'los k que se le escapan a A < 2^d son exactamente {6,9,10,15,21,33,35}',
      escapan == [6, 9, 10, 15, 21, 33, 35], u'salen: %s' % escapan)
comprueba(u'y todos ellos estan dentro del rango certificado k <= %d' % K_CERT,
      escapan and max(escapan) <= K_CERT)

# el caso de igualdad, que es el que se cuela por un "<" escrito donde iba un "<="
igualdad = [k for k in range(6, K_MIDE + 1)
            if Q[k].degree() > 0 and abs(int(Q[k].LC())) == 2 ** Q[k].degree()]
L.append(u'   con IGUALDAD A = 2^d (el filo del argumento): %s' % igualdad)
comprueba(u'los casos de igualdad A = 2^d estan identificados', igualdad == [6, 10],
      u'salen: %s' % igualdad)

# --------------------------------------------------------------------------
L.append(u'')
L.append(u'VEREDICTO.  El computo finito de prop:s3 esta hecho: P_k es')
L.append(u'irreducible sobre Q para todo 6 <= k <= %d con deg Q_k > 0, cada caso con un primo' % K_CERT)
L.append(u'testigo.  Eso cubre con holgura los siete k que se le escapan a la desigualdad')
L.append(u'(%s).  Para k > %d la desigualdad se aplica' % (escapan, K_CERT))
L.append(u'en todo el rango, apoyada en que A divide a k --- DEMOSTRADA en P3d_lc_Qk.py, no')
L.append(u'medida: el coeficiente principal crudo de Q_k es 2 si k es par y k si es impar,')
L.append(u'por Waring en el monomio b = k mod 2, y primitivar solo puede dividirlo.  Con la')
L.append(u'cota floja A <= k el conjunto excepcional crece de estos 7 a 29 valores, todos')
L.append(u'k <= 37, y siguen cayendo dentro del rango certificado.  NO queda hueco.')

# --------------------------------------------------------------------------
# LA TABLA DE LA NOTA, ESCRITA POR LA MEDIDA
# --------------------------------------------------------------------------
# Un certificado que no se imprime no es un certificado: es una afirmacion.  La tabla la
# escribe este guion y la nota la \input --- no se teclea a mano, que es como los numeros de
# una tabla se separan de lo que miden.
def guarda(nombre, contenido):
    u"""Escribe, junto a este guion, SOLO si el contenido cambia.

    Reescribir la tabla en cada pasada adelantaria su fecha sin que nada haya cambiado, y
    cualquier control de frescura ("el PDF es mas nuevo que sus fuentes") suspenderia.
    """
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), nombre)
    if os.path.exists(d) and io.open(d, encoding='utf-8').read() == contenido:
        return False
    io.open(d, 'w', encoding='utf-8', newline='').write(contenido)
    return True


TAB = []
TAB.append(u'% GENERADO POR gates/P3d_s3_certificados.py --- NO EDITAR A MANO.')
TAB.append(u'% Cada fila es una afirmacion comprobable en un renglon: P_k mod p irreducible.')
TAB.append(u'\\begin{tabular}{rrrr@{\\qquad}rrrr@{\\qquad}rrrr}')
TAB.append(u'\\toprule')
cab = u'\\cellcolor{azul!14}$k$ & \\cellcolor{azul!14}$d$ & \\cellcolor{azul!14}$A$'\
      u' & \\cellcolor{azul!14}$p$'
TAB.append(u' & '.join([cab] * 3) + u'\\\\')
TAB.append(u'\\midrule')

# EL COLOR DICE DONDE HACE FALTA EL CERTIFICADO, que es lo unico que esta tabla tiene que
# comunicar de un vistazo.  Dos codigos, y ninguno decorativo:
#
#   * la celda de A lleva CALOR segun el margen 2^d / A: oscuro donde la desigualdad va justa,
#     palido donde sobra.  Se ve que el margen crece a lo bruto con k --- exponencial contra
#     lineal --- y que por eso el conjunto de excepciones es finito;
#   * las SIETE filas donde la desigualdad NO se aplica van en filacaliente.  Son exactamente
#     las que el certificado salva, y estan todas arriba.
#
# El lector ve en la banda naranja la razon de que el computo finito cierre el argumento en vez
# de aplazarlo.  Un color que no se pueda leer asi no se pone.
def calor(A, d):
    u"""Del margen 2^d/A a un tono de azul: 1 (va justo) -> oscuro, >= 64 -> casi blanco."""
    m = (2.0 ** d) / A
    if m <= 1:
        return None                     # no se aplica: la fila va en filacaliente
    for umbral, tono in ((2, 40), (4, 30), (8, 22), (16, 15), (64, 9)):
        if m < umbral:
            return tono
    return 5


datos = []
for k, g3, cert in certs:
    d = g3 // 3
    A = abs(int(Q[k].LC()))
    tono = calor(A, d)
    fila = u'\\cellcolor{filacaliente}' if tono is None else u''
    cA = u'' if tono is None else u'\\cellcolor{azul!%d}' % tono
    marca = u'$^\\dagger$' if tono is None else u''
    datos.append(u'%s%d%s & %s%d & %s%s%d & %s%s' % (
        fila, k, marca, fila, d, fila, cA, A, fila, cert if cert else u'--'))
alto = (len(datos) + 2) // 3
for i in range(alto):
    tres = [datos[i + j * alto] if i + j * alto < len(datos) else u' & & & '
            for j in range(3)]
    TAB.append(u' & '.join(tres) + u'\\\\')
TAB.append(u'\\bottomrule')
TAB.append(u'\\end{tabular}')
# La explicacion del dagger NO va aqui: quedaba entre la tabla y su pie, que se lee como si
# fueran dos pies.  Va dentro del \caption, y la compuerta la escribe aparte para que la nota
# la incluya en el sitio bueno sin teclear ningun numero.
# Y se entrega como MACRO, no como texto suelto: un \input dentro de \caption{} rompe --- el
# pie es un argumento movil y LaTeX da "Argument of \@tempf has an extra }".  El PDF llega a
# salir, asi que una comprobacion que solo cuente Overfull y Warning no lo ve; hay que contar
# las lineas del log que empiezan por "! ".
# El conjunto excepcional de la COTA DEMOSTRADA (A <= k), que es el que usa la prueba.  No se
# teclea: se cuenta aqui, igual que la tabla.  Ver gates/P3d_lc_Qk.py para la demostracion.
_grado = lambda kk: kk // 6 - (1 if kk % 6 == 1 else 0)
escapan_cota = [kk for kk in range(3, 5000) if _grado(kk) > 0 and kk >= 2 ** _grado(kk)]

PIE = (u'$\\dagger$: the %d values of $k$ for which the inequality $A<2^{d}$ does \\emph{not} '
       u'apply at the measured value of $A$; for $k=%d$ and $k=%d$ it fails by equality. The '
       u'proof uses only $A\\le k$, which is weaker and leaves %d such values rather than %d; '
       u'all of those lie inside the certified range too, which is why the finite computation '
       u'closes the argument rather than postponing it.'
       % (len(escapan), igualdad[0], igualdad[1], len(escapan_cota), len(escapan)))
# El % de la cabecera de LaTeX es tambien el operador de formato de Python: los dos literales
# se concatenan ANTES del %, asi que "% GENERADO" se lee como la especificacion "% G" (un float)
# y revienta con "must be real number, not str".  Se formatea aparte y se concatena despues.
guarda(u'tabla_certificados_pie.tex',
       u'% GENERADO POR gates/P3d_s3_certificados.py --- NO EDITAR A MANO.\n'
       + (u'\\newcommand{\\pieCertificados}{%s}\n' % PIE))
guarda(u'tabla_certificados.tex', u'\n'.join(TAB) + u'\n')
L.append(u'')
L.append(u'   tabla en tabla_certificados.tex, junto a este guion (%d filas)' % len(datos))
L.append(u'')
L.append(u'=' * 70)
L.append(u'VERIFICACION prop:s3 --- CERTIFICADOS:  %d ok, %d MAL' % (ok, mal))
L.append(u'=' * 70)

print(u'\n'.join(L))
SAL = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   u'P3d_s3_certificados_OUT.txt')
io.open(SAL, 'w', encoding='utf-8', newline='').write(u'\n'.join(L) + u'\n')
sys.exit(1 if mal else 0)
